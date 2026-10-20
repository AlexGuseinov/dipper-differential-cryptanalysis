/* lowgraph: exact one-round transition graph of all differences with at most K active
 * nibbles, and its best iterative trails and clusters.
 *
 *   lowgraph -k K -m MIN_LOG2P [-R a,b,c,d] [-t THREADS] [-f EDGEFILE] [-c CACHEDIR]
 *
 * With -c the half distributions are cached in CACHEDIR (one file each), so an interrupted
 * run resumes; the cache depends on K, MIN_LOG2P and the rotations (use one directory per setting).
 *
 * Every half input difference (A,B) or (C,D) with at most K active nibbles is enumerated
 * exhaustively over its 2^32 values (as `onestep round`, using the pair symmetry to halve
 * the work); half transitions of probability >= 2^MIN_LOG2P are kept.  For every 64-bit
 * difference with 1..K active nibbles the half transitions are combined, permuted, and
 * those that land again on a difference with at most K active nibbles and have
 * probability >= 2^MIN_LOG2P become edges.  All probabilities are exact one-round
 * probabilities over a uniform state (equivalently a uniform round key).
 *
 * Output:
 *   - the minimum mean weight cycle (Karp), i.e. the best iterative trail in the graph,
 *     with its bits per round and its differences;
 *   - the best self-loops;
 *   - the decay per round of the whole graph as a cluster: log2 of the spectral radius of
 *     the transition matrix restricted to the graph (power iteration).  Under the
 *     independent-round-key model, for any start a and end b in the graph,
 *     Pr[a -> b, r rounds] >= (M^r)_{ab}, and max_{a,b} (M^r)_{ab} decays like rho^r.
 */
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dipper.h"

typedef struct { uint32_t out; uint64_t num; } hent_t;
typedef struct { uint16_t da, db; hent_t *e; size_t n; } half_t;

static int K = 2, THREADS = 4;
static const char *CACHE = NULL; /* -c DIR: half distributions are stored there and reused (resumable) */
static double MINLOG = -16;
static uint64_t MINNUM;

static int cmp_hent(const void *a, const void *b) {
    const hent_t *x = a, *y = b;
    return x->num < y->num ? 1 : x->num > y->num ? -1 : (x->out > y->out) - (x->out < y->out);
}

static int nib_active(uint64_t d) {
    int c = 0;
    for (int j = 0; j < 16; j++) c += ((d >> (4 * j)) & 15) != 0;
    return c;
}

/* exact distribution of the half output (outHi<<16 | outLo) for input (da, db) of the half
 * whose added word is wa and passed word wb; counts out of 2^32; keeps counts >= MINNUM.
 * Uses the symmetry (x, x^d) <-> (x^d, x): only one of each unordered pair is visited. */
static hent_t *half_dist(int wa, int wb, uint16_t da, uint16_t db, size_t *nout) {
    if (!da && !db) {
        hent_t *r = malloc(sizeof *r);
        r->out = 0; r->num = 1ULL << 32; *nout = 1;
        return r;
    }
    const uint16_t *TA = DIP_TW[wa], *TB = DIP_TW[wb];
    uint64_t *cnt = malloc(65536 * sizeof(uint64_t));
    /* group B by its output difference */
    uint32_t *ord = malloc(65536 * sizeof(uint32_t));
    uint32_t start[65537] = {0};
    static __thread uint32_t key[65536];
    for (uint32_t B = 0; B < 65536; B++) { key[B] = (uint16_t)(TB[B] ^ TB[B ^ db]); start[key[B] + 1]++; }
    for (int i = 0; i < 65536; i++) start[i + 1] += start[i];
    uint32_t pos[65536];
    memcpy(pos, start, sizeof pos);
    for (uint32_t B = 0; B < 65536; B++) ord[pos[key[B]]++] = B;
    int hbB = db ? 31 - __builtin_clz(db) : -1, hbA = da ? 31 - __builtin_clz(da) : -1;
    size_t n = 0, cap = 256;
    hent_t *res = malloc(cap * sizeof *res);
    for (uint32_t ob = 0; ob < 65536; ob++) {
        if (start[ob] == start[ob + 1]) continue;
        memset(cnt, 0, 65536 * sizeof(uint64_t));
        for (uint32_t g = start[ob]; g < start[ob + 1]; g++) {
            uint32_t B = ord[g];
            if (hbB >= 0 && ((B >> hbB) & 1)) continue; /* visit B < B^db only */
            uint16_t v = TB[B], v2 = TB[B ^ db];
            if (hbB >= 0) {
                for (uint32_t A = 0; A < 65536; A++) {
                    uint16_t u = TA[A], u2 = TA[A ^ da];
                    cnt[(uint16_t)((uint16_t)(u + v) ^ (uint16_t)(u2 + v2))]++;
                }
            } else { /* db == 0: halve over A */
                uint32_t lo = (1u << hbA) - 1;
                for (uint32_t a = 0; a < 32768; a++) {
                    uint32_t A = ((a & ~lo) << 1) | (a & lo);
                    uint16_t u = TA[A], u2 = TA[A ^ da];
                    cnt[(uint16_t)((uint16_t)(u + v) ^ (uint16_t)(u2 + v2))]++;
                }
            }
        }
        for (uint32_t o = 0; o < 65536; o++) {
            uint64_t c = 2 * cnt[o];
            if (c && c >= MINNUM) {
                if (n == cap) { cap *= 2; res = realloc(res, cap * sizeof *res); }
                res[n].out = (o << 16) | ob; res[n].num = c; n++;
            }
        }
    }
    free(cnt); free(ord);
    qsort(res, n, sizeof *res, cmp_hent);
    *nout = n;
    return res;
}

/* all 16-bit pairs (da, db) with at most K active nibbles in the 32-bit half */
static uint32_t *HIN; static size_t NHIN;
static void gen_half_inputs(void) {
    NHIN = 0;
    HIN = malloc(sizeof(uint32_t) * 70000);
    /* enumerate by active nibble sets of size <= K over 8 nibbles */
    for (int mask = 0; mask < 256; mask++) {
        int c = __builtin_popcount(mask);
        if (c > K) continue;
        int pos[8], np = 0;
        for (int j = 0; j < 8; j++) if (mask >> j & 1) pos[np++] = j;
        int tot = 1;
        for (int i = 0; i < np; i++) tot *= 15;
        for (int t = 0; t < tot; t++) {
            uint32_t v = 0; int q = t;
            for (int i = 0; i < np; i++) { v |= (uint32_t)(q % 15 + 1) << (4 * pos[i]); q /= 15; }
            HIN[NHIN++] = v;
        }
    }
}

static half_t *HL, *HR; /* indexed like HIN */
static _Atomic size_t NEXTJOB;
static void *worker(void *arg) {
    (void)arg;
    for (;;) {
        size_t j = atomic_fetch_add(&NEXTJOB, 1);
        if (j >= 2 * NHIN) break;
        size_t i = j >> 1;
        uint32_t v = HIN[i];
        half_t *h = (j & 1) ? &HR[i] : &HL[i];
        h->da = v >> 16; h->db = v & 0xFFFF;
        char path[512];
        FILE *f = NULL;
        if (CACHE) {
            snprintf(path, sizeof path, "%s/%c%04x%04x.bin", CACHE, (j & 1) ? 'R' : 'L', h->da, h->db);
            f = fopen(path, "rb");
        }
        if (f && fread(&h->n, sizeof h->n, 1, f) == 1 && (h->e = malloc((h->n ? h->n : 1) * sizeof *h->e)) &&
            fread(h->e, sizeof *h->e, h->n, f) == h->n) {
            fclose(f);
            continue;
        }
        if (f) fclose(f);
        h->e = (j & 1) ? half_dist(2, 3, h->da, h->db, &h->n) : half_dist(0, 1, h->da, h->db, &h->n);
        if (CACHE) { /* write to a temporary name, then rename: an interrupted write is never read back */
            char tmp[520];
            snprintf(tmp, sizeof tmp, "%s.tmp", path);
            FILE *g = fopen(tmp, "wb");
            if (g && fwrite(&h->n, sizeof h->n, 1, g) == 1 && fwrite(h->e, sizeof *h->e, h->n, g) == h->n && !fclose(g))
                rename(tmp, path);
        }
        if (j % 500 == 0) fprintf(stderr, "half job %zu / %zu\n", j, 2 * NHIN);
    }
    return NULL;
}

/* hash map 64-bit difference -> node index */
static uint64_t *NODES; static size_t NN, HCAP; static int64_t *HT;
static int64_t node_get(uint64_t d) {
    size_t h = (size_t)((d * 0x9E3779B97F4A7C15ULL) >> 20) & (HCAP - 1);
    while (HT[h] >= 0) { if (NODES[HT[h]] == d) return HT[h]; h = (h + 1) & (HCAP - 1); }
    return -1;
}
static void node_add(uint64_t d) {
    size_t h = (size_t)((d * 0x9E3779B97F4A7C15ULL) >> 20) & (HCAP - 1);
    while (HT[h] >= 0) h = (h + 1) & (HCAP - 1);
    HT[h] = (int64_t)NN; NODES[NN++] = d;
}

typedef struct { int32_t from, to; double p; } edge_t;
static size_t NALIVE; /* nodes on or leading to cycles */

int main(int argc, char **argv) {
    int rot[4] = {1, 4, 7, 11};
    const char *efile = NULL;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "-k")) K = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-m")) MINLOG = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "-R")) sscanf(argv[i + 1], "%d,%d,%d,%d", &rot[0], &rot[1], &rot[2], &rot[3]);
        else if (!strcmp(argv[i], "-t")) THREADS = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-f")) efile = argv[i + 1];
        else if (!strcmp(argv[i], "-c")) CACHE = argv[i + 1];
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    dip_init_rot(rot[0], rot[1], rot[2], rot[3]);
    MINNUM = (uint64_t)ldexp(1.0, 32 + (int)floor(MINLOG));
    if (MINNUM < 1) MINNUM = 1;
    gen_half_inputs();
    printf("rotations=%d,%d,%d,%d K=%d min_log2p=%.2f half_inputs=%zu per half\n", rot[0], rot[1], rot[2], rot[3],
           K, MINLOG, NHIN);
    HL = calloc(NHIN, sizeof *HL); HR = calloc(NHIN, sizeof *HR);
    pthread_t th[64];
    for (int t = 0; t < THREADS; t++) pthread_create(&th[t], NULL, worker, NULL);
    for (int t = 0; t < THREADS; t++) pthread_join(th[t], NULL);
    /* nodes: all 64-bit differences (L, R) with 1..K active nibbles */
    size_t maxn = 0;
    for (size_t i = 0; i < NHIN; i++)
        for (size_t j = 0; j < NHIN; j++) {
            uint64_t d = ((uint64_t)HIN[i] << 32) | HIN[j];
            if (d && nib_active(d) <= K) maxn++;
        }
    HCAP = 1; while (HCAP < 2 * maxn + 16) HCAP <<= 1;
    NODES = malloc(maxn * sizeof *NODES); HT = malloc(HCAP * sizeof *HT);
    memset(HT, 0xff, HCAP * sizeof *HT);
    size_t *NI = malloc(maxn * sizeof *NI), *NJ = malloc(maxn * sizeof *NJ);
    for (size_t i = 0; i < NHIN; i++)
        for (size_t j = 0; j < NHIN; j++) {
            uint64_t d = ((uint64_t)HIN[i] << 32) | HIN[j];
            if (d && nib_active(d) <= K) { NI[NN] = i; NJ[NN] = j; node_add(d); }
        }
    size_t ne = 0, ecap = 1 << 20;
    edge_t *E = malloc(ecap * sizeof *E);
    for (size_t x = 0; x < NN; x++) {
        const half_t *L = &HL[NI[x]], *R = &HR[NJ[x]];
        for (size_t a = 0; a < L->n; a++) {
            double pl = ldexp((double)L->e[a].num, -32);
            for (size_t b = 0; b < R->n; b++) {
                double p = pl * ldexp((double)R->e[b].num, -32);
                if (log2(p) < MINLOG) break;
                uint64_t o = dip_perm(((uint64_t)L->e[a].out << 32) | R->e[b].out);
                if (nib_active(o) > K) continue;
                int64_t t = node_get(o);
                if (t < 0) continue;
                if (ne == ecap) { ecap *= 2; E = realloc(E, ecap * sizeof *E); }
                E[ne].from = (int32_t)x; E[ne].to = (int32_t)t; E[ne].p = p; ne++;
            }
        }
    }
    printf("nodes=%zu edges=%zu\n", NN, ne);
    if (efile) {
        FILE *f = fopen(efile, "w");
        for (size_t e = 0; e < ne; e++)
            fprintf(f, "%016" PRIx64 " %016" PRIx64 " %.17g\n", NODES[E[e].from], NODES[E[e].to], E[e].p);
        fclose(f);
    }
    /* best self-loops */
    printf("best self-loops (log2 p):\n");
    {
        int shown = 0;
        double lastbest = 1;
        for (int it = 0; it < 8; it++) {
            double best = 0; size_t be = 0;
            for (size_t e = 0; e < ne; e++)
                if (E[e].from == E[e].to && E[e].p > best && E[e].p < lastbest) { best = E[e].p; be = e; }
            if (best == 0) break;
            int cntsame = 0;
            for (size_t e = 0; e < ne; e++)
                if (E[e].from == E[e].to && E[e].p == best) {
                    if (cntsame < 4) printf("  %016" PRIx64 " %.4f\n", NODES[E[e].from], log2(best));
                    cntsame++;
                }
            if (cntsame > 4) printf("  ... %d loops with %.4f\n", cntsame, log2(best));
            lastbest = best; shown++;
            (void)be;
        }
        (void)shown;
    }
    /* minimum mean weight cycle (w = -log2 p) by Howard's policy iteration on the subgraph
     * left after repeatedly removing nodes without outgoing edges */
    {
        size_t n = NN;
        int *alive = malloc(n * sizeof(int)), *outdeg = calloc(n, sizeof(int));
        for (size_t i = 0; i < n; i++) alive[i] = 1;
        for (int changed = 1; changed;) {
            changed = 0;
            memset(outdeg, 0, n * sizeof(int));
            for (size_t e = 0; e < ne; e++) if (alive[E[e].to] && alive[E[e].from]) outdeg[E[e].from]++;
            for (size_t i = 0; i < n; i++) if (alive[i] && !outdeg[i]) { alive[i] = 0; changed = 1; }
        }
        size_t nalive = 0;
        for (size_t i = 0; i < n; i++) nalive += alive[i];
        printf("nodes on or leading to cycles: %zu\n", nalive);
        NALIVE = nalive;
        int32_t *pol = malloc(n * sizeof(int32_t));
        double *pw = malloc(n * sizeof(double)), *eta = malloc(n * sizeof(double)), *x = malloc(n * sizeof(double));
        for (size_t i = 0; i < n; i++) { pol[i] = -1; pw[i] = INFINITY; }
        for (size_t e = 0; e < ne; e++) {
            if (!alive[E[e].from] || !alive[E[e].to]) continue;
            double w = -log2(E[e].p);
            if (w < pw[E[e].from]) { pw[E[e].from] = w; pol[E[e].from] = E[e].to; }
        }
        int *mark = malloc(n * sizeof(int));
        int32_t *stk = malloc(n * sizeof(int32_t));
        const double EPS = 1e-9;
        for (int iter = 0; iter < 10000; iter++) {
            /* value determination on the functional graph pol */
            for (size_t i = 0; i < n; i++) mark[i] = 0; /* 0 new, 1 on stack, 2 done */
            for (size_t s0 = 0; s0 < n; s0++) {
                if (!alive[s0] || mark[s0]) continue;
                int sp = 0; int32_t u = (int32_t)s0;
                while (!mark[u]) { mark[u] = 1; stk[sp++] = u; u = pol[u]; }
                if (mark[u] == 1) { /* new cycle through u */
                    double sum = 0; int len = 0; int32_t c = u;
                    do { sum += pw[c]; len++; c = pol[c]; } while (c != u);
                    double m = sum / len;
                    x[u] = 0; eta[u] = m;
                    /* nodes of the cycle after u, backwards: x(c) = pw(c) - m + x(pol c) */
                    int32_t cyc[4096]; int cl = 0; c = pol[u];
                    while (c != u && cl < 4096) { cyc[cl++] = c; c = pol[c]; }
                    for (int q = cl - 1; q >= 0; q--) { eta[cyc[q]] = m; x[cyc[q]] = pw[cyc[q]] - m + x[pol[cyc[q]]]; mark[cyc[q]] = 2; }
                    mark[u] = 2;
                }
                /* stack nodes not yet done: tree part, resolve from the top down */
                while (sp > 0) {
                    int32_t t = stk[--sp];
                    if (mark[t] == 2) continue;
                    eta[t] = eta[pol[t]]; x[t] = pw[t] - eta[t] + x[pol[t]]; mark[t] = 2;
                }
            }
            /* policy improvement */
            int ch = 0;
            for (size_t e = 0; e < ne; e++) {
                int32_t a2 = E[e].from, b2 = E[e].to;
                if (!alive[a2] || !alive[b2]) continue;
                if (eta[b2] < eta[a2] - EPS) { pol[a2] = b2; pw[a2] = -log2(E[e].p); eta[a2] = eta[b2]; ch = 1; }
            }
            if (!ch)
                for (size_t e = 0; e < ne; e++) {
                    int32_t a2 = E[e].from, b2 = E[e].to;
                    if (!alive[a2] || !alive[b2] || fabs(eta[b2] - eta[a2]) > EPS) continue;
                    double w = -log2(E[e].p);
                    if (w - eta[a2] + x[b2] < x[a2] - EPS) { pol[a2] = b2; pw[a2] = w; x[a2] = w - eta[a2] + x[b2]; ch = 1; }
                }
            if (!ch) { printf("Howard converged after %d iterations\n", iter + 1); break; }
        }
        double best = INFINITY; int32_t bv = -1;
        for (size_t i = 0; i < n; i++) if (alive[i] && eta[i] < best) { best = eta[i]; bv = (int32_t)i; }
        if (bv >= 0) {
            /* walk to the cycle of bv */
            int32_t u = bv;
            for (size_t q = 0; q < n; q++) u = pol[u];
            printf("minimum mean cycle: %.4f bits per round, cycle:", best);
            int32_t c = u; int len = 0; double sum = 0;
            do { printf(" %016" PRIx64, NODES[c]); sum += pw[c]; len++; c = pol[c]; } while (c != u && len < 64);
            printf("  (length %d, product 2^-%.4f)\n", len, sum);
        }
    }
    /* spectral radius by power iteration (0 when the graph has no cycle) */
    if (!NALIVE) printf("no cycle in the graph: spectral radius 0\n");
    else {
        double *v = malloc(NN * sizeof(double)), *w = malloc(NN * sizeof(double));
        for (size_t i = 0; i < NN; i++) v[i] = 1.0;
        double lam = 0;
        for (int it = 0; it < 2000; it++) {
            memset(w, 0, NN * sizeof(double));
            for (size_t e = 0; e < ne; e++) w[E[e].to] += v[E[e].from] * E[e].p;
            double s = 0;
            for (size_t i = 0; i < NN; i++) s = w[i] > s ? w[i] : s;
            if (s == 0) break;
            for (size_t i = 0; i < NN; i++) v[i] = w[i] / s;
            lam = s;
        }
        printf("spectral radius: log2 rho = %.4f bits per round\n", log2(lam));
        printf("largest components of the Perron vector:\n");
        for (int q = 0; q < 6; q++) {
            size_t bi = 0;
            for (size_t i = 0; i < NN; i++) if (v[i] > v[bi]) bi = i;
            printf("  %016" PRIx64 " %.4f\n", NODES[bi], v[bi]);
            v[bi] = -1;
        }
    }
    return 0;
}
