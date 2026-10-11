/* onestep: exact one-round analysis of Dipper by exhaustive enumeration.
 *
 *   onestep round -i DIN [-R a,b,c,d] [-m MIN_LOG2P]
 *       Exact distribution of the one-round output difference for input difference DIN,
 *       over a uniformly random state (equivalently a uniform round key).  The two halves
 *       (A,B) and (C,D) are independent, so each half is enumerated over all 2^32 inputs
 *       and the results are combined.  Prints every output difference with probability
 *       >= 2^MIN_LOG2P as "dout=<hex> num=<n> log2den=64 log2p=<x>".
 *
 *   onestep cluster -i DIN [-m MIN_LOG2P] [-N MAXNODES] [-r ROUNDS] [-f MATRIXFILE]
 *       Discovers up to MAXNODES differences reachable from DIN through exact one-round
 *       transitions of probability >= 2^MIN_LOG2P (best-path order), and iterates the
 *       resulting sub-stochastic matrix.  Every printed value is a lower bound, under the
 *       independent-round-key (Markov) model, on the probability of the r-round differential.
 *
 *   onestep bitgraph [-m MIN_LOG2P]
 *       For every rotation amount r in 0..15 and every input bit i of an added word
 *       (A or C), the exact number of (x, y) in 2^16 x 2^16 such that a one-bit input
 *       difference at bit i of the word gives the one-bit output difference e_j at the
 *       addition output.  Uses a carry DP over the uniform operand y.
 *       Prints "r=<r> i=<i> j=<j> num=<n> log2den=32".
 */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dipper.h"

/* ---- half enumeration ---- */
typedef struct { uint16_t outB; uint64_t *cnt; } bucket_t;

typedef struct { uint32_t out; uint64_t num; } hent_t; /* out = outHi<<16 | outLo, num / 2^32 */

static int cmp_hent(const void *a, const void *b) {
    const hent_t *x = a, *y = b;
    return x->num < y->num ? 1 : x->num > y->num ? -1 : (x->out > y->out) - (x->out < y->out);
}

/* hi word = added word (index wa), lo word = passed word (index wb).
 * Returns all half output differences with count >= minnum (count out of 2^32). */
static uint32_t BORD[65536];
static uint16_t BKEY[65536];
static int cmp_bkey(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (BKEY[x] > BKEY[y]) - (BKEY[x] < BKEY[y]);
}
static hent_t *half_dist(int wa, int wb, uint16_t da, uint16_t db, size_t *nout, uint64_t minnum) {
    hent_t *res;
    if (da == 0 && db == 0) {
        res = malloc(sizeof *res);
        res[0].out = 0; res[0].num = 1ULL << 32;
        *nout = 1;
        return res;
    }
    for (uint32_t B = 0; B < 65536; B++) { BORD[B] = B; BKEY[B] = DIP_TW[wb][B] ^ DIP_TW[wb][B ^ db]; }
    qsort(BORD, 65536, sizeof BORD[0], cmp_bkey);
    uint64_t *cnt = calloc(65536, sizeof(uint64_t));
    size_t n = 0, cap = 1024;
    res = malloc(cap * sizeof *res);
    size_t g0 = 0;
    while (g0 < 65536) {
        uint16_t ob = BKEY[BORD[g0]];
        size_t g1 = g0;
        while (g1 < 65536 && BKEY[BORD[g1]] == ob) g1++;
        memset(cnt, 0, 65536 * sizeof(uint64_t));
        for (size_t g = g0; g < g1; g++) {
            uint32_t B = BORD[g];
            uint16_t v = DIP_TW[wb][B], v2 = DIP_TW[wb][B ^ db];
            for (uint32_t A = 0; A < 65536; A++) {
                uint16_t u = DIP_TW[wa][A], u2 = DIP_TW[wa][A ^ da];
                cnt[(uint16_t)((uint16_t)(u + v) ^ (uint16_t)(u2 + v2))]++;
            }
        }
        for (int o = 0; o < 65536; o++)
            if (cnt[o] && cnt[o] >= minnum) {
                if (n == cap) { cap *= 2; res = realloc(res, cap * sizeof *res); }
                res[n].out = ((uint32_t)o << 16) | ob; res[n].num = cnt[o]; n++;
            }
        g0 = g1;
    }
    free(cnt);
    qsort(res, n, sizeof *res, cmp_hent);
    *nout = n;
    return res;
}

static int cmd_round(uint64_t din, double minlog) {
    size_t nl, nr;
    uint64_t mn = (uint64_t)ldexp(1.0, 32 + (int)floor(minlog)); if (mn < 1) mn = 1;
    hent_t *L = half_dist(0, 1, din >> 48, (din >> 32) & 0xFFFF, &nl, mn);
    hent_t *Rr = half_dist(2, 3, (din >> 16) & 0xFFFF, din & 0xFFFF, &nr, mn);
    uint64_t total = 0;
    for (size_t a = 0; a < nl; a++) {
        for (size_t b = 0; b < nr; b++) {
            double lp = log2((double)L[a].num) + log2((double)Rr[b].num) - 64.0;
            if (lp < minlog) break; /* Rr sorted descending */
            uint64_t st = ((uint64_t)L[a].out << 32) | Rr[b].out;
            uint64_t dout = dip_perm(st);
            printf("dout=%016" PRIx64 " numL=%" PRIu64 " numR=%" PRIu64 " log2den=64 log2p=%.6f\n", dout,
                   L[a].num, Rr[b].num, lp);
            total++;
        }
    }
    printf("din=%016" PRIx64 " left_outputs=%zu right_outputs=%zu printed=%" PRIu64 "\n", din, nl, nr, total);
    free(L); free(Rr);
    return 0;
}

/* ---- single-bit graph by carry DP ---- */
/* number of y in [0,2^16) with (u+y) xor (u2+y) == g */
static uint64_t add_count(uint16_t u, uint16_t u2, uint16_t g) {
    uint64_t st[2][2] = {{1, 0}, {0, 0}}; /* st[c][c2] */
    for (int i = 0; i < 16; i++) {
        uint64_t ns[2][2] = {{0, 0}, {0, 0}};
        int ui = (u >> i) & 1, u2i = (u2 >> i) & 1, gi = (g >> i) & 1;
        for (int c = 0; c < 2; c++)
            for (int c2 = 0; c2 < 2; c2++) {
                if (!st[c][c2]) continue;
                for (int y = 0; y < 2; y++) {
                    int s = ui ^ y ^ c, s2 = u2i ^ y ^ c2;
                    if ((s ^ s2) != gi) continue;
                    int nc = (ui & y) | (ui & c) | (y & c), nc2 = (u2i & y) | (u2i & c2) | (y & c2);
                    ns[nc][nc2] += st[c][c2];
                }
            }
        memcpy(st, ns, sizeof st);
    }
    return st[0][0] + st[0][1] + st[1][0] + st[1][1];
}

static int cmd_bitgraph(void) {
    for (int r = 0; r < 16; r++) {
        uint16_t T[65536];
        for (uint32_t x = 0; x < 65536; x++) T[x] = rotl16(dip_sub16((uint16_t)x), r);
        for (int i = 0; i < 16; i++) {
            uint64_t cnt[16] = {0};
            for (uint32_t A = 0; A < 65536; A++) {
                uint16_t u = T[A], u2 = T[A ^ (1u << i)];
                for (int j = 0; j < 16; j++) cnt[j] += add_count(u, u2, (uint16_t)(1u << j));
            }
            for (int j = 0; j < 16; j++)
                if (cnt[j]) printf("r=%d i=%d j=%d num=%" PRIu64 " log2den=32\n", r, i, j, cnt[j]);
        }
    }
    return 0;
}


/* ---- cluster: Markov-chain lower bound over a discovered set of differences ---- */
#define HC_SIZE 4096
typedef struct { int used; int side; uint16_t da, db; hent_t *e; size_t n; } hcache_t;
static hcache_t HC[HC_SIZE];
static uint64_t HC_MIN;
static const hent_t *half_cached(int side, uint16_t da, uint16_t db, size_t *n) {
    uint32_t h = ((uint32_t)side * 0x9E3779B1u) ^ ((uint32_t)da * 0x85EBCA6Bu) ^ ((uint32_t)db * 0xC2B2AE35u);
    for (uint32_t k = 0; k < HC_SIZE; k++) {
        hcache_t *c = &HC[(h + k) % HC_SIZE];
        if (!c->used) {
            c->used = 1; c->side = side; c->da = da; c->db = db;
            c->e = half_dist(side ? 2 : 0, side ? 3 : 1, da, db, &c->n, HC_MIN);
            *n = c->n;
            return c->e;
        }
        if (c->side == side && c->da == da && c->db == db) { *n = c->n; return c->e; }
    }
    fprintf(stderr, "half cache full\n");
    exit(1);
}

typedef struct { uint64_t d; double best; int done; } node_t;
typedef struct { int from, to; double p; } edge_t;

static int cmd_cluster(uint64_t din, double minlog, int N, int R, const char *mfile) {
    HC_MIN = (uint64_t)ldexp(1.0, 32 + (int)floor(minlog));
    if (HC_MIN < 1) HC_MIN = 1;
    node_t *nd = calloc((size_t)N, sizeof *nd);
    int nn = 1;
    nd[0].d = din; nd[0].best = 0;
    size_t ne = 0, ecap = 1 << 16;
    edge_t *E = malloc(ecap * sizeof *E);
    for (;;) {
        /* pick the unprocessed node with the best path probability */
        int pick = -1;
        for (int i = 0; i < nn; i++)
            if (!nd[i].done && (pick < 0 || nd[i].best > nd[pick].best)) pick = i;
        if (pick < 0) break;
        nd[pick].done = 1;
        uint64_t d = nd[pick].d;
        size_t nl, nr;
        const hent_t *L = half_cached(0, d >> 48, (d >> 32) & 0xFFFF, &nl);
        const hent_t *Rr = half_cached(1, (d >> 16) & 0xFFFF, d & 0xFFFF, &nr);
        for (size_t a = 0; a < nl; a++)
            for (size_t b = 0; b < nr; b++) {
                double lp = log2((double)L[a].num) + log2((double)Rr[b].num) - 64.0;
                if (lp < minlog) break;
                uint64_t o = dip_perm(((uint64_t)L[a].out << 32) | Rr[b].out);
                int t = -1;
                for (int i = 0; i < nn; i++) if (nd[i].d == o) { t = i; break; }
                if (t < 0 && nn < N) { t = nn++; nd[t].d = o; nd[t].best = -1e9; }
                if (t < 0) continue;
                if (nd[pick].best + lp > nd[t].best) nd[t].best = nd[pick].best + lp;
                if (ne == ecap) { ecap *= 2; E = realloc(E, ecap * sizeof *E); }
                E[ne].from = pick; E[ne].to = t;
                E[ne].p = ldexp((double)L[a].num, -32) * ldexp((double)Rr[b].num, -32);
                ne++;
            }
        fprintf(stderr, "processed %016" PRIx64 " nodes=%d edges=%zu\n", d, nn, ne);
    }
    if (mfile) {
        FILE *f = fopen(mfile, "w");
        fprintf(f, "# nodes %d edges %zu minlog %.2f start %016" PRIx64 "\n", nn, ne, minlog, din);
        for (int i = 0; i < nn; i++) fprintf(f, "n %d %016" PRIx64 "\n", i, nd[i].d);
        for (size_t e = 0; e < ne; e++) fprintf(f, "e %d %d %.17g\n", E[e].from, E[e].to, E[e].p);
        fclose(f);
    }
    double *v = calloc((size_t)nn, sizeof *v), *w = calloc((size_t)nn, sizeof *w);
    v[0] = 1;
    for (int r = 1; r <= R; r++) {
        memset(w, 0, (size_t)nn * sizeof *w);
        for (size_t e = 0; e < ne; e++) w[E[e].to] += v[E[e].from] * E[e].p;
        memcpy(v, w, (size_t)nn * sizeof *v);
        double tot = 0, mx = 0; int im = 0;
        for (int i = 0; i < nn; i++) { tot += v[i]; if (v[i] > mx) { mx = v[i]; im = i; } }
        printf("r=%d to_start_log2=%.4f best_out=%016" PRIx64 " best_log2=%.4f cluster_total_log2=%.4f\n", r,
               v[0] > 0 ? log2(v[0]) : -INFINITY, nd[im].d, log2(mx), log2(tot));
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: onestep round|bitgraph [options]\n"); return 1; }
    int rot[4] = {1, 4, 7, 11};
    uint64_t din = 0;
    double minlog = -20;
    int N = 256, R = 28;
    const char *mfile = NULL;
    for (int i = 2; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "-i")) din = strtoull(argv[i + 1], 0, 16);
        else if (!strcmp(argv[i], "-R")) sscanf(argv[i + 1], "%d,%d,%d,%d", &rot[0], &rot[1], &rot[2], &rot[3]);
        else if (!strcmp(argv[i], "-m")) minlog = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "-N")) N = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-r")) R = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-f")) mfile = argv[i + 1];
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    dip_init_rot(rot[0], rot[1], rot[2], rot[3]);
    if (!strcmp(argv[1], "round")) return cmd_round(din, minlog);
    if (!strcmp(argv[1], "bitgraph")) return cmd_bitgraph();
    if (!strcmp(argv[1], "cluster")) return cmd_cluster(din, minlog, N, R, mfile);
    fprintf(stderr, "unknown command\n");
    return 1;
}
