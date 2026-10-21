/* lowgraph3: exact one-round transition graph of all differences with at most K active
 * nibbles (K <= 3), as lowgraph.c, with the half distributions from halffast.h
 * (conditioning on the active nibbles instead of enumerating 2^32 values).
 *
 *   lowgraph3 -k K -m MIN_LOG2P [-R a,b,c,d] [-t THREADS] [-c CACHEDIR] [-r ROUNDS]
 *
 * Edges: exact one-round transitions with probability >= 2^MIN_LOG2P between differences
 * with 1..K active nibbles. Output: best self-loops; minimum mean cycle of the whole graph
 * and of the subgraph of differences with exactly K active nibbles; spectral radius of
 * both; and (M^r) from 0x800000 to 0x800000, 0x900000, r = 1..ROUNDS (default 27), a lower
 * bound on the r-round probability under independent round keys.
 * With -c, half distributions are stored in CACHEDIR in blocks of 1000 jobs (resumable). */
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "halffast.h"
#include "graphalg.h"

typedef struct { uint64_t img; uint64_t num; uint16_t mask; } ment; /* permuted half output */
typedef struct { ment *e; size_t n; } mhalf;

static int K = 2, THREADS = 4;
static double MINLOG = -16;
static uint64_t MINNUM;
static const char *CACHE = NULL;

static uint16_t nibmask(uint64_t d) {
    uint16_t m = 0;
    for (int j = 0; j < 16; j++) if ((d >> (4 * j)) & 15) m |= (uint16_t)(1u << j);
    return m;
}

static uint32_t *HIN; static size_t NHIN;
static int nact32(uint32_t v) { int c = 0; for (int j = 0; j < 8; j++) c += ((v >> (4 * j)) & 15) != 0; return c; }
static void gen_half_inputs(void) {
    size_t cap = 1;
    for (int c = 1; c <= K; c++) { size_t t = 1; for (int i = 0; i < c; i++) t *= 15; cap += t * 70; }
    HIN = malloc(cap * sizeof(uint32_t));
    NHIN = 0;
    for (int c = 0; c <= K; c++) /* ordered by active count */
        for (int mask = 0; mask < 256; mask++) {
            if (__builtin_popcount(mask) != c) continue;
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

static mhalf *HL, *HR;
static size_t NJOB, NBLK;
#define BLK 1000
static _Atomic size_t NEXTBLK, DONEBLK;

static void do_job(size_t j) {
    size_t i = j >> 1;
    int right = (int)(j & 1);
    uint32_t v = HIN[i];
    size_t n;
    hf_ent *h = half_fast(right ? 2 : 0, right ? 3 : 1, (uint16_t)(v >> 16), (uint16_t)v, MINNUM, &n);
    mhalf *o = right ? &HR[i] : &HL[i];
    o->e = malloc((n ? n : 1) * sizeof(ment));
    o->n = 0;
    for (size_t q = 0; q < n; q++) {
        uint64_t img = right ? dip_perm(h[q].out) : dip_perm((uint64_t)h[q].out << 32);
        uint16_t m = nibmask(img);
        if (__builtin_popcount(m) > K) continue;
        o->e[o->n].img = img; o->e[o->n].num = h[q].num; o->e[o->n].mask = m; o->n++;
    }
    free(h);
    o->e = realloc(o->e, (o->n ? o->n : 1) * sizeof(ment)); /* keep only the filtered entries */
}

static int load_blk(size_t b) {
    char path[512];
    snprintf(path, sizeof path, "%s/blk%05zu.bin", CACHE, b);
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    size_t lo = b * BLK, hi = lo + BLK < NJOB ? lo + BLK : NJOB, j;
    for (j = lo; j < hi; j++) {
        mhalf *o = (j & 1) ? &HR[j >> 1] : &HL[j >> 1];
        if (fread(&o->n, sizeof o->n, 1, f) != 1) break;
        o->e = malloc((o->n ? o->n : 1) * sizeof(ment));
        if (fread(o->e, sizeof(ment), o->n, f) != o->n) break;
    }
    fclose(f);
    return j == hi;
}
static void save_blk(size_t b) {
    char path[512], tmp[520];
    snprintf(path, sizeof path, "%s/blk%05zu.bin", CACHE, b);
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "wb");
    if (!f) return;
    size_t lo = b * BLK, hi = lo + BLK < NJOB ? lo + BLK : NJOB;
    int ok = 1;
    for (size_t j = lo; j < hi && ok; j++) {
        mhalf *o = (j & 1) ? &HR[j >> 1] : &HL[j >> 1];
        ok = fwrite(&o->n, sizeof o->n, 1, f) == 1 && fwrite(o->e, sizeof(ment), o->n, f) == o->n;
    }
    if (!fclose(f) && ok) rename(tmp, path);
}

static void *half_worker(void *arg) {
    (void)arg;
    for (;;) {
        size_t b = atomic_fetch_add(&NEXTBLK, 1);
        if (b >= NBLK) return NULL;
        if (!(CACHE && load_blk(b))) {
            size_t lo = b * BLK, hi = lo + BLK < NJOB ? lo + BLK : NJOB;
            for (size_t j = lo; j < hi; j++) do_job(j);
            if (CACHE) save_blk(b);
        }
        size_t d = atomic_fetch_add(&DONEBLK, 1) + 1;
        if (d % 10 == 0 || d == NBLK) { fprintf(stderr, "half blocks %zu / %zu\n", d, NBLK); fflush(stderr); }
    }
}

/* nodes */
static uint64_t *NODES; static size_t NN, HCAP; static int64_t *HT;
static uint32_t *NI, *NJ;
static int64_t node_get(uint64_t d) {
    size_t h = (size_t)((d * 0x9E3779B97F4A7C15ULL) >> 20) & (HCAP - 1);
    while (HT[h] >= 0) { if (NODES[HT[h]] == d) return HT[h]; h = (h + 1) & (HCAP - 1); }
    return -1;
}
static void node_add(uint64_t d, uint32_t i, uint32_t j) {
    size_t h = (size_t)((d * 0x9E3779B97F4A7C15ULL) >> 20) & (HCAP - 1);
    while (HT[h] >= 0) h = (h + 1) & (HCAP - 1);
    HT[h] = (int64_t)NN; NI[NN] = i; NJ[NN] = j; NODES[NN++] = d;
}

/* edges, built in parallel over nodes */
typedef struct { g_edge *e; size_t n, cap; } ebuf;
static ebuf *EB;
static _Atomic size_t NEXTNODE;
static void *edge_worker(void *arg) {
    ebuf *B = arg;
    const double thr = pow(2.0, MINLOG);
    for (;;) {
        size_t x0 = atomic_fetch_add(&NEXTNODE, 4096);
        if (x0 >= NN) return NULL;
        size_t x1 = x0 + 4096 < NN ? x0 + 4096 : NN;
        for (size_t x = x0; x < x1; x++) {
            const mhalf *L = &HL[NI[x]], *R = &HR[NJ[x]];
            for (size_t a = 0; a < L->n; a++) {
                double pl = ldexp((double)L->e[a].num, -32);
                if (pl < thr) break;
                for (size_t b = 0; b < R->n; b++) {
                    double p = pl * ldexp((double)R->e[b].num, -32);
                    if (p < thr) break;
                    if (__builtin_popcount(L->e[a].mask | R->e[b].mask) > K) continue;
                    uint64_t o = L->e[a].img ^ R->e[b].img;
                    if (!o) continue;
                    int64_t t = node_get(o);
                    if (t < 0) continue;
                    if (B->n == B->cap) { B->cap = B->cap ? 2 * B->cap : 1 << 16; B->e = realloc(B->e, B->cap * sizeof(g_edge)); }
                    B->e[B->n].from = (int32_t)x; B->e[B->n].to = (int32_t)t; B->e[B->n].p = p; B->n++;
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    int rot[4] = {1, 4, 7, 11}, R = 27;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "-k")) K = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-m")) MINLOG = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "-R")) sscanf(argv[i + 1], "%d,%d,%d,%d", &rot[0], &rot[1], &rot[2], &rot[3]);
        else if (!strcmp(argv[i], "-t")) THREADS = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-c")) CACHE = argv[i + 1];
        else if (!strcmp(argv[i], "-r")) R = atoi(argv[i + 1]);
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    if (K < 1 || K > 3 || THREADS < 1 || THREADS > 64) { fprintf(stderr, "need 1 <= K <= 3, 1 <= THREADS <= 64\n"); return 1; }
    dip_init_rot(rot[0], rot[1], rot[2], rot[3]);
    MINNUM = (uint64_t)ldexp(1.0, 32 + (int)floor(MINLOG));
    if (MINNUM < 1) MINNUM = 1;
    gen_half_inputs();
    NJOB = 2 * NHIN; NBLK = (NJOB + BLK - 1) / BLK;
    printf("rotations=%d,%d,%d,%d K=%d min_log2p=%.2f half_inputs=%zu per half, %zu blocks\n", rot[0], rot[1], rot[2],
           rot[3], K, MINLOG, NHIN, NBLK);
    fflush(stdout);
    HL = calloc(NHIN, sizeof *HL); HR = calloc(NHIN, sizeof *HR);
    pthread_t th[64];
    for (int t = 0; t < THREADS; t++) pthread_create(&th[t], NULL, half_worker, NULL);
    for (int t = 0; t < THREADS; t++) pthread_join(th[t], NULL);
    size_t tot = 0;
    for (size_t i = 0; i < NHIN; i++) tot += HL[i].n + HR[i].n;
    printf("half entries kept (p >= 2^%.0f, <= %d active nibbles after P): %zu\n", MINLOG, K, tot);

    /* nodes: (L, R) with 1..K active nibbles in total */
    size_t cnt[4] = {0}, start[5] = {0};
    for (size_t i = 0; i < NHIN; i++) cnt[nact32(HIN[i])]++;
    for (int c = 0; c < 4; c++) start[c + 1] = start[c] + cnt[c];
    size_t maxn = 0;
    for (int a = 0; a <= K; a++) for (int b = 0; a + b <= K; b++) maxn += cnt[a] * cnt[b];
    maxn -= 1;
    HCAP = 1; while (HCAP < 2 * maxn + 16) HCAP <<= 1;
    NODES = malloc(maxn * sizeof *NODES); NI = malloc(maxn * sizeof *NI); NJ = malloc(maxn * sizeof *NJ);
    HT = malloc(HCAP * sizeof *HT);
    memset(HT, 0xff, HCAP * sizeof *HT);
    for (int a = 0; a <= K; a++)
        for (int b = 0; a + b <= K; b++)
            for (size_t i = start[a]; i < start[a + 1]; i++)
                for (size_t j = start[b]; j < start[b + 1]; j++) {
                    uint64_t d = ((uint64_t)HIN[i] << 32) | HIN[j];
                    if (d) node_add(d, (uint32_t)i, (uint32_t)j);
                }
    EB = calloc((size_t)THREADS, sizeof *EB);
    for (int t = 0; t < THREADS; t++) pthread_create(&th[t], NULL, edge_worker, &EB[t]);
    for (int t = 0; t < THREADS; t++) pthread_join(th[t], NULL);
    size_t ne = 0;
    for (int t = 0; t < THREADS; t++) ne += EB[t].n;
    g_edge *E = malloc((ne ? ne : 1) * sizeof *E);
    ne = 0;
    for (int t = 0; t < THREADS; t++) { memcpy(E + ne, EB[t].e, EB[t].n * sizeof *E); ne += EB[t].n; free(EB[t].e); }
    printf("nodes=%zu edges=%zu\n", NN, ne);
    fflush(stdout);

    /* best self-loops, overall and among nodes with exactly K active nibbles */
    for (int pass = 0; pass < 2; pass++) {
        printf("best self-loops (log2 p)%s:\n", pass ? " with exactly K active nibbles" : "");
        double last = 2;
        for (int it = 0; it < 6; it++) {
            double best = 0;
            for (size_t e = 0; e < ne; e++)
                if (E[e].from == E[e].to && E[e].p > best && E[e].p < last &&
                    (!pass || __builtin_popcount(nibmask(NODES[E[e].from])) == K)) best = E[e].p;
            if (best == 0) break;
            int same = 0;
            for (size_t e = 0; e < ne; e++)
                if (E[e].from == E[e].to && E[e].p == best && (!pass || __builtin_popcount(nibmask(NODES[E[e].from])) == K)) {
                    if (same < 4) printf("  %016" PRIx64 " %.4f\n", NODES[E[e].from], log2(best));
                    same++;
                }
            if (same > 4) printf("  ... %d loops with %.4f\n", same, log2(best));
            last = best;
        }
    }
    fflush(stdout);
    uint8_t *exK = malloc(NN);
    for (size_t i = 0; i < NN; i++) exK[i] = __builtin_popcount(nibmask(NODES[i])) == K;
    g_min_mean_cycle(NN, E, ne, NULL, NODES, "all");
    g_min_mean_cycle(NN, E, ne, exK, NODES, "exactly K");
    fflush(stdout);
    g_spectral(NN, E, ne, NULL, NODES, "all", 6);
    g_spectral(NN, E, ne, exK, NODES, "exactly K", 3);
    fflush(stdout);
    int64_t a = node_get(0x800000), b = node_get(0x900000);
    if (a >= 0 && b >= 0) {
        int32_t tg[2] = {(int32_t)a, (int32_t)b};
        printf("(M^r) from 0x800000 to 0x800000, 0x900000 (log2) and their sum:\n");
        g_propagate(NN, E, ne, (int32_t)a, tg, 2, R, "prop");
    }
    return 0;
}
