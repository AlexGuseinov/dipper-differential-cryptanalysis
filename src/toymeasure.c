/* toymeasure: exact per-key differential counts of Dipper-32 over the full codebook.
 *
 *   toymeasure -s SEED -n NKEYS -d DMAX [-S FILE] [-t THREADS] [-o HITFILE]
 *   toymeasure indep -d DMAX -l LOG2PAIRS -s SEED [-t THREADS]
 *
 * Keys: NKEYS random 48-bit master keys (Dipper-32/48 key schedule) from SEED.
 * For each key all 2^32 plaintexts are encrypted in 2^28 structures of 16 (all values of
 * nibble 3, bits 12..15). Per round r = 1..DMAX the program counts
 *   c8[r]: pairs (x, x ^ 0x8000) (2^31 of them) with output difference 0x8000;
 * and for each d listed in FILE (lines "d dout dout ..."), the pairs of each structure whose
 * difference after d rounds is in that target set, per input nibble difference v = a ^ b.
 * With -o, every structure with at least one such pair is written as a binary record
 * (uint8 key, uint8 d, uint32 structure, uint16 all-v count, uint16 v=8 count).
 * Mode indep: independent uniform round keys for every pair (the Markov model), counts
 * Pr[0x8000 -> 0x8000] for r = 1..DMAX over 2^LOG2PAIRS random pairs. */
#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "toy32.h"

#define MAXD 32
#define MAXS 8
static int DMAX = 16, THREADS = 8, NSET = 0;
static int SD[MAXS], SN[MAXS];
static uint32_t SV[MAXS][256];
static uint8_t SBM[MAXS][4096]; /* hashed bitmap of each set */
static uint32_t RK[MAXD];

static inline uint32_t hsh(uint32_t x) { return (x * 0x9E3779B1u) >> 20; }
static inline int inset(int q, uint32_t x) {
    if (!SBM[q][hsh(x)]) return 0;
    for (int i = 0; i < SN[q]; i++) if (SV[q][i] == x) return 1;
    return 0;
}

typedef struct { uint8_t d; uint32_t s; uint16_t all, v8; } hrec;
typedef struct {
    uint64_t c8[MAXD + 1], hv[MAXS][16];
    hrec *h; size_t nh, ch;
} tres;
static _Atomic uint64_t NEXT;
#define CH (1u << 16)

static void *worker(void *arg) {
    tres *T = arg;
    int at[MAXD + 1];
    for (int r = 0; r <= MAXD; r++) at[r] = -1;
    for (int q = 0; q < NSET; q++) at[SD[q]] = q;
    for (;;) {
        uint64_t c = atomic_fetch_add(&NEXT, 1);
        if (c >= (1u << 28) / CH) return NULL;
        for (uint32_t s = (uint32_t)c * CH; s < (uint32_t)(c + 1) * CH; s++) {
            uint32_t base = ((s >> 12) << 16) | (s & 0xFFF), y[16];
            for (int j = 0; j < 16; j++) y[j] = base | ((uint32_t)j << 12);
            for (int r = 1; r <= DMAX; r++) {
                for (int j = 0; j < 16; j++) y[j] = toy_round(y[j], RK[r - 1]);
                for (int a = 0; a < 8; a++) T->c8[r] += (y[a] ^ y[a ^ 8]) == 0x8000;
                int q = at[r];
                if (q < 0) continue;
                unsigned all = 0, v8 = 0;
                for (int a = 0; a < 16; a++)
                    for (int b = a + 1; b < 16; b++)
                        if (inset(q, y[a] ^ y[b])) { T->hv[q][a ^ b]++; all++; v8 += (a ^ b) == 8; }
                if (all) {
                    if (T->nh == T->ch) { T->ch = T->ch ? 2 * T->ch : 4096; T->h = realloc(T->h, T->ch * sizeof(hrec)); }
                    T->h[T->nh].d = (uint8_t)r; T->h[T->nh].s = s; T->h[T->nh].all = (uint16_t)all; T->h[T->nh].v8 = (uint16_t)v8;
                    T->nh++;
                }
            }
        }
    }
}

/* independent round keys per pair */
static int LOG2P = 30;
static uint64_t SEED = 1;
static void *worker_indep(void *arg) {
    tres *T = arg;
    for (;;) {
        uint64_t c = atomic_fetch_add(&NEXT, 1);
        if (c >= (1ULL << LOG2P) / CH) return NULL;
        rng_t g;
        rng_seed(&g, SEED * 0x9E3779B97F4A7C15ULL ^ (c + 1) * 0xD1B54A32D192ED03ULL);
        for (uint32_t i = 0; i < CH; i += 2) {
            uint64_t z = rng_next(&g);
            uint32_t x = (uint32_t)z, y = x ^ 0x8000;
            for (int r = 1; r <= DMAX; r++) {
                uint32_t k = (uint32_t)rng_next(&g);
                x = toy_round(x, k); y = toy_round(y, k);
                T->c8[r] += (x ^ y) == 0x8000;
            }
            /* second pair from the upper half of z */
            x = (uint32_t)(z >> 32); y = x ^ 0x8000;
            for (int r = 1; r <= DMAX; r++) {
                uint32_t k = (uint32_t)rng_next(&g);
                x = toy_round(x, k); y = toy_round(y, k);
                T->c8[r] += (x ^ y) == 0x8000;
            }
        }
    }
}

int main(int argc, char **argv) {
    int indep = argc > 1 && !strcmp(argv[1], "indep");
    int NK = 1;
    const char *sfile = NULL, *ofile = NULL;
    for (int i = 1 + indep; i + 1 < argc; i += 2) {
        const char *o = argv[i], *v = argv[i + 1];
        if (!strcmp(o, "-s")) SEED = strtoull(v, 0, 10);
        else if (!strcmp(o, "-n")) NK = atoi(v);
        else if (!strcmp(o, "-d")) DMAX = atoi(v);
        else if (!strcmp(o, "-S")) sfile = v;
        else if (!strcmp(o, "-o")) ofile = v;
        else if (!strcmp(o, "-t")) THREADS = atoi(v);
        else if (!strcmp(o, "-l")) LOG2P = atoi(v);
        else { fprintf(stderr, "unknown option %s\n", o); return 1; }
    }
    if (DMAX < 1 || DMAX > MAXD || THREADS < 1 || THREADS > 256 || NK < 1 || NK > 255) { fprintf(stderr, "bad arguments\n"); return 1; }
    toy_init(1, 4, 7, 3, 1, 2);
    { /* table round = reference round; keyless inverse */
        rng_t g; rng_seed(&g, 77);
        int bad = 0;
        for (int i = 0; i < 1000000; i++) {
            uint32_t x = (uint32_t)rng_next(&g), k = (uint32_t)rng_next(&g);
            bad += toy_round(x, k) != toy_round_ref(x, k);
            bad += toy_round_nokey_inv(toy_round_nokey(x)) != x;
        }
        printf("selfcheck (10^6 random rounds against the reference, keyless inverse): %s\n", bad ? "FAILED" : "ok");
        if (bad) return 1;
    }
    pthread_t th[256];
    tres *T = calloc((size_t)THREADS, sizeof *T);
    if (indep) {
        atomic_store(&NEXT, 0);
        for (int t = 0; t < THREADS; t++) pthread_create(&th[t], NULL, worker_indep, &T[t]);
        for (int t = 0; t < THREADS; t++) pthread_join(th[t], NULL);
        printf("cmd=toymeasure indep seed=%" PRIu64 " pairs=2^%d\n", SEED, LOG2P);
        for (int r = 1; r <= DMAX; r++) {
            uint64_t h = 0;
            for (int t = 0; t < THREADS; t++) h += T[t].c8[r];
            printf("r=%d hits=%" PRIu64 " log2p=%.4f\n", r, h, h ? log2((double)h) - LOG2P : -INFINITY);
        }
        return 0;
    }
    if (sfile) {
        FILE *f = fopen(sfile, "r");
        if (!f) { perror(sfile); return 1; }
        char line[8192];
        while (NSET < MAXS && fgets(line, sizeof line, f)) {
            char *tok = strtok(line, " \t\n");
            if (!tok || *tok == '#') continue;
            SD[NSET] = atoi(tok); SN[NSET] = 0;
            while ((tok = strtok(NULL, " \t\n")) && SN[NSET] < 256) SV[NSET][SN[NSET]++] = (uint32_t)strtoul(tok, 0, 16);
            for (int i = 0; i < SN[NSET]; i++) SBM[NSET][hsh(SV[NSET][i])] = 1;
            if (SD[NSET] >= 1 && SD[NSET] <= DMAX) NSET++;
        }
        fclose(f);
    }
    FILE *of = ofile ? fopen(ofile, "wb") : NULL;
    printf("cmd=toymeasure seed=%" PRIu64 " keys=%d dmax=%d sets=%d\n", SEED, NK, DMAX, NSET);
    for (int q = 0; q < NSET; q++) printf("set %d: d=%d size=%d\n", q, SD[q], SN[q]);
    uint64_t ks = SEED ^ 0x544f593332ULL;
    for (int kk = 0; kk < NK; kk++) {
        toy_key K;
        uint64_t z = splitmix(&ks);
        for (int w = 0; w < 6; w++) K.k[w] = (uint8_t)(z >> (8 * w));
        toy_expand(K, DMAX, RK);
        for (int t = 0; t < THREADS; t++) { free(T[t].h); memset(&T[t], 0, sizeof T[t]); }
        atomic_store(&NEXT, 0);
        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        for (int t = 0; t < THREADS; t++) pthread_create(&th[t], NULL, worker, &T[t]);
        for (int t = 0; t < THREADS; t++) pthread_join(th[t], NULL);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        printf("key %d master=%02x%02x%02x%02x%02x%02x seconds=%.1f\n", kk, K.k[5], K.k[4], K.k[3], K.k[2], K.k[1], K.k[0],
               (double)(t1.tv_sec - t0.tv_sec) + 1e-9 * (double)(t1.tv_nsec - t0.tv_nsec));
        printf("key %d c8:", kk);
        for (int r = 1; r <= DMAX; r++) {
            uint64_t h = 0;
            for (int t = 0; t < THREADS; t++) h += T[t].c8[r];
            printf(" %" PRIu64, h);
        }
        printf("\n");
        for (int q = 0; q < NSET; q++) {
            uint64_t tot = 0, hv[16] = {0};
            for (int t = 0; t < THREADS; t++) for (int v = 0; v < 16; v++) hv[v] += T[t].hv[q][v];
            for (int v = 1; v < 16; v++) tot += hv[v];
            printf("key %d set %d d=%d hits=%" PRIu64 " by v:", kk, q, SD[q], tot);
            for (int v = 1; v < 16; v++) printf(" %" PRIu64, hv[v]);
            printf("\n");
        }
        if (of)
            for (int t = 0; t < THREADS; t++)
                for (size_t i = 0; i < T[t].nh; i++) {
                    uint8_t kb = (uint8_t)kk, db = T[t].h[i].d;
                    fwrite(&kb, 1, 1, of); fwrite(&db, 1, 1, of); fwrite(&T[t].h[i].s, 4, 1, of);
                    fwrite(&T[t].h[i].all, 2, 1, of); fwrite(&T[t].h[i].v8, 2, 1, of);
                }
        fflush(stdout);
    }
    if (of) fclose(of);
    return 0;
}
