/* krec: key recovery on R-round Dipper (real key schedule, one fixed random master key)
 * with the iterative differential 0x800000 -> 0x800000 over rounds 1..R-2.
 *
 *   krec -r R -k ks96|ks128 -s SEED -n LOG2STRUCT [-t THREADS] [-T 8[,9..]] [-o FILE]
 *
 * Data: 2^LOG2STRUCT structures of 16 chosen plaintexts that take all values of nibble 5
 * (bits 20..23) and agree elsewhere; all 120 pairs of a structure are used.
 * Round R has no key after it, so (keyless round R)^-1 is applied to every ciphertext:
 * y = x_{R-1} ^ rk_R, z = P^-1(y). A right pair (Delta x_{R-2} in nibble 5) has Delta z = 0
 * outside bits 27..31; this 59-bit filter needs no key.
 * Key guess: the 32 bits of P^-1(rk_R) in the C2 and D words (krank.h); the score of a
 * guess is the number of survivors whose Delta x_{R-2}, recomputed under the guess,
 * equals 0x800000 (or 0x100000*v for v in -T, default 8).
 * Output: survivors, how many are right pairs (checked with the known key), the score
 * and rank of the correct 32-bit guess among all 2^32, and the score histogram. */
#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef BS_LANES
#define BS_LANES 4
#endif
#include "dipper.h"
#include "bitslice.h"
#include "bsinv.h"
#include "krank.h"

static int R = 14, THREADS = 4, LOG2S = 20, KWORDS = 6;
static uint64_t SEED = 1, RK[64];
static uint16_t TMASK = 1u << 8;
static V KS[64][64];
static const char *OUTF = NULL;

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

/* 64x64 bit-matrix transpose: on return a[j] bit i = old a[i] bit j */
static void tr64(uint64_t a[64]) {
    uint64_t m = 0x00000000FFFFFFFFULL;
    for (int j = 32; j; j >>= 1, m ^= m << j)
        for (int k = 0; k < 64; k = ((k | j) + 1) & ~j) {
            uint64_t t = ((a[k] >> j) ^ a[k | j]) & m;
            a[k] ^= t << j;
            a[k | j] ^= t;
        }
}

/* every 16-bit group of the result is all ones or all zeros, from 4 random bits */
static inline uint64_t spread4(uint64_t x) {
    uint64_t o = 0;
    for (int g = 0; g < 4; g++) if ((x >> g) & 1) o |= 0xFFFFULL << (16 * g);
    return o;
}

static const uint64_t NIB5[4] = {0xAAAAAAAAAAAAAAAAULL, 0xCCCCCCCCCCCCCCCCULL, 0xF0F0F0F0F0F0F0F0ULL,
                                 0xFF00FF00FF00FF00ULL};

typedef struct { kr_pair *p; size_t n, cap; } pbuf;
static pbuf *PB;
static uint64_t NCHUNK, CHUNKB; /* chunks of CHUNKB batches; one batch = BS_BLOCKS/16 structures */
static _Atomic uint64_t NEXTC, DONEC;

static uint64_t chunk_seed(uint64_t k) {
    uint64_t x = SEED * 0x9E3779B97F4A7C15ULL ^ (k + 1) * 0xD1B54A32D192ED03ULL ^ 0x6b72656300000000ULL;
    return splitmix(&x);
}

static void *data_worker(void *arg) {
    pbuf *B = arg;
    V X[64], X2[64], P0[64], Z[64], M[16];
    for (;;) {
        uint64_t k = atomic_fetch_add(&NEXTC, 1);
        if (k >= NCHUNK) return NULL;
        vrng_t r;
        vrng_seed(&r, chunk_seed(k));
        for (uint64_t bt = 0; bt < CHUNKB; bt++) {
            for (int i = 0; i < 64; i++) {
                if (i >= 20 && i < 24) { X[i] = vbcast(NIB5[i - 20]); continue; }
                V v = vrng_next(&r);
                for (int l = 0; l < BS_LANES; l++) X[i][l] = spread4(v[l]);
            }
            memcpy(P0, X, sizeof X);
            int q = 0;
            for (; q + 1 < R; q += 2) { bsi_round(X, X2, KS[q]); bsi_round(X2, X, KS[q + 1]); }
            if (q < R) { bsi_round(X, X2, KS[q]); memcpy(X, X2, sizeof X); }
            bsi_unpeel(X, Z);
            if (!bsi_match(Z, M)) continue;
            for (int l = 0; l < BS_LANES; l++) {
                int hit = 0;
                for (int d = 1; d < 16; d++) hit |= M[d][l] != 0;
                if (!hit) continue;
                uint64_t z[64], pt[64];
                for (int i = 0; i < 64; i++) { z[i] = Z[i][l]; pt[i] = P0[i][l]; }
                tr64(z);
                tr64(pt);
                for (int d = 1; d < 16; d++)
                    for (uint64_t m = M[d][l]; m; m &= m - 1) {
                        int a = __builtin_ctzll(m), b = a ^ d;
                        if (b < a) continue;
                        if (B->n == B->cap) { B->cap = B->cap ? 2 * B->cap : 1024; B->p = realloc(B->p, B->cap * sizeof(kr_pair)); }
                        kr_pair *o = &B->p[B->n++];
                        o->p0 = pt[a]; o->p1 = pt[b]; o->z0 = z[a]; o->z1 = z[b];
                        o->truth = dip_encrypt_rk(o->p0, RK, R - 2) ^ dip_encrypt_rk(o->p1, RK, R - 2);
                    }
            }
        }
        uint64_t d = atomic_fetch_add(&DONEC, 1) + 1;
        if (d % 64 == 0 || d == NCHUNK) { fprintf(stderr, "chunks %" PRIu64 "/%" PRIu64 "\n", d, NCHUNK); fflush(stderr); }
    }
}

/* ranking */
static kr_pair *SP; static size_t NS;
static uint16_t *VT;
static uint16_t KC_TRUE, KD_TRUE, S_TRUE;
static _Atomic uint32_t NEXTKD;
typedef struct { uint64_t hist[65536]; uint64_t top_key[16]; uint16_t top_s[16]; } rres;

static void *rank_worker(void *arg) {
    rres *o = arg;
    uint16_t s[32];
    for (;;) {
        uint32_t kd = atomic_fetch_add(&NEXTKD, 1);
        if (kd >= 65536) return NULL;
        for (uint32_t lo = 0; lo < 2048; lo++) {
            kr_score32(SP, VT, NS, (uint16_t)lo, (uint16_t)kd, s);
            for (int h = 0; h < 32; h++) {
                o->hist[s[h]]++;
                if (s[h] > o->top_s[15]) { /* keep the 16 best keys */
                    int q = 15;
                    while (q > 0 && o->top_s[q - 1] < s[h]) { o->top_s[q] = o->top_s[q - 1]; o->top_key[q] = o->top_key[q - 1]; q--; }
                    o->top_s[q] = s[h]; o->top_key[q] = ((uint64_t)((h << 11) | lo) << 16) | kd;
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    for (int i = 1; i + 1 < argc; i += 2) {
        const char *o = argv[i], *v = argv[i + 1];
        if (!strcmp(o, "-r")) R = atoi(v);
        else if (!strcmp(o, "-k")) KWORDS = !strcmp(v, "ks128") ? 8 : 6;
        else if (!strcmp(o, "-s")) SEED = strtoull(v, 0, 10);
        else if (!strcmp(o, "-n")) LOG2S = atoi(v);
        else if (!strcmp(o, "-t")) THREADS = atoi(v);
        else if (!strcmp(o, "-o")) OUTF = v;
        else if (!strcmp(o, "-T")) {
            TMASK = 0;
            char *s = strdup(v), *tok = strtok(s, ",");
            while (tok) { TMASK |= (uint16_t)(1u << (strtoul(tok, 0, 16) & 15)); tok = strtok(0, ","); }
            free(s);
        } else { fprintf(stderr, "unknown option %s\n", o); return 1; }
    }
    if (R < 4 || R > 64 || THREADS < 1 || THREADS > 256 || LOG2S < 4 || LOG2S > 44) {
        fprintf(stderr, "usage: krec -r R -k ks96|ks128 -s SEED -n LOG2STRUCT [-t THREADS] [-T 8[,9]] [-o FILE]\n");
        return 1;
    }
    dip_init();
    bs_init();
    /* master key from the seed */
    dip_key K;
    memset(&K, 0, sizeof K);
    K.words = KWORDS;
    uint64_t ks = SEED ^ 0x4b45594b45594b45ULL;
    for (int w = 0; w < KWORDS; w++) K.k[w] = (uint16_t)splitmix(&ks);
    dip_expand(K, R, RK);
    for (int q = 0; q < R; q++) bs_key(RK[q], KS[q]);
    /* check: bitsliced encryption and transpose against the table implementation */
    {
        uint64_t blk[BS_BLOCKS], c[64];
        rng_t g; rng_seed(&g, 99);
        for (int b = 0; b < BS_BLOCKS; b++) blk[b] = rng_next(&g);
        V X[64];
        bs_load(blk, X);
        for (int q = 0; q < R; q++) bs_round(X, KS[q]);
        int bad = 0;
        for (int l = 0; l < BS_LANES; l++) {
            for (int i = 0; i < 64; i++) c[i] = X[i][l];
            tr64(c);
            for (int j = 0; j < 64; j++) bad += c[j] != dip_encrypt_rk(blk[64 * l + j], RK, R);
        }
        /* fast round and unpeel against the table implementation */
        V Y[64], Y2[64], Zs[64];
        bs_load(blk, Y);
        for (int q = 0; q < R; q++) { bsi_round(Y, Y2, KS[q]); memcpy(Y, Y2, sizeof Y); }
        bsi_unpeel(Y, Zs);
        for (int l = 0; l < BS_LANES; l++) {
            for (int i = 0; i < 64; i++) c[i] = Zs[i][l];
            tr64(c);
            for (int j = 0; j < 64; j++)
                bad += c[j] != dip_perm_inv(dip_round_nokey_inv(dip_encrypt_rk(blk[64 * l + j], RK, R)));
        }
        if (bad) { fprintf(stderr, "selfcheck FAILED: %d blocks\n", bad); return 1; }
    }
    uint64_t kz = dip_perm_inv(RK[R - 1]);
    KC_TRUE = (uint16_t)(kz >> 16); KD_TRUE = (uint16_t)kz;
    printf("cmd=krec rounds=%d key=%s seed=%" PRIu64 " log2structures=%d targets=0x%04x\n", R,
           KWORDS == 8 ? "ks128" : "ks96", SEED, LOG2S, TMASK);
    printf("master key (words, most significant first):");
    for (int w = KWORDS - 1; w >= 0; w--) printf(" %04x", K.k[w]);
    printf("\nrk_R=%016" PRIx64 " P^-1(rk_R)=%016" PRIx64 " true guess kc=%04x kd=%04x\n", RK[R - 1], kz, KC_TRUE,
           KD_TRUE);

    /* data */
    uint64_t nbatch = (1ULL << LOG2S) / (BS_BLOCKS / 16);
    CHUNKB = nbatch < 4096 ? nbatch : 4096;
    NCHUNK = nbatch / CHUNKB;
    PB = calloc((size_t)THREADS, sizeof *PB);
    pthread_t th[256];
    double t0 = now();
    for (int t = 0; t < THREADS; t++) pthread_create(&th[t], NULL, data_worker, &PB[t]);
    for (int t = 0; t < THREADS; t++) pthread_join(th[t], NULL);
    double tdata = now() - t0;
    NS = 0;
    for (int t = 0; t < THREADS; t++) NS += PB[t].n;
    SP = malloc((NS ? NS : 1) * sizeof *SP);
    NS = 0;
    for (int t = 0; t < THREADS; t++) { memcpy(SP + NS, PB[t].p, PB[t].n * sizeof *SP); NS += PB[t].n; }
    /* sort survivors by plaintext so that the output does not depend on thread timing */
    for (size_t i = 1; i < NS; i++) {
        kr_pair x = SP[i]; size_t j = i;
        while (j > 0 && (SP[j - 1].p0 > x.p0 || (SP[j - 1].p0 == x.p0 && SP[j - 1].p1 > x.p1))) { SP[j] = SP[j - 1]; j--; }
        SP[j] = x;
    }
    uint64_t right = 0, nib5 = 0, other = 0, byv[16] = {0};
    for (size_t i = 0; i < NS; i++) {
        uint64_t d = SP[i].truth;
        if (d == 0x800000) right++;
        if (d && !(d & ~0xF00000ULL)) { nib5++; byv[(d >> 20) & 15]++; } else other++;
    }
    double cps = ldexp(1.0, LOG2S + 4);
    printf("data: chosen plaintexts=2^%.3f pairs=2^%.3f seconds=%.1f\n", log2(cps), log2(ldexp(120.0, LOG2S)), tdata);
    printf("survivors=%zu right(Delta x_{R-2}=0x800000)=%" PRIu64 " nibble5=%" PRIu64 " other=%" PRIu64 "\n", NS, right,
           nib5, other);
    printf("survivors by Delta x_{R-2} nibble-5 value v (0x100000*v):");
    for (int v = 1; v < 16; v++) printf(" %x:%" PRIu64, v, byv[v]);
    printf("\n");
    if (OUTF) {
        FILE *f = fopen(OUTF, "w");
        for (size_t i = 0; i < NS; i++)
            fprintf(f, "%016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 "\n", SP[i].p0, SP[i].p1,
                    SP[i].z0, SP[i].z1, SP[i].truth);
        fclose(f);
    }
    if (!NS) { printf("no survivors\n"); return 0; }

    /* ranking of the 2^32 guesses */
    VT = malloc(NS * 1024 * sizeof(uint16_t));
    for (size_t i = 0; i < NS; i++) kr_table(&SP[i], TMASK, VT + i * 1024);
    S_TRUE = 0;
    for (size_t i = 0; i < NS; i++) S_TRUE += (uint16_t)kr_direct(&SP[i], KC_TRUE, KD_TRUE, TMASK);
    { /* table method = plain formula on the true key and 4096 random keys */
        rng_t g; rng_seed(&g, SEED + 7);
        int bad = 0;
        for (int it = 0; it <= 4096; it++) {
            uint16_t kc = it ? (uint16_t)rng_next(&g) : KC_TRUE, kd = it ? (uint16_t)rng_next(&g) : KD_TRUE, s[32];
            kr_score32(SP, VT, NS, kc & 0x7FF, kd, s);
            uint16_t d = 0;
            for (size_t i = 0; i < NS; i++) d += (uint16_t)kr_direct(&SP[i], kc, kd, TMASK);
            bad += s[kc >> 11] != d;
        }
        printf("ranking selfcheck (table vs formula, 4097 keys): %s\n", bad ? "FAILED" : "ok");
        if (bad) return 1;
    }
    rres *RR = calloc((size_t)THREADS, sizeof *RR);
    t0 = now();
    for (int t = 0; t < THREADS; t++) pthread_create(&th[t], NULL, rank_worker, &RR[t]);
    for (int t = 0; t < THREADS; t++) pthread_join(th[t], NULL);
    double trank = now() - t0;
    uint64_t hist[65536] = {0};
    for (int t = 0; t < THREADS; t++) for (int s = 0; s < 65536; s++) hist[s] += RR[t].hist[s];
    uint64_t above = 0, tie = hist[S_TRUE];
    for (int s = S_TRUE + 1; s < 65536; s++) above += hist[s];
    int smax = 0;
    for (int s = 0; s < 65536; s++) if (hist[s]) smax = s;
    printf("ranking: seconds=%.1f score of the correct guess=%u max score=%d\n", trank, S_TRUE, smax);
    printf("guesses with a higher score=%" PRIu64 " with the same score=%" PRIu64 " (including the correct one)\n", above, tie);
    printf("candidates kept (score >= correct)=2^%.3f, i.e. %.3f of the 32 bits determined\n", log2((double)(above + tie)),
           32 - log2((double)(above + tie)));
    printf("score histogram (score:count), top 12 scores:");
    for (int s = smax, c = 0; s >= 0 && c < 12; s--) if (hist[s]) { printf(" %d:%" PRIu64, s, hist[s]); c++; }
    printf("\n");
    return 0;
}
