/* Bitsliced Dipper: 256 blocks per batch (4 x 64-bit lanes per slice).
 * Slice i holds state bit i of every block.  Rotations and the bit
 * permutation are index maps; S-box and the two adders are boolean circuits. */
#ifndef BITSLICE_H
#define BITSLICE_H

#include <stdint.h>
#include "dipper.h"

typedef uint64_t V __attribute__((vector_size(32)));
#define BS_LANES 4
#define BS_BLOCKS 256

static int BS_ROTSRC[64]; /* after rotation, slice i comes from slice BS_ROTSRC[i] */
static int BS_PDST[64];   /* permutation: slice i moves to BS_PDST[i] */

static void bs_init(void) {
    for (int w = 0; w < 4; w++) {
        int base = 48 - 16 * w, r = DIP_ROT[w];
        for (int k = 0; k < 16; k++) BS_ROTSRC[base + ((k + r) & 15)] = base + k;
    }
    for (int i = 0; i < 64; i++) BS_PDST[i] = dip_perm_bit(i);
}

static inline V vbcast(uint64_t x) { V v = {x, x, x, x}; return v; }

/* GIFT S-box on 4 slices, x0 = least significant bit */
static inline void bs_sbox(V *x0, V *x1, V *x2, V *x3) {
    V s0 = *x0, s1 = *x1, s2 = *x2, s3 = *x3;
    s1 ^= s0 & s2;
    s0 ^= s1 & s3;
    s2 ^= s0 | s1;
    s3 ^= s2;
    s1 ^= s3;
    s3 = ~s3;
    s2 ^= s0 & s1;
    *x0 = s3; *x1 = s1; *x2 = s2; *x3 = s0;
}

/* one round with round-key slices K (K == NULL means no key) */
static inline void bs_round(V *X, const V *K) {
    V T[64], U[64];
    if (K)
        for (int i = 0; i < 64; i++) X[i] ^= K[i];
    for (int j = 0; j < 16; j++) bs_sbox(&X[4 * j], &X[4 * j + 1], &X[4 * j + 2], &X[4 * j + 3]);
    for (int i = 0; i < 64; i++) T[i] = X[BS_ROTSRC[i]];
    /* A2 = A1 + B1 (slices 48.. and 32..), C2 = C1 + D1 (slices 16.. and 0..) */
    for (int h = 0; h < 2; h++) {
        int a = h ? 16 : 48, b = h ? 0 : 32;
        V c = vbcast(0);
        for (int k = 0; k < 16; k++) {
            V x = T[a + k], y = T[b + k], t = x ^ y;
            U[a + k] = t ^ c;
            c = (x & y) | (c & t);
        }
        for (int k = 0; k < 16; k++) U[b + k] = T[b + k];
    }
    for (int i = 0; i < 64; i++) X[BS_PDST[i]] = U[i];
}

/* broadcast a 64-bit round key into slices */
static inline void bs_key(uint64_t rk, V *K) {
    for (int i = 0; i < 64; i++) K[i] = vbcast((rk >> i) & 1 ? ~0ULL : 0);
}

/* ---- vector xoshiro256++ ---- */
typedef struct { V s[4]; } vrng_t;
static inline V vrotl(V x, int k) { return (x << k) | (x >> (64 - k)); }
static inline V vrng_next(vrng_t *r) {
    V *s = r->s;
    V res = vrotl(s[0] + s[3], 23) + s[0], t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = vrotl(s[3], 45);
    return res;
}
static void vrng_seed(vrng_t *r, uint64_t seed) {
    for (int i = 0; i < 4; i++)
        for (int l = 0; l < BS_LANES; l++) r->s[i][l] = splitmix(&seed);
}

static inline uint64_t vpopcnt(V x) {
    uint64_t s = 0;
    for (int l = 0; l < BS_LANES; l++) s += (uint64_t)__builtin_popcountll(x[l]);
    return s;
}

/* transpose helpers (tests only) */
static void bs_load(const uint64_t *blk, V *X) {
    for (int i = 0; i < 64; i++) {
        V v = vbcast(0);
        for (int b = 0; b < BS_BLOCKS; b++) v[b >> 6] |= ((blk[b] >> i) & 1ULL) << (b & 63);
        X[i] = v;
    }
}
static void bs_store(const V *X, uint64_t *blk) {
    for (int b = 0; b < BS_BLOCKS; b++) {
        uint64_t o = 0;
        for (int i = 0; i < 64; i++) o |= ((X[i][b >> 6] >> (b & 63)) & 1ULL) << i;
        blk[b] = o;
    }
}

#endif
