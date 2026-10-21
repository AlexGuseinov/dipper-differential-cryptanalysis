/* bsinv.h: bitsliced pieces for krec: one round with constant index maps (one state),
 * the inverse of the keyless round followed by P^-1, and the in-structure equality test.
 * Requires bitslice.h and bs_const.h (rotations 1, 4, 7, 11). */
#ifndef BSINV_H
#define BSINV_H

#include "bitslice.h"
#include "bs_const.h"

/* one keyed round, X -> X2 (as bs_round, with the index maps as constants) */
static inline __attribute__((always_inline)) void bsi_round(const V *restrict X, V *restrict X2, const V *restrict K) {
    V S[64];
    for (int j = 0; j < 16; j++) {
        V a0 = X[4 * j] ^ K[4 * j], a1 = X[4 * j + 1] ^ K[4 * j + 1], a2 = X[4 * j + 2] ^ K[4 * j + 2],
          a3 = X[4 * j + 3] ^ K[4 * j + 3];
        bs_sbox(&a0, &a1, &a2, &a3);
        S[4 * j] = a0; S[4 * j + 1] = a1; S[4 * j + 2] = a2; S[4 * j + 3] = a3;
    }
    for (int h = 0; h < 2; h++) {
        const int a = h ? 16 : 48, b = h ? 0 : 32;
        V c = vbcast(0);
        for (int k = 0; k < 16; k++) {
            V x = S[BSC_ROTSRC[a + k]], y = S[BSC_ROTSRC[b + k]], t = x ^ y;
            X2[BSC_PDST[a + k]] = t ^ c;
            X2[BSC_PDST[b + k]] = y;
            c = (x & y) | (c & t);
        }
    }
}

/* inverse GIFT S-box on 4 slices (undoes bs_sbox) */
static inline void bsi_sbox_inv(V *x0, V *x1, V *x2, V *x3) {
    V s3 = *x0, s1 = *x1, s2 = *x2, s0 = *x3;
    s2 ^= s0 & s1;
    s3 = ~s3;
    s1 ^= s3;
    s3 ^= s2;
    s2 ^= s0 | s1;
    s0 ^= s1 & s3;
    s1 ^= s0 & s2;
    *x0 = s0; *x1 = s1; *x2 = s2; *x3 = s3;
}

/* Z = P^-1(keyless round^-1 (C)) */
static inline void bsi_unpeel(const V *restrict C, V *restrict Z) {
    V W[64], T[64], Y[64];
    for (int i = 0; i < 64; i++) W[i] = C[BSC_PDST[i]];
    for (int h = 0; h < 2; h++) { /* A = A2 - B, C = C2 - D: x - y = x + ~y + 1 */
        const int a = h ? 16 : 48, b = h ? 0 : 32;
        V c = vbcast(~0ULL);
        for (int k = 0; k < 16; k++) {
            V x = W[a + k], y = ~W[b + k], t = x ^ y;
            T[a + k] = t ^ c;
            c = (x & y) | (c & t);
            T[b + k] = W[b + k];
        }
    }
    for (int i = 0; i < 64; i++) Y[BSC_ROTSRC[i]] = T[i];
    for (int j = 0; j < 16; j++) bsi_sbox_inv(&Y[4 * j], &Y[4 * j + 1], &Y[4 * j + 2], &Y[4 * j + 3]);
    for (int i = 0; i < 64; i++) Z[i] = Y[BSC_PDST[i]];
}

/* within 16-bit groups: x with bit q moved to bit q ^ d */
static inline V bsi_xswap(V x, int d) {
    if (d & 1) x = ((x & vbcast(0x5555555555555555ULL)) << 1) | ((x >> 1) & vbcast(0x5555555555555555ULL));
    if (d & 2) x = ((x & vbcast(0x3333333333333333ULL)) << 2) | ((x >> 2) & vbcast(0x3333333333333333ULL));
    if (d & 4) x = ((x & vbcast(0x0F0F0F0F0F0F0F0FULL)) << 4) | ((x >> 4) & vbcast(0x0F0F0F0F0F0F0F0FULL));
    if (d & 8) x = ((x & vbcast(0x00FF00FF00FF00FFULL)) << 8) | ((x >> 8) & vbcast(0x00FF00FF00FF00FFULL));
    return x;
}

/* M[d] bit q set iff blocks q and q ^ d (same 16-group) agree on every slice of Z except 27..31;
 * returns 1 if any M[d] is nonzero. Slices are tested in an order that kills random
 * coincidences quickly and stops early. */
static inline int bsi_match(const V *Z, V *M) {
    int any = 0;
    for (int d = 1; d < 16; d++) {
        V m = vbcast(~0ULL);
        for (int i = 0; i < 64; i++) {
            if (i >= 27 && i < 32) continue;
            m &= ~(Z[i] ^ bsi_xswap(Z[i], d));
            if ((i & 7) == 7) {
                V z = vbcast(0);
                if (!memcmp(&m, &z, sizeof m)) break;
            }
        }
        M[d] = m;
        for (int l = 0; l < BS_LANES; l++) any |= m[l] != 0;
    }
    return any;
}

#endif
