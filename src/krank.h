/* krank.h: ranking of the 32 last-round key bits guessed by the key-recovery attack (krec.c).
 *
 * A surviving pair gives z0, z1 = P^-1(y) with y = (keyless round R)^-1 (ciphertext),
 * so z = P^-1(x_{R-1}) ^ P^-1(rk_R). Survivors have equal z except in bits 27..31
 * (the top 5 bits of the C2 word). For a guess (kc, kd) of the C2 and D words of
 * P^-1(rk_R), round R-1 is inverted on the C half:
 *     T = (zc ^ kc) - (zd ^ kd),  C = S^-1(rotr(T, rC)),
 * and the pair counts if C0 ^ C1 = v << 4 for a target nibble value v (v = 8 is
 * Delta x_{R-2} = 0x800000). Both members share the low 11 bits of zc and all of zd, so
 * T0, T1 share their low 11 bits and the test depends only on h = kc >> 11 and
 * e = ((zd ^ kd) >> 11) + b (mod 32), b = [(zc ^ kc) & 0x7ff < (zd ^ kd) & 0x7ff].
 * Per pair a table V[e][h] is built once; the 2^32 keys are then scored as
 * (kd, kc & 0x7ff) x 32 values of h. kr_direct() is the plain formula (used for checks). */
#ifndef KRANK_H
#define KRANK_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "dipper.h"

typedef struct { uint64_t p0, p1, z0, z1, truth; } kr_pair; /* truth = Delta x_{R-2} (known-key check) */

static int kr_target(uint16_t dc, uint16_t tmask) {
    if (dc & ~0x00F0) return 0;
    int v = (dc >> 4) & 15;
    return v && ((tmask >> v) & 1);
}

static uint16_t kr_invC(uint16_t T) { return dip_sub16_inv(rotr16(T, DIP_ROT[2])); }

/* plain formula: does pair q count for key (kc, kd)? */
static int kr_direct(const kr_pair *q, uint16_t kc, uint16_t kd, uint16_t tmask) {
    uint16_t zd = (uint16_t)q->z0, zc0 = (uint16_t)(q->z0 >> 16), zc1 = (uint16_t)(q->z1 >> 16);
    uint16_t D = zd ^ kd;
    uint16_t C0 = kr_invC((uint16_t)((zc0 ^ kc) - D)), C1 = kr_invC((uint16_t)((zc1 ^ kc) - D));
    return kr_target(C0 ^ C1, tmask);
}

/* V[e*32 + h] for one pair (low part of T taken as 0: the test only involves the top 5 bits
 * of T as long as T0 and T1 agree on the low 11 bits, which survivors guarantee) */
static void kr_table(const kr_pair *q, uint16_t tmask, uint16_t *V) {
    uint16_t c0 = (uint16_t)(q->z0 >> 27) & 31, c1 = (uint16_t)(q->z1 >> 27) & 31;
    for (int e = 0; e < 32; e++)
        for (int h = 0; h < 32; h++) {
            uint16_t t0 = (uint16_t)((((c0 ^ h) - e) & 31) << 11), t1 = (uint16_t)((((c1 ^ h) - e) & 31) << 11);
            V[e * 32 + h] = (uint16_t)kr_target(kr_invC(t0) ^ kr_invC(t1), tmask);
        }
}

/* scores of the 32 keys (h << 11 | kclo, kd), h = 0..31, into s[32] */
static inline void kr_score32(const kr_pair *P, const uint16_t *V, size_t n, uint16_t kclo, uint16_t kd,
                              uint16_t *s) {
    uint16_t acc[32] = {0};
    for (size_t i = 0; i < n; i++) {
        uint16_t zd = (uint16_t)P[i].z0, zc = (uint16_t)(P[i].z0 >> 16);
        uint16_t D = zd ^ kd;
        int b = ((zc ^ kclo) & 0x7FF) < (D & 0x7FF);
        int e = ((D >> 11) + b) & 31;
        const uint16_t *v = V + (i * 1024) + e * 32;
        for (int h = 0; h < 32; h++) acc[h] += v[h];
    }
    memcpy(s, acc, sizeof acc);
}

#endif
