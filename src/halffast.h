/* halffast.h: exact output distribution of one 32-bit half of the Dipper round
 * (added word a, passed word b) for an input difference (da, db), without
 * enumerating all 2^32 values.
 *
 * Half: u = TA[A] = rotl(S(A), rA), v = TB[B]; outputs (u + v, v).
 * For a pair (A, A^da), (B, B^db) the output difference is
 *     ((u + v) ^ (u2 + v2), v ^ v2).
 * Condition on the values of the active nibbles of A and B (16^(kA+kB) cases).
 * The inactive nibbles contribute free bits: u = uc + uf, v = vc + vf (disjoint bit
 * sets, so + is exact), u2 = u2c + uf, v2 = v2c + vf. With s = u + v:
 *     u2 + v2 = s + D,  D = (u2c - uc) + (v2c - vc)   (mod 2^16)
 * so the sum difference is s ^ (s + D), and s = c + t with c = uc + vc and
 * t = uf + vf distributed as G(t) = #{(uf, vf) : uf + vf = t}, which depends only on
 * the free-bit masks. Counts are out of 2^32, as in lowgraph's half_dist. */
#ifndef HALFFAST_H
#define HALFFAST_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "dipper.h"

typedef struct { uint32_t out; uint64_t num; } hf_ent;

/* fixed parts uc (value) and u2c (value ^ diff) for active-nibble values packed in x */
static void hf_fixed(int w, uint16_t d, uint32_t x, uint16_t *uc, uint16_t *u2c) {
    uint16_t a = 0, b = 0;
    int q = 0;
    for (int j = 0; j < 4; j++) {
        int dj = (d >> (4 * j)) & 15;
        if (!dj) continue;
        int xj = (x >> (4 * q++)) & 15;
        a |= (uint16_t)(DIP_SBOX[xj] << (4 * j));
        b |= (uint16_t)(DIP_SBOX[xj ^ dj] << (4 * j));
    }
    *uc = rotl16(a, DIP_ROT[w]);
    *u2c = rotl16(b, DIP_ROT[w]);
}

static int hf_nact(uint16_t d) {
    int c = 0;
    for (int j = 0; j < 4; j++) c += ((d >> (4 * j)) & 15) != 0;
    return c;
}

/* free-bit mask of word w for input difference d (bits of the inactive nibbles, rotated) */
static uint16_t hf_free(int w, uint16_t d) {
    uint16_t m = 0;
    for (int j = 0; j < 4; j++)
        if (!((d >> (4 * j)) & 15)) m |= (uint16_t)(0xF << (4 * j));
    return rotl16(m, DIP_ROT[w]);
}

/* G[t] = #{uf subset-of FU, vf subset-of FV : uf + vf = t mod 2^16}; returns 1 if uniform */
static int hf_G(uint16_t FU, uint16_t FV, uint64_t *G) {
    int nu = __builtin_popcount(FU), nv = __builtin_popcount(FV);
    if (nu == 16 || nv == 16) {
        uint64_t c = 1ULL << (nu + nv - 16);
        for (int t = 0; t < 65536; t++) G[t] = c;
        return 1;
    }
    memset(G, 0, 65536 * sizeof(uint64_t));
    /* enumerate submasks of FV once, then of FU */
    uint16_t *vs = malloc(sizeof(uint16_t) << nv);
    size_t n = 0;
    uint16_t s = 0;
    do { vs[n++] = s; s = (uint16_t)((s - FV) & FV); } while (s);
    uint16_t a = 0;
    do {
        for (size_t i = 0; i < n; i++) G[(uint16_t)(a + vs[i])]++;
        a = (uint16_t)((a - FU) & FU);
    } while (a);
    free(vs);
    return 0;
}

static int hf_cmp(const void *x, const void *y) {
    const hf_ent *a = x, *b = y;
    return a->num < b->num ? 1 : a->num > b->num ? -1 : (a->out > b->out) - (a->out < b->out);
}

typedef struct { uint16_t c, D; uint32_t w; } hf_cond;
static int hf_cmp_cond(const void *x, const void *y) {
    const hf_cond *a = x, *b = y;
    if (a->D != b->D) return a->D < b->D ? -1 : 1;
    return (a->c > b->c) - (a->c < b->c);
}

/* distribution for words (wa added, wb passed), input (da, db); keeps counts >= minnum.
 * Output entries: out = (sum difference << 16) | (passed-word difference), sorted by count. */
static hf_ent *half_fast(int wa, int wb, uint16_t da, uint16_t db, uint64_t minnum, size_t *nout) {
    if (!da && !db) {
        hf_ent *r = malloc(sizeof *r);
        r->out = 0; r->num = 1ULL << 32; *nout = 1;
        return r;
    }
    int kA = hf_nact(da), kB = hf_nact(db);
    uint64_t *G = malloc(65536 * sizeof(uint64_t));
    int uni = hf_G(hf_free(wa, da), hf_free(wb, db), G);
    uint32_t nA = 1u << (4 * kA), nB = 1u << (4 * kB);
    /* group conditions by the passed-word difference ob */
    uint16_t *obv = malloc(nB * sizeof(uint16_t)), *vcv = malloc(nB * sizeof(uint16_t)),
             *dvv = malloc(nB * sizeof(uint16_t));
    for (uint32_t xb = 0; xb < nB; xb++) {
        uint16_t vc, v2c;
        hf_fixed(wb, db, xb, &vc, &v2c);
        obv[xb] = vc ^ v2c; vcv[xb] = vc; dvv[xb] = (uint16_t)(v2c - vc);
    }
    uint16_t *ucv = malloc(nA * sizeof(uint16_t)), *duv = malloc(nA * sizeof(uint16_t));
    for (uint32_t xa = 0; xa < nA; xa++) {
        uint16_t uc, u2c;
        hf_fixed(wa, da, xa, &uc, &u2c);
        ucv[xa] = uc; duv[xa] = (uint16_t)(u2c - uc);
    }
    uint8_t *seen = calloc(65536, 1);
    uint64_t *cnt = malloc(65536 * sizeof(uint64_t));
    hf_cond *cd = malloc((size_t)nA * nB * sizeof(hf_cond));
    size_t n = 0, cap = 256;
    hf_ent *res = malloc(cap * sizeof *res);
    for (uint32_t xb0 = 0; xb0 < nB; xb0++) {
        uint16_t ob = obv[xb0];
        if (seen[ob]) continue;
        seen[ob] = 1;
        size_t nc = 0;
        for (uint32_t xb = xb0; xb < nB; xb++) {
            if (obv[xb] != ob) continue;
            for (uint32_t xa = 0; xa < nA; xa++) {
                cd[nc].c = uni ? 0 : (uint16_t)(ucv[xa] + vcv[xb]);
                cd[nc].D = (uint16_t)(duv[xa] + dvv[xb]);
                cd[nc].w = 1; nc++;
            }
        }
        /* merge equal (D, c) */
        qsort(cd, nc, sizeof *cd, hf_cmp_cond);
        size_t m = 0;
        for (size_t i = 0; i < nc; i++) {
            if (m && cd[m - 1].D == cd[i].D && cd[m - 1].c == cd[i].c) cd[m - 1].w += cd[i].w;
            else cd[m++] = cd[i];
        }
        memset(cnt, 0, 65536 * sizeof(uint64_t));
        for (size_t i = 0; i < m; i++) {
            uint16_t c = cd[i].c, D = cd[i].D;
            uint64_t w = cd[i].w;
            for (uint32_t t = 0; t < 65536; t++) {
                uint16_t s = (uint16_t)(c + t);
                cnt[(uint16_t)(s ^ (uint16_t)(s + D))] += w * G[t];
            }
        }
        for (uint32_t o = 0; o < 65536; o++)
            if (cnt[o] && cnt[o] >= minnum) {
                if (n == cap) { cap *= 2; res = realloc(res, cap * sizeof *res); }
                res[n].out = (o << 16) | ob; res[n].num = cnt[o]; n++;
            }
    }
    free(G); free(obv); free(vcv); free(dvv); free(ucv); free(duv); free(seen); free(cnt); free(cd);
    qsort(res, n, sizeof *res, hf_cmp);
    *nout = n;
    return res;
}

#endif
