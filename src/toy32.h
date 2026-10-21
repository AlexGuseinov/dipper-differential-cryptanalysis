/* toy32.h: Dipper-32, a 32-bit small-scale variant of Dipper for direct measurements.
 *
 * State: four 8-bit words A | B | C | D (A = bits 31..24, D = bits 7..0), 8 nibbles.
 * Round (as Dipper, words of 8 bits instead of 16):
 *   x ^= rk; S-box (Dipper's 4-bit S-box) on every nibble; rotate word w left by TROT[w]
 *   (w = 0..3 for A, B, C, D); A = A + B, C = C + D (mod 2^8); bit permutation TP.
 * Bit permutation (GIFT-like, as Dipper's: the 4 bits of a nibble go to 4 different words):
 *   bit i = 8q + 4a + b (q = word index from D = 0 to A = 3, a = nibble in word, b = bit)
 *   goes to word ((TPC * b + 2a + (q >> 1) + TPS) mod 4), nibble (q & 1), bit b.
 * Key schedule Dipper-32/48 (scaled from Dipper-64/96): 48-bit register of six 8-bit words
 *   k5..k0, round key = k3 k2 k1 k0; update: (k5..k0) <- (k3, k2, rotr(k1, 1), rotr(k0, 6), k5, k4),
 *   then S-box on the high nibble of k5, the low nibble of k4, the high nibble of k2, the low
 *   nibble of k1, and the 5-bit LFSR round constant (Dipper's) XORed into the low bits of k3. */
#ifndef TOY32_H
#define TOY32_H

#include <stdint.h>
#include <string.h>
#include "dipper.h"

static int TROT[4] = {1, 4, 7, 3};
static int TPC = 1, TPS = 0;
static uint16_t THT[2][65536]; /* half tables: (A,B) -> (A2,B) and (C,D) -> (C2,D) after S-box, rotations, add */
static uint32_t TPT[4][256];    /* permutation by bytes */
static uint32_t TPINV_T[4][256];

static inline uint8_t rotl8(uint8_t x, int r) { r &= 7; return r ? (uint8_t)((x << r) | (x >> (8 - r))) : x; }
static inline uint8_t rotr8(uint8_t x, int r) { return rotl8(x, 8 - (r & 7)); }
static inline uint8_t tsub8(uint8_t x) { return (uint8_t)(DIP_SBOX[x & 15] | (DIP_SBOX[x >> 4] << 4)); }

static inline int toy_perm_bit(int i) {
    int q = i >> 3, a = (i >> 2) & 1, b = i & 3;
    int w = (TPC * b + 2 * a + (q >> 1) + TPS) & 3;
    return 8 * w + 4 * (q & 1) + b;
}
static uint32_t toy_perm_slow(uint32_t s) {
    uint32_t o = 0;
    for (int i = 0; i < 32; i++) o |= ((s >> i) & 1u) << toy_perm_bit(i);
    return o;
}

static void toy_init(int ra, int rb, int rc, int rd, int pc, int ps) {
    TROT[0] = ra; TROT[1] = rb; TROT[2] = rc; TROT[3] = rd; TPC = pc; TPS = ps;
    dip_init();
    for (int h = 0; h < 2; h++)
        for (uint32_t x = 0; x < 65536; x++) {
            uint8_t a = rotl8(tsub8((uint8_t)(x >> 8)), TROT[2 * h]), b = rotl8(tsub8((uint8_t)x), TROT[2 * h + 1]);
            THT[h][x] = (uint16_t)(((uint8_t)(a + b) << 8) | b);
        }
    for (int k = 0; k < 4; k++)
        for (int v = 0; v < 256; v++) {
            uint32_t s = (uint32_t)v << (8 * k), o = toy_perm_slow(s), inv = 0;
            TPT[k][v] = o;
            for (int i = 0; i < 32; i++) if ((s >> toy_perm_bit(i)) & 1) inv |= 1u << i;
            TPINV_T[k][v] = inv;
        }
}

static inline uint32_t toy_perm(uint32_t s) {
    return TPT[0][s & 255] ^ TPT[1][(s >> 8) & 255] ^ TPT[2][(s >> 16) & 255] ^ TPT[3][s >> 24];
}
static inline uint32_t toy_perm_inv(uint32_t s) {
    return TPINV_T[0][s & 255] ^ TPINV_T[1][(s >> 8) & 255] ^ TPINV_T[2][(s >> 16) & 255] ^ TPINV_T[3][s >> 24];
}
static inline uint32_t toy_round_nokey(uint32_t s) {
    return toy_perm(((uint32_t)THT[0][s >> 16] << 16) | THT[1][s & 0xFFFF]);
}
static inline uint32_t toy_round(uint32_t s, uint32_t rk) { return toy_round_nokey(s ^ rk); }

/* bit-level reference round, written from the description above */
static uint32_t toy_round_ref(uint32_t s, uint32_t rk) {
    s ^= rk;
    uint8_t w[4];
    for (int q = 0; q < 4; q++) w[q] = (uint8_t)(s >> (24 - 8 * q)); /* w[0] = A .. w[3] = D */
    for (int q = 0; q < 4; q++) w[q] = rotl8(tsub8(w[q]), TROT[q]);
    w[0] = (uint8_t)(w[0] + w[1]);
    w[2] = (uint8_t)(w[2] + w[3]);
    uint32_t t = ((uint32_t)w[0] << 24) | ((uint32_t)w[1] << 16) | ((uint32_t)w[2] << 8) | w[3];
    return toy_perm_slow(t);
}

/* keyless round inverse (used to peel the last round) */
static inline uint32_t toy_round_nokey_inv(uint32_t s) {
    s = toy_perm_inv(s);
    uint8_t A2 = (uint8_t)(s >> 24), B = (uint8_t)(s >> 16), C2 = (uint8_t)(s >> 8), D = (uint8_t)s;
    uint8_t w[4] = {(uint8_t)(A2 - B), B, (uint8_t)(C2 - D), D};
    uint32_t o = 0;
    for (int q = 0; q < 4; q++) {
        uint8_t v = rotr8(w[q], TROT[q]);
        v = (uint8_t)(DIP_SBOX_INV[v & 15] | (DIP_SBOX_INV[v >> 4] << 4));
        o |= (uint32_t)v << (24 - 8 * q);
    }
    return o;
}

/* key schedule Dipper-32/48; k[0] = least significant word */
typedef struct { uint8_t k[6]; } toy_key;
static void toy_expand(toy_key K, int R, uint32_t *rk) {
    uint8_t l = 1;
    for (int i = 0; i < R; i++) {
        uint8_t *k = K.k;
        rk[i] = ((uint32_t)k[3] << 24) | ((uint32_t)k[2] << 16) | ((uint32_t)k[1] << 8) | k[0];
        uint8_t n5 = k[3], n4 = k[2], n3 = rotr8(k[1], 1), n2 = rotr8(k[0], 6), n1 = k[5], n0 = k[4];
        k[5] = n5; k[4] = n4; k[3] = n3; k[2] = n2; k[1] = n1; k[0] = n0;
        k[5] = (uint8_t)((k[5] & 0x0F) | (DIP_SBOX[k[5] >> 4] << 4));
        k[4] = (uint8_t)((k[4] & 0xF0) | DIP_SBOX[k[4] & 15]);
        k[2] = (uint8_t)((k[2] & 0x0F) | (DIP_SBOX[k[2] >> 4] << 4));
        k[1] = (uint8_t)((k[1] & 0xF0) | DIP_SBOX[k[1] & 15]);
        k[3] ^= (uint8_t)(l & 0x1F);
        l = dip_lfsr_step(l);
    }
}
static inline uint32_t toy_encrypt(uint32_t p, const uint32_t *rk, int R) {
    for (int i = 0; i < R; i++) p = toy_round(p, rk[i]);
    return p;
}

#endif
