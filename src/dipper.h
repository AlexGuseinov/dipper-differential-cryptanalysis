/* Dipper v1.1: fast table-based implementation, a slow bit-level reference,
 * and the 96/128-bit key schedules. Header-only. */
#ifndef DIPPER_H
#define DIPPER_H
#pragma GCC diagnostic ignored "-Wunused-function"

#include <stdint.h>
#include <string.h>

static const uint8_t DIP_SBOX[16] = {0x1, 0xA, 0x4, 0xC, 0x6, 0xF, 0x3, 0x9,
                                     0x2, 0xD, 0xB, 0x7, 0x5, 0x0, 0x8, 0xE};
static const uint8_t DIP_SBOX_INV_INIT = 0; /* filled by dip_init */
static uint8_t DIP_SBOX_INV[16];

/* rotation constants of the round function; can be changed for step 10 */
static int DIP_ROT[4] = {1, 4, 7, 11};

static inline uint16_t rotl16(uint16_t x, int r) {
    r &= 15;
    return r ? (uint16_t)((x << r) | (x >> (16 - r))) : x;
}
static inline uint16_t rotr16(uint16_t x, int r) { return rotl16(x, 16 - (r & 15)); }

static inline int dip_perm_bit(int i) {
    return 4 * (i / 16) + 16 * ((3 * ((i % 16) / 4) + (i % 4)) % 4) + (i % 4);
}

/* tables: TW[w][x] = rotl16(sub16(x), DIP_ROT[w]); PT[k][b] = P applied to byte b at byte k */
static uint16_t DIP_TW[4][65536];
static uint64_t DIP_PT[8][256];
static uint64_t DIP_PINV_T[8][256];

static inline uint16_t dip_sub16(uint16_t x) {
    uint16_t o = 0;
    for (int j = 0; j < 4; j++) o |= (uint16_t)DIP_SBOX[(x >> (4 * j)) & 15] << (4 * j);
    return o;
}
static inline uint16_t dip_sub16_inv(uint16_t x) {
    uint16_t o = 0;
    for (int j = 0; j < 4; j++) o |= (uint16_t)DIP_SBOX_INV[(x >> (4 * j)) & 15] << (4 * j);
    return o;
}

static uint64_t dip_perm_slow(uint64_t s) {
    uint64_t o = 0;
    for (int i = 0; i < 64; i++) o |= ((s >> i) & 1ULL) << dip_perm_bit(i);
    return o;
}
static uint64_t dip_perm_inv_slow(uint64_t s) {
    uint64_t o = 0;
    for (int i = 0; i < 64; i++) o |= ((s >> dip_perm_bit(i)) & 1ULL) << i;
    return o;
}

static void dip_init_rot(int r0, int r1, int r2, int r3) {
    DIP_ROT[0] = r0; DIP_ROT[1] = r1; DIP_ROT[2] = r2; DIP_ROT[3] = r3;
    for (int x = 0; x < 16; x++) DIP_SBOX_INV[DIP_SBOX[x]] = (uint8_t)x;
    for (int w = 0; w < 4; w++)
        for (uint32_t x = 0; x < 65536; x++) DIP_TW[w][x] = rotl16(dip_sub16((uint16_t)x), DIP_ROT[w]);
    for (int k = 0; k < 8; k++)
        for (int b = 0; b < 256; b++) {
            DIP_PT[k][b] = dip_perm_slow((uint64_t)b << (8 * k));
            DIP_PINV_T[k][b] = dip_perm_inv_slow((uint64_t)b << (8 * k));
        }
    (void)DIP_SBOX_INV_INIT;
}
static __attribute__((unused)) void dip_init(void) { dip_init_rot(1, 4, 7, 11); }

static inline uint64_t dip_perm(uint64_t s) {
    return DIP_PT[0][s & 255] ^ DIP_PT[1][(s >> 8) & 255] ^ DIP_PT[2][(s >> 16) & 255] ^
           DIP_PT[3][(s >> 24) & 255] ^ DIP_PT[4][(s >> 32) & 255] ^ DIP_PT[5][(s >> 40) & 255] ^
           DIP_PT[6][(s >> 48) & 255] ^ DIP_PT[7][(s >> 56) & 255];
}
static inline uint64_t dip_perm_inv(uint64_t s) {
    return DIP_PINV_T[0][s & 255] ^ DIP_PINV_T[1][(s >> 8) & 255] ^ DIP_PINV_T[2][(s >> 16) & 255] ^
           DIP_PINV_T[3][(s >> 24) & 255] ^ DIP_PINV_T[4][(s >> 32) & 255] ^
           DIP_PINV_T[5][(s >> 40) & 255] ^ DIP_PINV_T[6][(s >> 48) & 255] ^
           DIP_PINV_T[7][(s >> 56) & 255];
}

/* keyless part of a round: SubCells, rotations, ARX, permutation */
static inline uint64_t dip_round_nokey(uint64_t s) {
    uint16_t A = DIP_TW[0][s >> 48], B = DIP_TW[1][(s >> 32) & 0xFFFF];
    uint16_t C = DIP_TW[2][(s >> 16) & 0xFFFF], D = DIP_TW[3][s & 0xFFFF];
    uint16_t A2 = (uint16_t)(A + B), C2 = (uint16_t)(C + D);
    return dip_perm(((uint64_t)A2 << 48) | ((uint64_t)B << 32) | ((uint64_t)C2 << 16) | D);
}
static inline uint64_t dip_round(uint64_t s, uint64_t rk) { return dip_round_nokey(s ^ rk); }

/* inverse of the keyless part (used to peel the last round) */
static inline uint64_t dip_round_nokey_inv(uint64_t s) {
    s = dip_perm_inv(s);
    uint16_t A2 = s >> 48, B = (s >> 32) & 0xFFFF, C2 = (s >> 16) & 0xFFFF, D = s & 0xFFFF;
    uint16_t A = (uint16_t)(A2 - B), C = (uint16_t)(C2 - D);
    A = dip_sub16_inv(rotr16(A, DIP_ROT[0]));
    B = dip_sub16_inv(rotr16(B, DIP_ROT[1]));
    C = dip_sub16_inv(rotr16(C, DIP_ROT[2]));
    D = dip_sub16_inv(rotr16(D, DIP_ROT[3]));
    return ((uint64_t)A << 48) | ((uint64_t)B << 32) | ((uint64_t)C << 16) | D;
}

/* bit-level reference round written directly from the specification text */
static uint64_t dip_round_ref(uint64_t s, uint64_t rk) {
    s ^= rk;
    uint64_t t = 0;
    for (int j = 0; j < 16; j++) t |= (uint64_t)DIP_SBOX[(s >> (4 * j)) & 15] << (4 * j);
    uint16_t A = rotl16(t >> 48, DIP_ROT[0]), B = rotl16((t >> 32) & 0xFFFF, DIP_ROT[1]);
    uint16_t C = rotl16((t >> 16) & 0xFFFF, DIP_ROT[2]), D = rotl16(t & 0xFFFF, DIP_ROT[3]);
    uint32_t A2 = ((uint32_t)A + B) % 65536u, C2 = ((uint32_t)C + D) % 65536u;
    uint64_t st = ((uint64_t)A2 << 48) | ((uint64_t)B << 32) | ((uint64_t)C2 << 16) | D;
    return dip_perm_slow(st);
}

/* ---------------- key schedules ---------------- */
/* key register as 16-bit words, k[0] = least significant */
typedef struct { uint16_t k[8]; int words; } dip_key;

static inline uint64_t dip_rk(const dip_key *K) {
    return ((uint64_t)K->k[3] << 48) | ((uint64_t)K->k[2] << 32) | ((uint64_t)K->k[1] << 16) | K->k[0];
}
static inline uint16_t sbox_top(uint16_t w) {
    return (uint16_t)((w & 0x0FFF) | ((uint16_t)DIP_SBOX[w >> 12] << 12));
}
static void dip_ks_update(dip_key *K, uint8_t rc) {
    uint16_t *k = K->k;
    if (K->words == 6) {
        uint16_t t1 = rotr16(k[1], 2), t0 = rotr16(k[0], 12);
        uint16_t n5 = k[3], n4 = k[2], n3 = t1, n2 = t0, n1 = k[5], n0 = k[4];
        k[5] = n5; k[4] = n4; k[3] = n3; k[2] = n2; k[1] = n1; k[0] = n0;
        /* S-boxes at K[95:92], K[71:68], K[47:44], K[23:20] */
        k[5] = sbox_top(k[5]);
        k[4] = (uint16_t)((k[4] & 0xFF0F) | ((uint16_t)DIP_SBOX[(k[4] >> 4) & 15] << 4)); /* bits 71:68 */
        k[2] = (uint16_t)((k[2] & 0x0FFF) | ((uint16_t)DIP_SBOX[(k[2] >> 12) & 15] << 12)); /* 47:44 */
        k[1] = (uint16_t)((k[1] & 0xFF0F) | ((uint16_t)DIP_SBOX[(k[1] >> 4) & 15] << 4)); /* 23:20 */
        k[3] ^= (uint16_t)(rc & 0x1F) << 0; /* K[52:48] -> bits 4:0 of k3 */
    } else {
        uint16_t t1 = rotr16(k[1], 2), t0 = rotr16(k[0], 12);
        uint16_t n7 = k[3], n6 = k[2], n5 = k[7], n4 = k[6], n3 = k[5], n2 = k[4], n1 = t1, n0 = t0;
        k[7] = n7; k[6] = n6; k[5] = n5; k[4] = n4; k[3] = n3; k[2] = n2; k[1] = n1; k[0] = n0;
        /* S-boxes at K[127:124], K[63:60], K[31:28] */
        k[7] = sbox_top(k[7]);
        k[3] = sbox_top(k[3]);
        k[1] = sbox_top(k[1]);
        k[4] ^= (uint16_t)(rc & 0x1F); /* K[68:64] */
    }
}
static inline uint8_t dip_lfsr_step(uint8_t l) {
    uint8_t b = ((l >> 4) ^ (l >> 1)) & 1;
    return (uint8_t)(((l << 1) | b) & 0x1F);
}
/* expand round keys RK_1..RK_R into rk[0..R-1] */
static void dip_expand(dip_key K, int R, uint64_t *rk) {
    uint8_t l = 1;
    for (int i = 0; i < R; i++) {
        rk[i] = dip_rk(&K);
        dip_ks_update(&K, l);
        l = dip_lfsr_step(l);
    }
}
static inline uint64_t dip_encrypt_rk(uint64_t p, const uint64_t *rk, int R) {
    for (int i = 0; i < R; i++) p = dip_round(p, rk[i]);
    return p;
}

/* ---------------- PRNG: xoshiro256** ---------------- */
typedef struct { uint64_t s[4]; } rng_t;
static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rng_next(rng_t *r) {
    uint64_t *s = r->s;
    uint64_t res = rotl64(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl64(s[3], 45);
    return res;
}
static inline uint64_t splitmix(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static void rng_seed(rng_t *r, uint64_t seed) {
    for (int i = 0; i < 4; i++) r->s[i] = splitmix(&seed);
}

#endif
