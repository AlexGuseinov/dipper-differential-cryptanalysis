/* carrycheck: exhaustive check of the paper's Proposition on 0x800000 -> 0x800000.
 *
 * Input difference 0x80 on word C (bit 7), so the half is
 *   u = S(a) <<< 7, u' = S(a ^ 0x80) <<< 7, v uniform (the rotated S-box output of word D).
 * Over all 2^32 values of (a, v): the output difference of the addition is 1 << 11 exactly when
 * (i) the S-box output difference of nibble 1 of a is 0x3 and (ii) v_11 != c_11, the carry
 * into bit 11 of u + v. Also prints the S-box outputs per output difference. */
#include <inttypes.h>
#include <stdio.h>
#include <stdint.h>

static const int SB[16] = {1, 10, 4, 12, 6, 15, 3, 9, 2, 13, 11, 7, 5, 0, 8, 14};

static uint16_t sw(uint16_t x) {
    uint16_t r = 0;
    for (int i = 0; i < 4; i++) r |= (uint16_t)(SB[(x >> (4 * i)) & 15] << (4 * i));
    return r;
}

static uint16_t rotl(uint16_t x, int r) { return (uint16_t)((x << r) | (x >> (16 - r))); }

int main(void) {
    uint64_t hits = 0, mismatches = 0, by[16] = {0};
    for (uint32_t a = 0; a < 65536; a++) {
        uint16_t u = rotl(sw((uint16_t)a), 7), w = rotl(sw((uint16_t)(a ^ 0x80)), 7);
        int n = (a >> 4) & 15, sg = SB[n] ^ SB[n ^ 8];
        for (uint32_t v = 0; v < 65536; v++) {
            int hit = (((uint16_t)(u + v)) ^ ((uint16_t)(w + v))) == (1u << 11);
            uint32_t c11 = ((u & 0x7FFu) + (v & 0x7FFu)) >> 11;
            int pred = sg == 3 && (((v >> 11) ^ c11) & 1);
            hits += (uint64_t)hit;
            by[sg] += (uint64_t)hit;
            mismatches += (uint64_t)(hit != pred);
        }
    }
    printf("inputs 2^32, output difference 1<<11: %" PRIu64 " (2^32/%.0f)\n", hits, 4294967296.0 / (double)hits);
    printf("mismatches with (i) and (ii): %" PRIu64 "\n", mismatches);
    for (int s = 0; s < 16; s++) {
        int any = 0;
        for (int n = 0; n < 16; n++) any |= (SB[n] ^ SB[n ^ 8]) == s;
        if (!any) continue;
        printf("sigma=%x hits=%" PRIu64 " S(n):", s, by[s]);
        for (int n = 0; n < 16; n++) if ((SB[n] ^ SB[n ^ 8]) == s) printf(" %x", SB[n]);
        printf("\n");
    }
    return mismatches != 0;
}
