/* hfcheck: half_fast (halffast.h) against the brute-force 2^32 enumeration of lowgraph.
 *   hfcheck N SEED [K] [fast]   N random (da, db) with at most K (default 3) active nibbles, both halves;
 *   with "fast" only half_fast is run (timing).
 * Prints one line per case and "mismatches=M". */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "halffast.h"
typedef hf_ent hent_t;
static uint64_t bminnum = 1;
static int cmp_hent(const void *a, const void *b) { return hf_cmp(a, b); }
static hent_t *brute_dist(int wa, int wb, uint16_t da, uint16_t db, size_t *nout) {
    if (!da && !db) {
        hent_t *r = malloc(sizeof *r);
        r->out = 0; r->num = 1ULL << 32; *nout = 1;
        return r;
    }
    const uint16_t *TA = DIP_TW[wa], *TB = DIP_TW[wb];
    uint64_t *cnt = malloc(65536 * sizeof(uint64_t));
    /* group B by its output difference */
    uint32_t *ord = malloc(65536 * sizeof(uint32_t));
    uint32_t start[65537] = {0};
    static __thread uint32_t key[65536];
    for (uint32_t B = 0; B < 65536; B++) { key[B] = (uint16_t)(TB[B] ^ TB[B ^ db]); start[key[B] + 1]++; }
    for (int i = 0; i < 65536; i++) start[i + 1] += start[i];
    uint32_t pos[65536];
    memcpy(pos, start, sizeof pos);
    for (uint32_t B = 0; B < 65536; B++) ord[pos[key[B]]++] = B;
    int hbB = db ? 31 - __builtin_clz(db) : -1, hbA = da ? 31 - __builtin_clz(da) : -1;
    size_t n = 0, cap = 256;
    hent_t *res = malloc(cap * sizeof *res);
    for (uint32_t ob = 0; ob < 65536; ob++) {
        if (start[ob] == start[ob + 1]) continue;
        memset(cnt, 0, 65536 * sizeof(uint64_t));
        for (uint32_t g = start[ob]; g < start[ob + 1]; g++) {
            uint32_t B = ord[g];
            if (hbB >= 0 && ((B >> hbB) & 1)) continue; /* visit B < B^db only */
            uint16_t v = TB[B], v2 = TB[B ^ db];
            if (hbB >= 0) {
                for (uint32_t A = 0; A < 65536; A++) {
                    uint16_t u = TA[A], u2 = TA[A ^ da];
                    cnt[(uint16_t)((uint16_t)(u + v) ^ (uint16_t)(u2 + v2))]++;
                }
            } else { /* db == 0: halve over A */
                uint32_t lo = (1u << hbA) - 1;
                for (uint32_t a = 0; a < 32768; a++) {
                    uint32_t A = ((a & ~lo) << 1) | (a & lo);
                    uint16_t u = TA[A], u2 = TA[A ^ da];
                    cnt[(uint16_t)((uint16_t)(u + v) ^ (uint16_t)(u2 + v2))]++;
                }
            }
        }
        for (uint32_t o = 0; o < 65536; o++) {
            uint64_t c = 2 * cnt[o];
            if (c && c >= bminnum) {
                if (n == cap) { cap *= 2; res = realloc(res, cap * sizeof *res); }
                res[n].out = (o << 16) | ob; res[n].num = c; n++;
            }
        }
    }
    free(cnt); free(ord);
    qsort(res, n, sizeof *res, cmp_hent);
    *nout = n;
    return res;
}


static int nibs(uint32_t v) { int c = 0; for (int j = 0; j < 8; j++) c += ((v >> (4 * j)) & 15) != 0; return c; }
int main(int argc, char **argv) {
    int N = argc > 1 ? atoi(argv[1]) : 4, K = argc > 3 ? atoi(argv[3]) : 3;
    uint64_t seed = argc > 2 ? strtoull(argv[2], 0, 10) : 1;
    dip_init();
    rng_t r; rng_seed(&r, seed);
    int bad = 0, fastonly = argc > 4 && !strcmp(argv[4], "fast");
    for (int i = 0; i < N; i++) {
        uint32_t v;
        do { v = 0; int k = 1 + (int)(rng_next(&r) % (uint64_t)K);
             for (int q = 0; q < k; q++) { int j = (int)(rng_next(&r) % 8); v &= ~(0xFu << (4 * j)); v |= (uint32_t)(1 + rng_next(&r) % 15) << (4 * j); }
        } while (!v || nibs(v) > K);
        int h = (int)(rng_next(&r) & 1), wa = h ? 2 : 0, wb = h ? 3 : 1;
        uint16_t da = v >> 16, db = v & 0xFFFF;
        size_t n1, n2;
        hf_ent *a = half_fast(wa, wb, da, db, 1, &n1);
        if (fastonly) { printf("case %d da=%04x db=%04x entries=%zu\n", i, da, db, n1); free(a); continue; }
        hf_ent *b = brute_dist(wa, wb, da, db, &n2);
        int ok = n1 == n2;
        for (size_t q = 0; ok && q < n1; q++) ok = a[q].out == b[q].out && a[q].num == b[q].num;
        printf("case %d half=%c da=%04x db=%04x entries=%zu/%zu %s\n", i, h ? 'R' : 'L', da, db, n1, n2, ok ? "equal" : "DIFFERENT");
        fflush(stdout);
        bad += !ok;
        free(a); free(b);
    }
    printf("mismatches=%d\n", bad);
    return bad != 0;
}
