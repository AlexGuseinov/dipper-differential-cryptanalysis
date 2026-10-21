/* toysearch: exact single-nibble difference graph of Dipper-32 (toy32.h) for every
 * permutation variant (TPC in {1,3}, TPS in 0..3) and every rotation 4-tuple in [0,8)^4.
 *
 *   toysearch [TOP]
 *
 * For each configuration: the exact one-round transitions between the 120 differences with
 * one active nibble (exhaustive over each 16-bit half), the minimum mean weight cycle
 * (Karp), the best self-loop, and log2 of the spectral radius (cluster decay).
 * Prints the TOP (default 20) configurations with the slowest cluster decay and every
 * configuration whose 1-nibble graph matches Dipper's: self-loop 2^-3 and decay
 * log2(32/7) = 2.1926 bits per round. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "toy32.h"

/* HD[r1][r2][j][o]: count over 2^16 half values of output difference o (16 bits) for the
 * half input difference with one active nibble j = 4*pos + (v-1)... stored sparsely */
typedef struct { uint16_t o; uint32_t c; } hd_ent;
static hd_ent *HD[8][8][64]; static int HN[8][8][64];

static uint16_t hin(int j) { int pos = j / 15, v = j % 15 + 1; return (uint16_t)(v << (4 * pos)); }

static void half_tables(void) {
    static uint32_t cnt[65536];
    for (int r1 = 0; r1 < 8; r1++)
        for (int r2 = 0; r2 < 8; r2++) {
            uint16_t T[65536];
            for (uint32_t x = 0; x < 65536; x++) {
                uint8_t a = rotl8(tsub8((uint8_t)(x >> 8)), r1), b = rotl8(tsub8((uint8_t)x), r2);
                T[x] = (uint16_t)(((uint8_t)(a + b) << 8) | b);
            }
            for (int j = 0; j < 60; j++) {
                uint16_t d = hin(j);
                memset(cnt, 0, sizeof cnt);
                for (uint32_t x = 0; x < 65536; x++) cnt[(uint16_t)(T[x] ^ T[x ^ d])]++;
                int n = 0;
                for (int o = 0; o < 65536; o++) n += cnt[o] != 0;
                HD[r1][r2][j] = malloc(n * sizeof(hd_ent)); HN[r1][r2][j] = n; n = 0;
                for (int o = 0; o < 65536; o++) if (cnt[o]) { HD[r1][r2][j][n].o = (uint16_t)o; HD[r1][r2][j][n].c = cnt[o]; n++; }
            }
        }
}

static int nib1(uint32_t d) { /* index of the single active nibble difference, or -1 */
    int c = 0, idx = -1;
    for (int k = 0; k < 8; k++) { int v = (d >> (4 * k)) & 15; if (v) { c++; idx = 15 * k + v - 1; } }
    return c == 1 ? idx : -1;
}

typedef struct { int pc, ps, r[4]; double rho, cyc, loop; uint32_t loopd; } cfg_res;

int main(int argc, char **argv) {
    int TOP = argc > 1 ? atoi(argv[1]) : 20;
    dip_init();
    half_tables();
    static double M[120][120];
    cfg_res *res = malloc(sizeof(cfg_res) * 8 * 4096);
    int nr = 0;
    for (int pc = 1; pc <= 3; pc += 2)
        for (int ps = 0; ps < 4; ps++)
            for (int rr = 0; rr < 4096; rr++) {
                int r[4] = {rr >> 9, (rr >> 6) & 7, (rr >> 3) & 7, rr & 7};
                TPC = pc; TPS = ps;
                memset(M, 0, sizeof M);
                for (int u = 0; u < 120; u++) {
                    int k = u / 15, v = u % 15 + 1, h = k >= 4; /* nibbles 4..7 = A,B half */
                    int j = 15 * (k & 3) + v - 1;
                    int r1 = h ? r[0] : r[2], r2 = h ? r[1] : r[3];
                    for (int e = 0; e < HN[r1][r2][j]; e++) {
                        uint32_t o = (uint32_t)HD[r1][r2][j][e].o << (h ? 16 : 0);
                        uint32_t po = 0;
                        for (int i = 0; i < 32; i++) po |= ((o >> i) & 1u) << toy_perm_bit(i);
                        int t = nib1(po);
                        if (t >= 0) M[u][t] += HD[r1][r2][j][e].c / 65536.0;
                    }
                }
                /* best self-loop */
                double loop = 0; uint32_t ld = 0;
                for (int u = 0; u < 120; u++) if (M[u][u] > loop) { loop = M[u][u]; ld = (uint32_t)(u % 15 + 1) << (4 * (u / 15)); }
                /* edges */
                static int EU[14400], EV[14400]; static double EW[14400], EP[14400];
                int ne = 0;
                for (int u = 0; u < 120; u++) for (int t = 0; t < 120; t++)
                    if (M[u][t] > 0) { EU[ne] = u; EV[ne] = t; EP[ne] = M[u][t]; EW[ne] = -log2(M[u][t]); ne++; }
                /* Karp minimum mean cycle on weights -log2 p */
                static double Dk[121][120];
                double best = INFINITY;
                for (int s0 = 0; s0 < 120; s0++) Dk[0][s0] = 0;
                for (int k = 1; k <= 120; k++) {
                    for (int v = 0; v < 120; v++) Dk[k][v] = INFINITY;
                    for (int e = 0; e < ne; e++) {
                        double w = Dk[k - 1][EU[e]] + EW[e];
                        if (w < Dk[k][EV[e]]) Dk[k][EV[e]] = w;
                    }
                }
                for (int v = 0; v < 120; v++) {
                    if (Dk[120][v] == INFINITY) continue;
                    double mx = -INFINITY;
                    for (int k = 0; k < 120; k++)
                        if (Dk[k][v] < INFINITY) { double q = (Dk[120][v] - Dk[k][v]) / (120 - k); if (q > mx) mx = q; }
                    if (mx < best) best = mx;
                }
                /* spectral radius */
                double vv[120], ww[120], lam = 0;
                for (int i = 0; i < 120; i++) vv[i] = 1;
                for (int it = 0; it < 1000; it++) {
                    for (int i = 0; i < 120; i++) ww[i] = 0;
                    for (int e = 0; e < ne; e++) ww[EV[e]] += vv[EU[e]] * EP[e];
                    double sm = 0;
                    for (int i = 0; i < 120; i++) sm = ww[i] > sm ? ww[i] : sm;
                    if (sm == 0) { lam = 0; break; }
                    for (int i = 0; i < 120; i++) vv[i] = ww[i] / sm;
                    lam = sm;
                }
                cfg_res *c = &res[nr++];
                c->pc = pc; c->ps = ps; memcpy(c->r, r, sizeof r);
                c->rho = lam > 0 ? log2(lam) : -INFINITY; c->cyc = best; c->loop = loop > 0 ? log2(loop) : -INFINITY; c->loopd = ld;
            }
    /* sort by rho descending */
    for (int i = 1; i < nr; i++) { cfg_res x = res[i]; int j = i; while (j > 0 && res[j - 1].rho < x.rho) { res[j] = res[j - 1]; j--; } res[j] = x; }
    printf("configurations: %d (TPC, TPS, rotations A B C D: log2 rho, min mean cycle, best self-loop)\n", nr);
    printf("top %d by cluster decay:\n", TOP);
    for (int i = 0; i < TOP && i < nr; i++)
        printf("  pc=%d ps=%d rot=%d,%d,%d,%d rho=%.4f cycle=%.4f loop=%.4f at %08x\n", res[i].pc, res[i].ps, res[i].r[0],
               res[i].r[1], res[i].r[2], res[i].r[3], res[i].rho, res[i].cyc, res[i].loop, res[i].loopd);
    printf("Dipper's rotations mod 8 (1,4,7,3), all permutation variants:\n");
    for (int i = 0; i < nr; i++)
        if (res[i].r[0] == 1 && res[i].r[1] == 4 && res[i].r[2] == 7 && res[i].r[3] == 3)
            printf("  pc=%d ps=%d rho=%.4f cycle=%.4f loop=%.4f at %08x\n", res[i].pc, res[i].ps, res[i].rho, res[i].cyc,
                   res[i].loop, res[i].loopd);
    printf("configurations matching Dipper's 1-nibble graph (loop -3, rho -2.1926):\n");
    int m = 0;
    for (int i = 0; i < nr; i++)
        if (fabs(res[i].loop + 3) < 1e-9 && fabs(res[i].rho + 2.1926) < 5e-5) {
            if (m < 200)
                printf("  pc=%d ps=%d rot=%d,%d,%d,%d rho=%.4f cycle=%.4f loop at %08x\n", res[i].pc, res[i].ps, res[i].r[0],
                       res[i].r[1], res[i].r[2], res[i].r[3], res[i].rho, res[i].cyc, res[i].loopd);
            m++;
        }
    printf("matching: %d\n", m);
    return 0;
}
