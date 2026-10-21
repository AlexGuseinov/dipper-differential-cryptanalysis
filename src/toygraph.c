/* toygraph: Markov model of Dipper-32 (toy32.h) for the full-codebook validation (step 15).
 *
 *   toygraph -k K -m MIN_LOG2P -d DMAX -o PREFIX [-t TOP]
 *
 * Exact one-round transitions (each 16-bit half enumerated over its 2^16 values) between all
 * differences with 1..K active nibbles, edges with probability >= 2^MIN_LOG2P.
 * Then, under independent round keys, lower bounds on Pr[din -> dout, d rounds] by
 * propagation in this graph, for din = 0x8000 and for the 15 nibble-3 differences v << 12.
 * Writes PREFIX_d<d>.txt for d = 1..DMAX: the TOP outputs ranked by the structure signal
 *   sig(dout) = sum_v Pr[v<<12 -> dout],
 * one line each: "dout sig p8000" with p8000 = Pr[0x8000 -> dout].
 * Prints min mean cycle and spectral radius (graphalg.h). */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "toy32.h"
#include "graphalg.h"

typedef struct { uint32_t img; uint32_t num; uint8_t mask; } tent; /* permuted half output */
static int K = 3;
static double MINLOG = -20;

static int nact16(uint16_t v) { int c = 0; for (int j = 0; j < 4; j++) c += ((v >> (4 * j)) & 15) != 0; return c; }
static uint8_t nmask(uint32_t d) { uint8_t m = 0; for (int j = 0; j < 8; j++) if ((d >> (4 * j)) & 15) m |= (uint8_t)(1u << j); return m; }

static uint16_t *HIN; static int NHIN;
static tent *HE[2][65536]; static int HNE[2][65536];

static int cmp_tent(const void *a, const void *b) {
    const tent *x = a, *y = b;
    return x->num < y->num ? 1 : x->num > y->num ? -1 : (x->img > y->img) - (x->img < y->img);
}

static void half_dist(int h, uint16_t d) {
    static uint32_t cnt[65536];
    memset(cnt, 0, sizeof cnt);
    if (!d) cnt[0] = 65536;
    else for (uint32_t x = 0; x < 65536; x++) cnt[(uint16_t)(THT[h][x] ^ THT[h][x ^ d])]++;
    int n = 0;
    tent *e = malloc(65536 * sizeof(tent));
    uint32_t minnum = (uint32_t)ceil(ldexp(1.0, 16) * pow(2.0, MINLOG));
    if (minnum < 1) minnum = 1;
    for (uint32_t o = 0; o < 65536; o++) {
        if (!cnt[o] || cnt[o] < minnum) continue;
        uint32_t img = toy_perm(h ? o : o << 16);
        uint8_t m = nmask(img);
        if (__builtin_popcount(m) > K) continue;
        e[n].img = img; e[n].num = cnt[o]; e[n].mask = m; n++;
    }
    qsort(e, (size_t)n, sizeof *e, cmp_tent);
    HE[h][d] = e; HNE[h][d] = n;
}

int main(int argc, char **argv) {
    int DMAX = 20, TOP = 64;
    const char *pre = "toy";
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "-k")) K = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-m")) MINLOG = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "-d")) DMAX = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-o")) pre = argv[i + 1];
        else if (!strcmp(argv[i], "-t")) TOP = atoi(argv[i + 1]);
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    toy_init(1, 4, 7, 3, 1, 2);
    HIN = malloc(65536 * sizeof(uint16_t));
    NHIN = 0;
    for (uint32_t v = 0; v < 65536; v++) if (nact16((uint16_t)v) <= K) HIN[NHIN++] = (uint16_t)v;
    for (int h = 0; h < 2; h++) for (int i = 0; i < NHIN; i++) half_dist(h, HIN[i]);
    /* nodes */
    size_t NN = 0, cap = 1 << 18;
    uint64_t *NODES = malloc(cap * sizeof(uint64_t));
    for (int i = 0; i < NHIN; i++)
        for (int j = 0; j < NHIN; j++) {
            if (nact16(HIN[i]) + nact16(HIN[j]) > K || (!HIN[i] && !HIN[j])) continue;
            if (NN == cap) { cap *= 2; NODES = realloc(NODES, cap * sizeof(uint64_t)); }
            NODES[NN++] = ((uint64_t)HIN[i] << 16) | HIN[j];
        }
    /* index by a dense 2^32 -> node map is too large; use a hash */
    size_t HC = 1; while (HC < 2 * NN) HC <<= 1;
    int64_t *HT = malloc(HC * sizeof(int64_t));
    memset(HT, 0xff, HC * sizeof(int64_t));
    for (size_t n = 0; n < NN; n++) {
        size_t hh = (size_t)((NODES[n] * 0x9E3779B97F4A7C15ULL) >> 20) & (HC - 1);
        while (HT[hh] >= 0) hh = (hh + 1) & (HC - 1);
        HT[hh] = (int64_t)n;
    }
#define NODE_GET(d, out) do { size_t hh_ = (size_t)(((uint64_t)(d) * 0x9E3779B97F4A7C15ULL) >> 20) & (HC - 1); (out) = -1; \
        while (HT[hh_] >= 0) { if (NODES[HT[hh_]] == (uint64_t)(d)) { (out) = HT[hh_]; break; } hh_ = (hh_ + 1) & (HC - 1); } } while (0)
    size_t ne = 0, ecap = 1 << 20;
    g_edge *E = malloc(ecap * sizeof *E);
    const double thr = pow(2.0, MINLOG);
    for (size_t x = 0; x < NN; x++) {
        uint16_t dl = (uint16_t)(NODES[x] >> 16), dr = (uint16_t)NODES[x];
        const tent *L = HE[0][dl], *R = HE[1][dr];
        for (int a = 0; a < HNE[0][dl]; a++) {
            double pl = L[a].num / 65536.0;
            if (pl < thr) break;
            for (int b = 0; b < HNE[1][dr]; b++) {
                double p = pl * (R[b].num / 65536.0);
                if (p < thr) break;
                if (__builtin_popcount(L[a].mask | R[b].mask) > K) continue;
                uint32_t o = L[a].img ^ R[b].img;
                int64_t t;
                NODE_GET(o, t);
                if (t < 0) continue;
                if (ne == ecap) { ecap *= 2; E = realloc(E, ecap * sizeof *E); }
                E[ne].from = (int32_t)x; E[ne].to = (int32_t)t; E[ne].p = p; ne++;
            }
        }
    }
    printf("Dipper-32 rot=1,4,7,3 perm pc=1 ps=2; K=%d min_log2p=%.1f nodes=%zu edges=%zu\n", K, MINLOG, NN, ne);
    g_min_mean_cycle(NN, E, ne, NULL, NODES, "toy");
    g_spectral(NN, E, ne, NULL, NODES, "toy", 6);
    fflush(stdout);
    /* propagation: 16 start vectors (0x8000 and v << 12, v = 1..15; 0x8000 is v = 8) */
    double *P[16], *W = malloc(NN * sizeof(double));
    for (int v = 1; v < 16; v++) {
        P[v] = calloc(NN, sizeof(double));
        int64_t s;
        NODE_GET((uint32_t)v << 12, s);
        P[v][s] = 1.0;
    }
    double *sig = malloc(NN * sizeof(double));
    int32_t *ord = malloc(NN * sizeof(int32_t));
    int64_t a8;
    NODE_GET(0x8000, a8);
    for (int d = 1; d <= DMAX; d++) {
        for (int v = 1; v < 16; v++) {
            memset(W, 0, NN * sizeof(double));
            for (size_t e = 0; e < ne; e++) W[E[e].to] += P[v][E[e].from] * E[e].p;
            memcpy(P[v], W, NN * sizeof(double));
        }
        for (size_t n = 0; n < NN; n++) { sig[n] = 0; for (int v = 1; v < 16; v++) sig[n] += P[v][n]; }
        /* top outputs by sig (partial selection) */
        int ntop = 0;
        for (size_t n = 0; n < NN; n++) {
            if (sig[n] <= 0) continue;
            if (ntop < TOP) ord[ntop++] = (int32_t)n;
            else if (sig[n] > sig[ord[TOP - 1]]) ord[TOP - 1] = (int32_t)n;
            else continue;
            for (int q = ntop - 1; q > 0 && sig[ord[q]] > sig[ord[q - 1]]; q--) { int32_t t = ord[q]; ord[q] = ord[q - 1]; ord[q - 1] = t; }
        }
        char path[512];
        snprintf(path, sizeof path, "%s_d%d.txt", pre, d);
        FILE *f = fopen(path, "w");
        for (int q = 0; q < ntop; q++)
            fprintf(f, "%08" PRIx64 " %.10e %.10e\n", NODES[ord[q]], sig[ord[q]], P[8][ord[q]]);
        fclose(f);
        printf("d=%d Pr[0x8000->0x8000] >= 2^%.4f, top output %08" PRIx64 " sig 2^%.4f\n", d, log2(P[8][a8]),
               NODES[ord[0]], log2(sig[ord[0]]));
        fflush(stdout);
    }
    return 0;
}
