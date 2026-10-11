/* lintool: exact one-round linear correlations of Dipper and a search for
 * iterative linear approximations among masks with few active nibbles.
 *
 * One round (without the key, which only changes signs) maps x to y = P(ARX(rot(S(x)))).
 * For an input mask a and output mask b, cor(a, b) = 2^-64 sum_x (-1)^(a.x + b.y).
 * The two halves (A,B) and (C,D) are independent, so cor = cor_L * cor_R, and for one half
 *   cor((aA,aB),(bA,bB)) = sum_{u,v} cadd(u,v -> bA) * cS(aA, rotr(u,rA)) * cS(aB, rotr(v^bB,rB))
 * where cS is the correlation of 4 parallel S-boxes and cadd the correlation of 16-bit
 * addition (computed exactly by a 2x2 carry transfer matrix).  Since cS(0, g) = 0 for g != 0,
 * u and v^bB are supported on the rotated active nibbles of aA and aB, so the sum is short.
 * Values are dyadic rationals and are exact in double precision at these sizes.
 *
 *   lintool corr -a A -b B [-R rots]     exact cor(A, B) by the formula
 *   lintool enum -a A -b B [-R rots]     exact cor(A, B) by enumeration of each half (2^32 each)
 *   lintool search -k MAXNIB -m MINLOG2 [-R rots] [-L MAXLEN]
 *       builds the graph of one-round approximations between masks with at most MAXNIB
 *       active nibbles and |cor| >= 2^MINLOG2, and prints the best cycles (iterative
 *       approximations) up to length MAXLEN.
 */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dipper.h"

static double LAT[16][16]; /* LAT[a][b] = correlation of a.x + b.S(x) */

static void lat_init(void) {
    for (int a = 0; a < 16; a++)
        for (int b = 0; b < 16; b++) {
            int s = 0;
            for (int x = 0; x < 16; x++) s += (__builtin_popcount((a & x) ^ (b & DIP_SBOX[x])) & 1) ? -1 : 1;
            LAT[a][b] = s / 16.0;
        }
}

/* correlation of w.(x+y) + u.x + v.y over uniform 16-bit x, y */
static double cadd(uint16_t u, uint16_t v, uint16_t w) {
    double st[2] = {1, 0};
    for (int i = 0; i < 16; i++) {
        int ui = (u >> i) & 1, vi = (v >> i) & 1, wi = (w >> i) & 1;
        double ns[2] = {0, 0};
        for (int c = 0; c < 2; c++) {
            if (st[c] == 0) continue;
            for (int x = 0; x < 2; x++)
                for (int y = 0; y < 2; y++) {
                    int s = x ^ y ^ c, nc = (x & y) | (x & c) | (y & c);
                    int e = (wi & s) ^ (ui & x) ^ (vi & y);
                    ns[nc] += st[c] * (e ? -0.25 : 0.25);
                }
        }
        st[0] = ns[0]; st[1] = ns[1];
    }
    return st[0] + st[1];
}

static double cS16(uint16_t a, uint16_t g) {
    double c = 1;
    for (int j = 0; j < 4; j++) {
        c *= LAT[(a >> (4 * j)) & 15][(g >> (4 * j)) & 15];
        if (c == 0) return 0;
    }
    return c;
}

static uint16_t nibmask(uint16_t a) {
    uint16_t m = 0;
    for (int j = 0; j < 4; j++) if ((a >> (4 * j)) & 15) m |= (uint16_t)(15u << (4 * j));
    return m;
}

/* enumerate all submasks of support m (as 16-bit) */
static int submasks(uint16_t m, uint16_t *out) {
    int n = 0;
    uint16_t s = 0;
    do { out[n++] = s; s = (uint16_t)((s - m) & m); } while (s);
    return n;
}

/* exact half correlation; wa = added word index (0 or 2), wb = passed word index */
static double half_corr(int wa, int wb, uint16_t aA, uint16_t aB, uint16_t bA, uint16_t bB) {
    int ra = DIP_ROT[wa], rb = DIP_ROT[wb];
    uint16_t mu = rotl16(nibmask(aA), ra), mv = rotl16(nibmask(aB), rb);
    static uint16_t U[65536], Vs[65536];
    int nu = submasks(mu, U), nv = submasks(mv, Vs);
    double tot = 0;
    for (int i = 0; i < nu; i++) {
        double s1 = cS16(aA, rotr16(U[i], ra));
        if (s1 == 0) continue;
        for (int k = 0; k < nv; k++) {
            uint16_t g = Vs[k]; /* g = v ^ bB, supported on mv */
            double s2 = cS16(aB, rotr16(g, rb));
            if (s2 == 0) continue;
            tot += s1 * s2 * cadd(U[i], (uint16_t)(g ^ bB), bA);
        }
    }
    return tot;
}

static double corr(uint64_t a, uint64_t b) {
    uint64_t bb = dip_perm_inv(b); /* mask on the ARX output */
    double cl = half_corr(0, 1, a >> 48, (a >> 32) & 0xFFFF, bb >> 48, (bb >> 32) & 0xFFFF);
    if (cl == 0) return 0;
    return cl * half_corr(2, 3, (a >> 16) & 0xFFFF, a & 0xFFFF, (bb >> 16) & 0xFFFF, bb & 0xFFFF);
}

static double half_enum(int wa, int wb, uint16_t aA, uint16_t aB, uint16_t bA, uint16_t bB) {
    int64_t s = 0;
    for (uint32_t B = 0; B < 65536; B++) {
        uint16_t v = DIP_TW[wb][B];
        int pb = __builtin_parity((aB & B) ^ (bB & v));
        for (uint32_t A = 0; A < 65536; A++) {
            uint16_t z = (uint16_t)(DIP_TW[wa][A] + v);
            s += (__builtin_parity((aA & A) ^ (bA & z)) ^ pb) ? -1 : 1;
        }
    }
    return ldexp((double)s, -32);
}
static double corr_enum(uint64_t a, uint64_t b) {
    uint64_t bb = dip_perm_inv(b);
    return half_enum(0, 1, a >> 48, (a >> 32) & 0xFFFF, bb >> 48, (bb >> 32) & 0xFFFF) *
           half_enum(2, 3, (a >> 16) & 0xFFFF, a & 0xFFFF, (bb >> 16) & 0xFFFF, bb & 0xFFFF);
}

/* ---- search ---- */
typedef struct { uint64_t m; } mask_t;
static uint64_t *MS; static int NM;
static int nnib(uint64_t m) { int n = 0; for (int j = 0; j < 16; j++) n += ((m >> (4 * j)) & 15) != 0; return n; }

static void gen_masks(int k) {
    /* all nonzero masks with 1..k active nibbles (k <= 2) */
    int cap = 16 * 15 + (k >= 2 ? 120 * 225 : 0);
    MS = malloc((size_t)cap * sizeof *MS);
    NM = 0;
    for (int p = 0; p < 16; p++)
        for (int x = 1; x < 16; x++) MS[NM++] = (uint64_t)x << (4 * p);
    if (k >= 2)
        for (int p = 0; p < 16; p++)
            for (int q = p + 1; q < 16; q++)
                for (int x = 1; x < 16; x++)
                    for (int y = 1; y < 16; y++) MS[NM++] = ((uint64_t)x << (4 * p)) | ((uint64_t)y << (4 * q));
}
static int find_mask(uint64_t m) {
    int lo = 0, hi = NM - 1;
    while (lo <= hi) { int mid = (lo + hi) / 2; if (MS[mid] == m) return mid; if (MS[mid] < m) lo = mid + 1; else hi = mid - 1; }
    return -1;
}
static int cmpu64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y); }

typedef struct { int to; double c; } ledge_t;
static ledge_t **ADJ; static int *NADJ;

static int cmd_search(int k, double minlog, int L) {
    gen_masks(k);
    qsort(MS, (size_t)NM, sizeof *MS, cmpu64);
    double thr = ldexp(1.0, (int)floor(minlog));
    ADJ = calloc((size_t)NM, sizeof *ADJ);
    NADJ = calloc((size_t)NM, sizeof *NADJ);
    int *cap = calloc((size_t)NM, sizeof *cap);
    long ne = 0;
    /* edges a -> b with a, b in the family: loop over a, then over b in the family.
     * Pruning: the output mask after P^-1 must touch only words reachable from a's words. */
    for (int ia = 0; ia < NM; ia++) {
        uint64_t a = MS[ia];
        for (int ib = 0; ib < NM; ib++) {
            uint64_t b = MS[ib], bb = dip_perm_inv(b);
            /* quick zero tests: a half with zero input mask needs zero output mask */
            if (!(a >> 32) && (bb >> 32)) continue;
            if (!(a & 0xFFFFFFFFULL) && (bb & 0xFFFFFFFFULL)) continue;
            if ((a >> 32) && !(bb >> 32)) continue;
            if ((a & 0xFFFFFFFFULL) && !(bb & 0xFFFFFFFFULL)) continue;
            double c = corr(a, b);
            if (fabs(c) < thr) continue;
            if (NADJ[ia] == cap[ia]) { cap[ia] = cap[ia] ? 2 * cap[ia] : 8; ADJ[ia] = realloc(ADJ[ia], (size_t)cap[ia] * sizeof(ledge_t)); }
            ADJ[ia][NADJ[ia]].to = ib; ADJ[ia][NADJ[ia]].c = c; NADJ[ia]++; ne++;
            if (getenv("LINTOOL_EDGES")) printf("edge a=%016" PRIx64 " b=%016" PRIx64 " cor=%.10g\n", a, b, c);
        }
    }
    fprintf(stderr, "masks=%d edges=%ld\n", NM, ne);
    /* best cycles up to length L by DFS from each node (node index must be minimal in cycle) */
    typedef struct { double w; int len; int nodes[8]; double c; } cyc_t;
    cyc_t best[20]; int nb = 0;
    int path[8]; double cpath[8];
    for (int s = 0; s < NM; s++) {
        /* iterative DFS */
        int depth = 0, idx[8] = {0};
        path[0] = s; cpath[0] = 1;
        while (depth >= 0) {
            int u = path[depth];
            if (idx[depth] >= NADJ[u]) { depth--; continue; }
            ledge_t e = ADJ[u][idx[depth]++];
            double c = cpath[depth] * e.c;
            if (e.to == s) {
                cyc_t cy; cy.len = depth + 1; cy.c = c; cy.w = log2(fabs(c)) / cy.len;
                for (int q = 0; q <= depth; q++) cy.nodes[q] = path[q];
                if (nb < 20 || cy.w > best[nb - 1].w) {
                    int pos = nb < 20 ? nb++ : 19;
                    best[pos] = cy;
                    for (int q = pos; q > 0 && best[q].w > best[q - 1].w; q--) { cyc_t t = best[q]; best[q] = best[q - 1]; best[q - 1] = t; }
                }
                continue;
            }
            if (e.to < s || depth + 1 >= L) continue;
            int seen = 0;
            for (int q = 0; q <= depth; q++) if (path[q] == e.to) seen = 1;
            if (seen) continue;
            depth++; path[depth] = e.to; cpath[depth] = c; idx[depth] = 0;
        }
    }
    for (int i = 0; i < nb; i++) {
        printf("cycle len=%d cor=%.10g log2abs_per_round=%.4f masks=", best[i].len, best[i].c, best[i].w);
        for (int q = 0; q < best[i].len; q++) printf("%s%016" PRIx64, q ? "," : "", MS[best[i].nodes[q]]);
        printf("\n");
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: lintool corr|enum|search ...\n"); return 1; }
    int rot[4] = {1, 4, 7, 11}, k = 1, L = 4;
    uint64_t a = 0, b = 0;
    double minlog = -6;
    for (int i = 2; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "-a")) a = strtoull(argv[i + 1], 0, 16);
        else if (!strcmp(argv[i], "-b")) b = strtoull(argv[i + 1], 0, 16);
        else if (!strcmp(argv[i], "-k")) k = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-m")) minlog = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "-L")) L = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-R")) sscanf(argv[i + 1], "%d,%d,%d,%d", &rot[0], &rot[1], &rot[2], &rot[3]);
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    if (L > 8) L = 8;
    dip_init_rot(rot[0], rot[1], rot[2], rot[3]);
    lat_init();
    if (!strcmp(argv[1], "corr")) { double c = corr(a, b); printf("a=%016" PRIx64 " b=%016" PRIx64 " cor=%.12g log2abs=%.4f\n", a, b, c, c ? log2(fabs(c)) : -INFINITY); return 0; }
    if (!strcmp(argv[1], "enum")) { double c = corr_enum(a, b); printf("a=%016" PRIx64 " b=%016" PRIx64 " cor=%.12g log2abs=%.4f\n", a, b, c, c ? log2(fabs(c)) : -INFINITY); return 0; }
    if (!strcmp(argv[1], "search")) return cmd_search(k, minlog, L);
    return 1;
}
