/* graphalg.h: algorithms on a weighted difference graph (edges with probabilities):
 * minimum mean weight cycle (Howard policy iteration, as in lowgraph.c), spectral
 * radius (power iteration) and r-round propagation of a start vector.
 * Weights are w = -log2 p. */
#ifndef GRAPHALG_H
#define GRAPHALG_H

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int32_t from, to; double p; } g_edge;

/* nodes that lie on or lead to a cycle (alive[i] = 1); returns their number */
static size_t g_alive(size_t n, const g_edge *E, size_t ne, const uint8_t *keep, int *alive) {
    int *outdeg = calloc(n, sizeof(int));
    for (size_t i = 0; i < n; i++) alive[i] = keep ? keep[i] : 1;
    for (int changed = 1; changed;) {
        changed = 0;
        memset(outdeg, 0, n * sizeof(int));
        for (size_t e = 0; e < ne; e++)
            if (alive[E[e].to] && alive[E[e].from]) outdeg[E[e].from]++;
        for (size_t i = 0; i < n; i++)
            if (alive[i] && !outdeg[i]) { alive[i] = 0; changed = 1; }
    }
    free(outdeg);
    size_t c = 0;
    for (size_t i = 0; i < n; i++) c += alive[i];
    return c;
}

/* minimum mean weight cycle on the nodes with keep[i] (NULL = all); prints it.
 * Returns the mean weight (INFINITY if there is no cycle). */
static double g_min_mean_cycle(size_t n, const g_edge *E, size_t ne, const uint8_t *keep,
                               const uint64_t *NODES, const char *label) {
    int *alive = malloc(n * sizeof(int));
    size_t na = g_alive(n, E, ne, keep, alive);
    printf("[%s] nodes on or leading to cycles: %zu\n", label, na);
    if (!na) { printf("[%s] no cycle\n", label); free(alive); return INFINITY; }
    int32_t *pol = malloc(n * sizeof(int32_t));
    double *pw = malloc(n * sizeof(double)), *eta = malloc(n * sizeof(double)), *x = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) { pol[i] = -1; pw[i] = INFINITY; eta[i] = INFINITY; x[i] = 0; }
    for (size_t e = 0; e < ne; e++) {
        if (!alive[E[e].from] || !alive[E[e].to]) continue;
        double w = -log2(E[e].p);
        if (w < pw[E[e].from]) { pw[E[e].from] = w; pol[E[e].from] = E[e].to; }
    }
    int *mark = malloc(n * sizeof(int));
    int32_t *stk = malloc(n * sizeof(int32_t)), *cyc = malloc(n * sizeof(int32_t));
    const double EPS = 1e-9;
    int conv = 0;
    for (int iter = 0; iter < 100000; iter++) {
        memset(mark, 0, n * sizeof(int));
        for (size_t s0 = 0; s0 < n; s0++) {
            if (!alive[s0] || mark[s0]) continue;
            size_t sp = 0;
            int32_t u = (int32_t)s0;
            while (!mark[u]) { mark[u] = 1; stk[sp++] = u; u = pol[u]; }
            if (mark[u] == 1) {
                double sum = 0; size_t len = 0; int32_t c = u;
                do { sum += pw[c]; len++; c = pol[c]; } while (c != u);
                double m = sum / (double)len;
                x[u] = 0; eta[u] = m;
                size_t cl = 0; c = pol[u];
                while (c != u) { cyc[cl++] = c; c = pol[c]; }
                for (size_t q = cl; q-- > 0;) { eta[cyc[q]] = m; x[cyc[q]] = pw[cyc[q]] - m + x[pol[cyc[q]]]; mark[cyc[q]] = 2; }
                mark[u] = 2;
            }
            while (sp > 0) {
                int32_t t = stk[--sp];
                if (mark[t] == 2) continue;
                eta[t] = eta[pol[t]]; x[t] = pw[t] - eta[t] + x[pol[t]]; mark[t] = 2;
            }
        }
        int ch = 0;
        for (size_t e = 0; e < ne; e++) {
            int32_t a = E[e].from, b = E[e].to;
            if (!alive[a] || !alive[b]) continue;
            if (eta[b] < eta[a] - EPS) { pol[a] = b; pw[a] = -log2(E[e].p); eta[a] = eta[b]; ch = 1; }
        }
        if (!ch)
            for (size_t e = 0; e < ne; e++) {
                int32_t a = E[e].from, b = E[e].to;
                if (!alive[a] || !alive[b] || fabs(eta[b] - eta[a]) > EPS) continue;
                double w = -log2(E[e].p);
                if (w - eta[a] + x[b] < x[a] - EPS) { pol[a] = b; pw[a] = w; x[a] = w - eta[a] + x[b]; ch = 1; }
            }
        if (!ch) { printf("[%s] Howard converged after %d iterations\n", label, iter + 1); conv = 1; break; }
    }
    if (!conv) printf("[%s] Howard did NOT converge\n", label);
    double best = INFINITY; int32_t bv = -1;
    for (size_t i = 0; i < n; i++) if (alive[i] && eta[i] < best) { best = eta[i]; bv = (int32_t)i; }
    if (bv >= 0) {
        int32_t u = bv;
        for (size_t q = 0; q < n && q < 10000000; q++) u = pol[u];
        printf("[%s] minimum mean cycle: %.4f bits per round, cycle:", label, best);
        int32_t c = u; int len = 0; double sum = 0;
        do { printf(" %016" PRIx64, NODES[c]); sum += pw[c]; len++; c = pol[c]; } while (c != u && len < 64);
        printf("  (length %d, product 2^-%.4f)\n", len, sum);
    }
    free(alive); free(pol); free(pw); free(eta); free(x); free(mark); free(stk); free(cyc);
    return best;
}

/* log2 of the spectral radius of the matrix restricted to keep (NULL = all), by power iteration */
static double g_spectral(size_t n, const g_edge *E, size_t ne, const uint8_t *keep, const uint64_t *NODES,
                         const char *label, int show) {
    double *v = malloc(n * sizeof(double)), *w = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) v[i] = keep && !keep[i] ? 0.0 : 1.0;
    double lam = 0, prev = -1;
    for (int it = 0; it < 5000; it++) {
        memset(w, 0, n * sizeof(double));
        for (size_t e = 0; e < ne; e++)
            if (!keep || (keep[E[e].from] && keep[E[e].to])) w[E[e].to] += v[E[e].from] * E[e].p;
        double s = 0;
        for (size_t i = 0; i < n; i++) s = w[i] > s ? w[i] : s;
        if (s == 0) { lam = 0; break; }
        for (size_t i = 0; i < n; i++) v[i] = w[i] / s;
        lam = s;
        if (it > 200 && fabs(lam - prev) < 1e-15 * lam) break;
        prev = lam;
    }
    if (lam == 0) printf("[%s] no cycle: spectral radius 0\n", label);
    else printf("[%s] spectral radius: log2 rho = %.4f bits per round\n", label, log2(lam));
    for (int q = 0; q < show && lam > 0; q++) {
        size_t bi = 0;
        for (size_t i = 0; i < n; i++) if (v[i] > v[bi]) bi = i;
        printf("[%s]   %016" PRIx64 " %.4f\n", label, NODES[bi], v[bi]);
        v[bi] = -1;
    }
    free(v); free(w);
    return lam > 0 ? log2(lam) : -INFINITY;
}

/* r-round propagation from node a: prints log2 (M^r)_{a,b} for each target b, r = 1..R */
static void g_propagate(size_t n, const g_edge *E, size_t ne, int32_t a, const int32_t *tg, int ntg, int R,
                        const char *label) {
    double *v = calloc(n, sizeof(double)), *w = malloc(n * sizeof(double));
    v[a] = 1.0;
    for (int r = 1; r <= R; r++) {
        memset(w, 0, n * sizeof(double));
        for (size_t e = 0; e < ne; e++) w[E[e].to] += v[E[e].from] * E[e].p;
        double *t = v; v = w; w = t;
        printf("[%s] r=%d", label, r);
        double sum = 0;
        for (int q = 0; q < ntg; q++) { printf(" %.4f", v[tg[q]] > 0 ? log2(v[tg[q]]) : -INFINITY); sum += v[tg[q]]; }
        printf(" | set %.4f\n", sum > 0 ? log2(sum) : -INFINITY);
    }
    free(v); free(w);
}

#endif
