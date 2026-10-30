/* krsim: ranking of the 32 last-round key bits at data volumes that cannot be run.
 *
 *   krsim -f SURVIVORS -c KC -d KD -l LAM1 -w LAM0 [-n TRIALS] [-s SEED] [-t THREADS] [-T 8,9]
 *
 * SURVIVORS is a krec -o file for one key; KC, KD (hex) are its true guess. Each trial
 * draws K ~ Poisson(LAM1) pairs without replacement from the file (the pairs with
 * Delta x_{R-2} in the targets) and adds W ~ Poisson(LAM0) noise pairs: the wrong pairs
 * that pass the filter, z0 uniform and z1 = z0 ^ (delta << 27) with delta in 1..31.
 * All 2^32 guesses are then ranked as in krec, and the trial prints K, W, the score of
 * the true guess, and the numbers of guesses with a higher and with the same score. */
#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dipper.h"
#include "krank.h"

#define MAXS 4096

static int THREADS = 8;
static uint16_t TMASK = (1u << 8) | (1u << 9);
static kr_pair SP[MAXS];
static size_t NS;
static uint16_t *VT;
static _Atomic uint32_t NEXTKD;
typedef struct { uint64_t hist[MAXS + 1]; } rres;

static void *rank_worker(void *arg) {
    rres *o = arg;
    uint16_t s[32];
    for (;;) {
        uint32_t kd = atomic_fetch_add(&NEXTKD, 1);
        if (kd >= 65536) return NULL;
        for (uint32_t lo = 0; lo < 2048; lo++) {
            kr_score32(SP, VT, NS, (uint16_t)lo, (uint16_t)kd, s);
            for (int h = 0; h < 32; h++) o->hist[s[h]]++;
        }
    }
}

static int poisson(rng_t *g, double lam) {
    double l = exp(-lam), p = 1.0;
    int k = 0;
    do { k++; p *= (double)(rng_next(g) >> 11) * 0x1p-53; } while (p > l);
    return k - 1;
}

int main(int argc, char **argv) {
    const char *file = NULL;
    uint16_t kc = 0, kd = 0;
    double lam1 = 0, lam0 = 0;
    int trials = 1;
    uint64_t seed = 1;
    for (int i = 1; i + 1 < argc; i += 2) {
        const char *o = argv[i], *v = argv[i + 1];
        if (!strcmp(o, "-f")) file = v;
        else if (!strcmp(o, "-c")) kc = (uint16_t)strtoul(v, 0, 16);
        else if (!strcmp(o, "-d")) kd = (uint16_t)strtoul(v, 0, 16);
        else if (!strcmp(o, "-l")) lam1 = atof(v);
        else if (!strcmp(o, "-w")) lam0 = atof(v);
        else if (!strcmp(o, "-n")) trials = atoi(v);
        else if (!strcmp(o, "-s")) seed = strtoull(v, 0, 10);
        else if (!strcmp(o, "-t")) THREADS = atoi(v);
        else if (!strcmp(o, "-T")) {
            TMASK = 0;
            char *s = strdup(v), *tok = strtok(s, ",");
            while (tok) { TMASK |= (uint16_t)(1u << (strtoul(tok, 0, 16) & 15)); tok = strtok(0, ","); }
            free(s);
        } else { fprintf(stderr, "unknown option %s\n", o); return 1; }
    }
    if (!file || THREADS < 1 || THREADS > 256) {
        fprintf(stderr, "usage: krsim -f FILE -c KC -d KD -l LAM1 -w LAM0 [-n TRIALS] [-s SEED] [-t THREADS] [-T 8,9]\n");
        return 1;
    }
    dip_init();
    /* pool: the survivors of the file with a target difference */
    size_t cap = 4096, np = 0;
    kr_pair *pool = malloc(cap * sizeof *pool);
    FILE *f = fopen(file, "r");
    if (!f) { perror(file); return 1; }
    kr_pair q;
    while (fscanf(f, "%" SCNx64 " %" SCNx64 " %" SCNx64 " %" SCNx64 " %" SCNx64, &q.p0, &q.p1, &q.z0, &q.z1, &q.truth) == 5) {
        if (q.truth & ~0xF00000ULL) continue;
        int v = (int)(q.truth >> 20) & 15;
        if (!v || !((TMASK >> v) & 1)) continue;
        if (np == cap) { cap *= 2; pool = realloc(pool, cap * sizeof *pool); }
        pool[np++] = q;
    }
    fclose(f);
    /* every pooled pair must count for the true guess */
    for (size_t i = 0; i < np; i++)
        if (!kr_direct(&pool[i], kc, kd, TMASK)) { fprintf(stderr, "pool pair %zu does not count for the true guess\n", i); return 1; }
    printf("cmd=krsim file=%s pool=%zu lam1=%.4f lam0=%.4f trials=%d seed=%" PRIu64 "\n", file, np, lam1, lam0, trials, seed);
    rng_t g;
    rng_seed(&g, seed ^ 0x6b7273696d000000ULL);
    VT = malloc(MAXS * 1024 * sizeof(uint16_t));
    rres *RR = malloc((size_t)THREADS * sizeof *RR);
    pthread_t th[256];
    for (int it = 0; it < trials; it++) {
        int K = poisson(&g, lam1), W = poisson(&g, lam0);
        if ((size_t)K > np || K + W > MAXS) { fprintf(stderr, "trial too large: K=%d W=%d pool=%zu\n", K, W, np); return 1; }
        for (int i = 0; i < K; i++) { /* partial Fisher-Yates */
            size_t j = (size_t)i + (size_t)(rng_next(&g) % (np - (size_t)i));
            kr_pair t = pool[i]; pool[i] = pool[j]; pool[j] = t;
            SP[i] = pool[i];
        }
        for (int i = 0; i < W; i++) {
            kr_pair *o = &SP[K + i];
            uint64_t d = 1 + rng_next(&g) % 31;
            o->p0 = o->p1 = o->truth = 0;
            o->z0 = rng_next(&g);
            o->z1 = o->z0 ^ (d << 27);
        }
        NS = (size_t)(K + W);
        unsigned strue = 0;
        for (size_t i = 0; i < NS; i++) {
            kr_table(&SP[i], TMASK, VT + i * 1024);
            strue += (unsigned)kr_direct(&SP[i], kc, kd, TMASK);
        }
        memset(RR, 0, (size_t)THREADS * sizeof *RR);
        atomic_store(&NEXTKD, 0);
        for (int t = 0; t < THREADS; t++) pthread_create(&th[t], NULL, rank_worker, &RR[t]);
        for (int t = 0; t < THREADS; t++) pthread_join(th[t], NULL);
        uint64_t above = 0, tie = 0;
        for (int t = 0; t < THREADS; t++) {
            tie += RR[t].hist[strue];
            for (size_t s = strue + 1; s <= MAXS; s++) above += RR[t].hist[s];
        }
        printf("trial %d K=%d W=%d score=%u higher=%" PRIu64 " same=%" PRIu64 " bits=%.3f\n", it, K, W, strue, above, tie,
               32 - log2((double)(above + tie)));
        fflush(stdout);
    }
    return 0;
}
