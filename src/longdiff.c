/* longdiff: long differential measurement with independent uniform round keys
 * per pair (the Markov model, as `dltool diff -k indep`), split into chunks so
 * that it can be interrupted and resumed.
 *
 *   longdiff -r R -i DIN -o DOUT[,DOUT..] -n LOG2PAIRS -c LOG2CHUNK -s SEED -f STATE
 *            [-t THREADS] [-l LOG] [-p SECONDS]
 *   longdiff bench -r R [-t THREADS] [-n LOG2PAIRS]
 *
 * Chunk k (2^LOG2CHUNK pairs) draws its plaintexts and round keys from its own
 * generator seeded by (SEED, k), so the counts do not depend on the number of
 * threads, on the order in which chunks finish, or on interruptions.
 * Every finished chunk is appended to STATE ("chunk k hits... set seconds")
 * and flushed to disk; a restart with the same arguments skips those chunks.
 * Progress (chunks, hits so far, rate, ETA) goes to LOG (default stderr) every
 * SECONDS (default 60). The final totals are printed in the format of dltool.
 */
#define _DARWIN_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifndef BS_LANES
#define BS_LANES 4
#endif
#include "dipper.h"
#include "bitslice.h"
#include "bs_const.h"

#define MAXOUT 16
#define FORMAT 1

typedef struct {
    int rounds, nout, log2n, log2c, threads, period;
    uint64_t din, dout[MAXOUT], seed;
    const char *state, *log;
} cfg_t;

typedef struct { uint64_t hits[MAXOUT], set; double sec; } chunk_t;

static cfg_t C;
static uint64_t NCHUNK;
static uint8_t *DONE;          /* chunk finished (from the state file or this run) */
static chunk_t *RES;
static _Atomic uint64_t NEXT;   /* next chunk index to try */
static _Atomic uint64_t NDONE, NDONE_RUN;
static pthread_mutex_t LOCK = PTHREAD_MUTEX_INITIALIZER;
static FILE *STATE, *LOG;

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

static uint64_t chunk_seed(uint64_t seed, uint64_t k) {
    uint64_t x = seed * 0x9E3779B97F4A7C15ULL ^ (k + 1) * 0xD1B54A32D192ED03ULL;
    return splitmix(&x);
}

/* one round on both states of 256*BS_LANES/4 pairs, same round key: bs_round with the
 * index maps as compile-time constants, writing to separate output arrays, so the
 * rotations and the permutation cost nothing. Checked against bs_round in main. */
static inline __attribute__((always_inline)) void round2(const V *restrict X, const V *restrict Y,
                                                          V *restrict X2, V *restrict Y2, const V *restrict K) {
    V SX[64], SY[64];
#pragma clang loop unroll(full)
    for (int j = 0; j < 16; j++) {
        V a0 = X[4 * j] ^ K[4 * j], a1 = X[4 * j + 1] ^ K[4 * j + 1], a2 = X[4 * j + 2] ^ K[4 * j + 2],
          a3 = X[4 * j + 3] ^ K[4 * j + 3];
        V b0 = Y[4 * j] ^ K[4 * j], b1 = Y[4 * j + 1] ^ K[4 * j + 1], b2 = Y[4 * j + 2] ^ K[4 * j + 2],
          b3 = Y[4 * j + 3] ^ K[4 * j + 3];
        bs_sbox(&a0, &a1, &a2, &a3);
        bs_sbox(&b0, &b1, &b2, &b3);
        SX[4 * j] = a0; SX[4 * j + 1] = a1; SX[4 * j + 2] = a2; SX[4 * j + 3] = a3;
        SY[4 * j] = b0; SY[4 * j + 1] = b1; SY[4 * j + 2] = b2; SY[4 * j + 3] = b3;
    }
#pragma clang loop unroll(full)
    for (int h = 0; h < 2; h++) {
        const int a = h ? 16 : 48, b = h ? 0 : 32;
        V cx = vbcast(0), cy = vbcast(0);
#pragma clang loop unroll(full)
        for (int k = 0; k < 16; k++) {
            V x = SX[BSC_ROTSRC[a + k]], y = SX[BSC_ROTSRC[b + k]], t = x ^ y;
            X2[BSC_PDST[a + k]] = t ^ cx;
            X2[BSC_PDST[b + k]] = y;
            cx = (x & y) | (cx & t);
            V u = SY[BSC_ROTSRC[a + k]], v = SY[BSC_ROTSRC[b + k]], w = u ^ v;
            Y2[BSC_PDST[a + k]] = w ^ cy;
            Y2[BSC_PDST[b + k]] = v;
            cy = (u & v) | (cy & w);
        }
    }
}

/* 2^log2c pairs (X, X ^ din) through R rounds with fresh round keys per pair */
static void run_chunk(uint64_t k, chunk_t *out) {
    const int R = C.rounds;
    V X[64], Y[64], X2[64], Y2[64], K[64], DI[64];
    vrng_t r;
    vrng_seed(&r, chunk_seed(C.seed, k));
    for (int i = 0; i < 64; i++) DI[i] = vbcast(((C.din >> i) & 1) ? ~0ULL : 0);
    memset(out, 0, sizeof *out);
    uint64_t batches = (1ULL << C.log2c) / BS_BLOCKS;
    for (uint64_t b = 0; b < batches; b++) {
        for (int i = 0; i < 64; i++) {
            X[i] = vrng_next(&r);
            Y[i] = X[i] ^ DI[i];
        }
        for (int q = 0; q < R; q += 2) {
            for (int i = 0; i < 64; i++) K[i] = vrng_next(&r);
            round2(X, Y, X2, Y2, K);
            if (q + 1 == R) {
                memcpy(X, X2, sizeof X);
                memcpy(Y, Y2, sizeof Y);
                break;
            }
            for (int i = 0; i < 64; i++) K[i] = vrng_next(&r);
            round2(X2, Y2, X, Y, K);
        }
        V any = vbcast(0);
        for (int o = 0; o < C.nout; o++) {
            V m = vbcast(~0ULL);
            for (int i = 0; i < 64; i++) m &= ((C.dout[o] >> i) & 1) ? X[i] ^ Y[i] : ~(X[i] ^ Y[i]);
            out->hits[o] += vpopcnt(m);
            any |= m;
        }
        out->set += vpopcnt(any);
    }
}

static void *worker(void *arg) {
    (void)arg;
    for (;;) {
        uint64_t k = atomic_fetch_add(&NEXT, 1);
        if (k >= NCHUNK) return NULL;
        if (DONE[k]) continue;
        chunk_t res;
        double t0 = now();
        run_chunk(k, &res);
        res.sec = now() - t0;
        pthread_mutex_lock(&LOCK);
        RES[k] = res;
        DONE[k] = 1;
        fprintf(STATE, "chunk %" PRIu64, k);
        for (int o = 0; o < C.nout; o++) fprintf(STATE, " %" PRIu64, res.hits[o]);
        fprintf(STATE, " %" PRIu64 " %.3f\n", res.set, res.sec);
        fflush(STATE);
        fsync(fileno(STATE));
        pthread_mutex_unlock(&LOCK);
        atomic_fetch_add(&NDONE, 1);
        atomic_fetch_add(&NDONE_RUN, 1);
    }
}

static void header(char *buf, size_t n) {
    int p = snprintf(buf, n, "# longdiff format=%d rounds=%d din=%016" PRIx64 " log2pairs=%d log2chunk=%d seed=%" PRIu64
                     " lanes=%d douts=", FORMAT, C.rounds, C.din, C.log2n, C.log2c, C.seed, BS_LANES);
    for (int o = 0; o < C.nout; o++) p += snprintf(buf + p, n - (size_t)p, "%s%016" PRIx64, o ? "," : "", C.dout[o]);
}

/* read the state file if it exists; returns the number of chunks recovered */
static uint64_t load_state(void) {
    char want[1024], line[1024];
    header(want, sizeof want);
    FILE *f = fopen(C.state, "r");
    if (!f) {
        STATE = fopen(C.state, "w");
        if (!STATE) { perror(C.state); exit(1); }
        fprintf(STATE, "%s\n", want);
        fflush(STATE);
        return 0;
    }
    if (!fgets(line, sizeof line, f) || strncmp(line, want, strlen(want)) || line[strlen(want)] != '\n') {
        fprintf(stderr, "%s: header does not match these arguments\n  file: %s  want: %s\n", C.state, line, want);
        exit(1);
    }
    uint64_t n = 0;
    long good = ftell(f);
    while (fgets(line, sizeof line, f)) {
        size_t L = strlen(line);
        if (L == 0 || line[L - 1] != '\n') break; /* torn last line: dropped and recomputed */
        chunk_t c;
        memset(&c, 0, sizeof c);
        uint64_t k;
        char *s = line;
        int used;
        if (sscanf(s, "chunk %" SCNu64 "%n", &k, &used) != 1 || k >= NCHUNK) break;
        s += used;
        int ok = 1;
        for (int o = 0; o < C.nout && ok; o++, s += used) ok = sscanf(s, "%" SCNu64 "%n", &c.hits[o], &used) == 1;
        if (!ok || sscanf(s, "%" SCNu64 " %lf", &c.set, &c.sec) != 2) break;
        if (!DONE[k]) n++;
        DONE[k] = 1;
        RES[k] = c;
        good = ftell(f);
    }
    fclose(f);
    if (truncate(C.state, good)) { perror("truncate"); exit(1); }
    STATE = fopen(C.state, "a");
    if (!STATE) { perror(C.state); exit(1); }
    return n;
}

static void totals(uint64_t *hits, uint64_t *set, double *sec, uint64_t *nd) {
    memset(hits, 0, sizeof(uint64_t) * MAXOUT);
    *set = 0; *sec = 0; *nd = 0;
    for (uint64_t k = 0; k < NCHUNK; k++)
        if (DONE[k]) {
            for (int o = 0; o < C.nout; o++) hits[o] += RES[k].hits[o];
            *set += RES[k].set; *sec += RES[k].sec; (*nd)++;
        }
}

static void progress(double t0) {
    uint64_t hits[MAXOUT], set, nd;
    double sec;
    pthread_mutex_lock(&LOCK);
    totals(hits, &set, &sec, &nd);
    pthread_mutex_unlock(&LOCK);
    uint64_t run = atomic_load(&NDONE_RUN);
    double el = now() - t0, rate = run ? (double)run * (double)(1ULL << C.log2c) / el : 0;
    double eta = run ? el / (double)run * (double)(NCHUNK - nd) : 0;
    time_t w = time(NULL);
    char ts[32];
    strftime(ts, sizeof ts, "%Y-%m-%dT%H:%M:%SZ", gmtime(&w));
    fprintf(LOG, "%s chunks %" PRIu64 "/%" PRIu64 " (%.1f%%) pairs 2^%.3f hits", ts, nd, NCHUNK,
            100.0 * (double)nd / (double)NCHUNK, nd ? log2((double)nd) + C.log2c : 0.0);
    for (int o = 0; o < C.nout; o++) fprintf(LOG, "%s%" PRIu64, o ? "," : " ", hits[o]);
    fprintf(LOG, " set %" PRIu64 " rate 2^%.2f pairs/s ETA %.0f min\n", set, rate > 0 ? log2(rate) : 0.0, eta / 60);
    fflush(LOG);
}

/* constant index maps equal bs_init's; round2 equals bs_round on random states */
static int selfcheck(void) {
    int bad = 0;
    for (int i = 0; i < 4; i++) bad += BSC_ROT[i] != DIP_ROT[i];
    for (int i = 0; i < 64; i++) bad += (BSC_ROTSRC[i] != BS_ROTSRC[i]) + (BSC_PDST[i] != BS_PDST[i]);
    vrng_t r;
    vrng_seed(&r, 12345);
    for (int n = 0; n < 100; n++) {
        V X[64], Y[64], X2[64], Y2[64], K[64];
        for (int i = 0; i < 64; i++) { X[i] = vrng_next(&r); Y[i] = vrng_next(&r); K[i] = vrng_next(&r); }
        round2(X, Y, X2, Y2, K);
        bs_round(X, K);
        bs_round(Y, K);
        for (int i = 0; i < 64; i++)
            for (int l = 0; l < BS_LANES; l++) bad += (X[i][l] != X2[i][l]) + (Y[i][l] != Y2[i][l]);
    }
    if (bad) fprintf(stderr, "selfcheck FAILED (%d)\n", bad);
    return bad != 0;
}

static int bench(void) {
    /* throughput of run_chunk on THREADS threads, chunks of 2^log2c pairs */
    NCHUNK = 1ULL << (C.log2n - C.log2c);
    DONE = calloc(NCHUNK, 1);
    RES = calloc(NCHUNK, sizeof(chunk_t));
    STATE = fopen("/dev/null", "w");
    pthread_t th[256];
    double t0 = now();
    for (int t = 0; t < C.threads; t++) pthread_create(&th[t], NULL, worker, NULL);
    for (int t = 0; t < C.threads; t++) pthread_join(th[t], NULL);
    double el = now() - t0, rate = (double)(1ULL << C.log2n) / el;
    printf("bench rounds=%d threads=%d lanes=%d pairs=2^%d seconds=%.2f rate=2^%.3f pairs/s 2^40_pairs_hours=%.2f\n",
           C.rounds, C.threads, BS_LANES, C.log2n, el, log2(rate), ldexp(1.0, 40) / rate / 3600);
    return 0;
}

int main(int argc, char **argv) {
    int is_bench = argc > 1 && !strcmp(argv[1], "bench");
    memset(&C, 0, sizeof C);
    C.threads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    C.period = 60;
    C.log2c = 28;
    C.seed = 1;
    C.log2n = 30;
    for (int i = 1 + is_bench; i < argc; i += 2) {
        const char *o = argv[i], *v = i + 1 < argc ? argv[i + 1] : "";
        if (!strcmp(o, "-r")) C.rounds = atoi(v);
        else if (!strcmp(o, "-i")) C.din = strtoull(v, 0, 16);
        else if (!strcmp(o, "-o")) {
            char *s = strdup(v), *tok = strtok(s, ",");
            while (tok && C.nout < MAXOUT) { C.dout[C.nout++] = strtoull(tok, 0, 16); tok = strtok(0, ","); }
            free(s);
        } else if (!strcmp(o, "-n")) C.log2n = atoi(v);
        else if (!strcmp(o, "-c")) C.log2c = atoi(v);
        else if (!strcmp(o, "-s")) C.seed = strtoull(v, 0, 10);
        else if (!strcmp(o, "-t")) C.threads = atoi(v);
        else if (!strcmp(o, "-f")) C.state = v;
        else if (!strcmp(o, "-l")) C.log = v;
        else if (!strcmp(o, "-p")) C.period = atoi(v);
        else { fprintf(stderr, "unknown option %s\n", o); return 1; }
    }
    if (C.threads < 1 || C.threads > 256) C.threads = 4;
    if (is_bench && C.log2c > C.log2n - 4) C.log2c = C.log2n - 4 > 9 ? C.log2n - 4 : 9;
    if (C.rounds < 1 || C.rounds > 64 || C.log2c < 9 || C.log2c > C.log2n || C.log2n > 62 ||
        (!is_bench && (C.nout < 1 || !C.state))) {
        fprintf(stderr, "usage: longdiff -r R -i DIN -o DOUT[,..] -n LOG2PAIRS -c LOG2CHUNK -s SEED -f STATE "
                        "[-t THREADS] [-l LOG] [-p SECONDS]\n       longdiff bench -r R [-t THREADS] [-n LOG2PAIRS]\n");
        return 1;
    }
    dip_init();
    bs_init();
    if (selfcheck()) return 1;
    if (is_bench) return bench();

    LOG = stderr;
    if (C.log && !(LOG = fopen(C.log, "a"))) { perror(C.log); return 1; }
    NCHUNK = 1ULL << (C.log2n - C.log2c);
    DONE = calloc(NCHUNK, 1);
    RES = calloc(NCHUNK, sizeof(chunk_t));
    uint64_t have = load_state();
    atomic_store(&NDONE, have);
    char hdr[1024];
    header(hdr, sizeof hdr);
    fprintf(LOG, "%s\nstart: %d threads, %" PRIu64 "/%" PRIu64 " chunks already in %s\n", hdr, C.threads, have,
            NCHUNK, C.state);
    fflush(LOG);

    double t0 = now();
    pthread_t th[256];
    for (int t = 0; t < C.threads; t++) pthread_create(&th[t], NULL, worker, NULL);
    double last = t0;
    while (atomic_load(&NDONE) < NCHUNK) {
        usleep(200000);
        if (now() - last >= C.period) { progress(t0); last = now(); }
    }
    for (int t = 0; t < C.threads; t++) pthread_join(th[t], NULL);
    progress(t0);
    fclose(STATE);

    uint64_t hits[MAXOUT], set, nd;
    double sec;
    totals(hits, &set, &sec, &nd);
    uint64_t used = NCHUNK << C.log2c;
    printf("cmd=longdiff rounds=%d mode=indep seed=%" PRIu64 " rot=%d,%d,%d,%d chunks=%" PRIu64 " log2chunk=%d\n",
           C.rounds, C.seed, DIP_ROT[0], DIP_ROT[1], DIP_ROT[2], DIP_ROT[3], NCHUNK, C.log2c);
    printf("din=%016" PRIx64 " pairs=%" PRIu64 " log2pairs=%.3f\n", C.din, used, log2((double)used));
    for (int o = 0; o < C.nout; o++)
        printf("dout=%016" PRIx64 " hits=%" PRIu64 " log2p=%.4f\n", C.dout[o], hits[o],
               hits[o] ? log2((double)hits[o] / (double)used) : -INFINITY);
    printf("set hits=%" PRIu64 " log2p=%.4f\n", set, set ? log2((double)set / (double)used) : -INFINITY);
    printf("thread_seconds=%.0f\n", sec);
    return 0;
}
