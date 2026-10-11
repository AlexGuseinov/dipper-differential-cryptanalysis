/* dltool: measurement engine for the differential, differential-linear and
 * linear experiments on Dipper.
 *
 *   dltool test
 *   dltool enc   < lines "BITS KEYHEX PTHEX"    (prints ciphertext per line, 28 rounds)
 *   dltool diff  -r R -i DIN -o DOUT[,DOUT..] -n LOG2PAIRS [-k MODE] [-K NKEYS] [-s SEED] [-t THREADS]
 *   dltool dl    -r R -i DIN -n LOG2PAIRS [-k MODE] [-K NKEYS] [-s SEED] [-t THREADS]
 *   dltool lin   -r R -a INMASK -b OUTMASK -n LOG2SAMPLES -K NKEYS [-k MODE] [-s SEED] [-t THREADS]
 *
 * MODE: indep  fresh independent uniform round keys for every pair (Markov model)
 *       rk     independent uniform round keys, fixed per key (NKEYS keys)
 *       ks96   Dipper-64/96 key schedule, NKEYS random master keys
 *       ks128  Dipper-64/128 key schedule, NKEYS random master keys
 * Global option -R a,b,c,d replaces the rotation constants.
 * Output lines are "key=value" so that the step scripts can parse them.
 */
#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dipper.h"
#include "bitslice.h"

#define MAXOUT 64
#define MAXKEYS 4096

enum { M_INDEP, M_RK, M_KS96, M_KS128 };

typedef struct {
    int cmd, rounds, mode, nkeys, threads, nout;
    uint64_t din, dout[MAXOUT], amask, bmask, seed;
    int log2n;
} cfg_t;

typedef struct {
    const cfg_t *c;
    int tid;
    int key0, key1;          /* key range for this thread (rk/ks modes) */
    uint64_t pairs;          /* pairs for indep mode */
    uint64_t cnt[MAXOUT];    /* diff: hits per target */
    uint64_t setcnt;         /* diff: hits on the union of targets */
    uint64_t *keyset;        /* diff: union hits per key */
    uint64_t *keycnt;        /* diff: hits on target 0 per key (size nkeys) */
    uint64_t bitcnt[64];     /* dl: pairs with output-difference bit j set */
    double *keycorr;         /* lin: correlation per key */
    uint64_t done;
} job_t;

static uint64_t keyseed(const cfg_t *c, int key) { return c->seed * 1000003ULL + (uint64_t)key * 7919ULL + 17; }

static void make_rk(const cfg_t *c, int key, uint64_t *rk) {
    rng_t r;
    rng_seed(&r, keyseed(c, key));
    if (c->mode == M_RK) {
        for (int i = 0; i < c->rounds; i++) rk[i] = rng_next(&r);
        return;
    }
    dip_key K;
    memset(&K, 0, sizeof K);
    K.words = c->mode == M_KS96 ? 6 : 8;
    for (int w = 0; w < K.words; w++) K.k[w] = (uint16_t)rng_next(&r);
    dip_expand(K, c->rounds, rk);
}

static inline void record_bs(job_t *j, const V *D, int key) {
    const cfg_t *c = j->c;
    if (c->cmd == 1) {
        V any = vbcast(0);
        for (int o = 0; o < c->nout; o++) {
            V m = vbcast(~0ULL);
            uint64_t t = c->dout[o];
            for (int i = 0; i < 64; i++) m &= ((t >> i) & 1) ? D[i] : ~D[i];
            uint64_t h = vpopcnt(m);
            j->cnt[o] += h;
            any |= m;
            if (o == 0 && key >= 0) j->keycnt[key] += h;
        }
        uint64_t h = vpopcnt(any);
        j->setcnt += h;
        if (key >= 0) j->keyset[key] += h;
    } else {
        for (int i = 0; i < 64; i++) j->bitcnt[i] += vpopcnt(D[i]);
    }
}

static void *worker(void *arg) {
    job_t *j = arg;
    const cfg_t *c = j->c;
    int R = c->rounds;
    uint64_t rk[64];
    static __thread V KS[64][64];
    V X[64], Y[64], D[64], K[64], X0[64];
    vrng_t r;
    vrng_seed(&r, c->seed ^ (0xA5A5A5A5ULL + (uint64_t)j->tid * 0x1234567ULL));
    if (c->cmd == 1 || c->cmd == 2) {
        if (c->mode == M_INDEP) {
            uint64_t batches = j->pairs / BS_BLOCKS;
            for (uint64_t n = 0; n < batches; n++) {
                for (int i = 0; i < 64; i++) {
                    X[i] = vrng_next(&r);
                    Y[i] = X[i] ^ vbcast(((c->din >> i) & 1) ? ~0ULL : 0);
                }
                for (int q = 0; q < R; q++) {
                    for (int i = 0; i < 64; i++) K[i] = vrng_next(&r);
                    bs_round(X, K);
                    bs_round(Y, K);
                }
                for (int i = 0; i < 64; i++) D[i] = X[i] ^ Y[i];
                record_bs(j, D, -1);
            }
        } else {
            uint64_t per = (1ULL << c->log2n) / (uint64_t)c->nkeys / BS_BLOCKS;
            for (int key = j->key0; key < j->key1; key++) {
                make_rk(c, key, rk);
                for (int q = 0; q < R; q++) bs_key(rk[q], KS[q]);
                vrng_t rx;
                vrng_seed(&rx, keyseed(c, key) ^ 0x5DEECE66DULL);
                for (uint64_t n = 0; n < per; n++) {
                    for (int i = 0; i < 64; i++) {
                        X[i] = vrng_next(&rx);
                        Y[i] = X[i] ^ vbcast(((c->din >> i) & 1) ? ~0ULL : 0);
                    }
                    for (int q = 0; q < R; q++) {
                        bs_round(X, KS[q]);
                        bs_round(Y, KS[q]);
                    }
                    for (int i = 0; i < 64; i++) D[i] = X[i] ^ Y[i];
                    record_bs(j, D, key);
                }
            }
        }
    } else { /* linear */
        uint64_t per = (1ULL << c->log2n) / BS_BLOCKS;
        if (c->mode == M_INDEP) {
            fprintf(stderr, "lin needs a fixed-key mode\n");
            exit(1);
        }
        for (int key = j->key0; key < j->key1; key++) {
            make_rk(c, key, rk);
            for (int q = 0; q < R; q++) bs_key(rk[q], KS[q]);
            vrng_t rx;
            vrng_seed(&rx, keyseed(c, key) ^ 0xC0FFEEULL);
            int64_t s = 0;
            for (uint64_t n = 0; n < per; n++) {
                for (int i = 0; i < 64; i++) X0[i] = X[i] = vrng_next(&rx);
                for (int q = 0; q < R; q++) bs_round(X, KS[q]);
                V p = vbcast(0);
                for (int i = 0; i < 64; i++) {
                    if ((c->amask >> i) & 1) p ^= X0[i];
                    if ((c->bmask >> i) & 1) p ^= X[i];
                }
                s += BS_BLOCKS - 2 * (int64_t)vpopcnt(p);
            }
            j->keycorr[key] = (double)s / (double)(per * BS_BLOCKS);
        }
    }
    return NULL;
}

static int parse_mode(const char *s) {
    if (!strcmp(s, "indep")) return M_INDEP;
    if (!strcmp(s, "rk")) return M_RK;
    if (!strcmp(s, "ks96")) return M_KS96;
    if (!strcmp(s, "ks128")) return M_KS128;
    fprintf(stderr, "unknown mode %s\n", s);
    exit(1);
}
static const char *mode_name(int m) {
    return m == M_INDEP ? "indep" : m == M_RK ? "rk" : m == M_KS96 ? "ks96" : "ks128";
}

/* ---------------- self-test ---------------- */
static int run_test(void) {
    int fail = 0;
    struct { int bits; uint64_t p; uint16_t k[8]; uint64_t c; } tv[4] = {
        {96, 0x0123456789ABCDEFULL, {0xAABB, 0x8899, 0x6677, 0x4455, 0x2233, 0x0011}, 0x5F189202D3152B8FULL},
        {96, 0, {0}, 0xD873355EC63B1B4BULL},
        {128, 0x0123456789ABCDEFULL, {0x0E0F, 0x0C0D, 0x0A0B, 0x0809, 0x0607, 0x0405, 0x0203, 0x0001}, 0x4B46284387969060ULL},
        {128, 0xFFFFFFFFFFFFFFFFULL, {0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF}, 0x6ABB4518063706B0ULL},
    };
    for (int t = 0; t < 4; t++) {
        dip_key K;
        memset(&K, 0, sizeof K);
        K.words = tv[t].bits / 16;
        memcpy(K.k, tv[t].k, sizeof K.k);
        uint64_t rk[28];
        dip_expand(K, 28, rk);
        uint64_t c = dip_encrypt_rk(tv[t].p, rk, 28);
        int ok = c == tv[t].c;
        fail |= !ok;
        printf("TV%d Dipper-64/%d: got %016" PRIX64 " expected %016" PRIX64 " %s\n", t + 1, tv[t].bits, c,
               tv[t].c, ok ? "OK" : "FAIL");
    }
    rng_t r;
    rng_seed(&r, 1);
    int bad = 0;
    for (int n = 0; n < 1000000; n++) {
        uint64_t s = rng_next(&r), k = rng_next(&r);
        if (dip_round(s, k) != dip_round_ref(s, k)) bad++;
        uint64_t t = dip_round_nokey(s);
        if (dip_round_nokey_inv(t) != s) bad++;
    }
    printf("fast round vs bit-level reference and inverse, 10^6 random inputs: %d mismatches %s\n", bad,
           bad ? "FAIL" : "OK");
    fail |= bad != 0;
    {
        uint64_t blk[BS_BLOCKS], ref[BS_BLOCKS], rk2[28];
        V X[64], KS[64];
        int bb = 0;
        for (int q = 0; q < 28; q++) rk2[q] = rng_next(&r);
        for (int b = 0; b < BS_BLOCKS; b++) { blk[b] = rng_next(&r); ref[b] = dip_encrypt_rk(blk[b], rk2, 28); }
        bs_load(blk, X);
        for (int q = 0; q < 28; q++) { bs_key(rk2[q], KS); bs_round(X, KS); }
        bs_store(X, blk);
        for (int b = 0; b < BS_BLOCKS; b++) bb += blk[b] != ref[b];
        printf("bitsliced vs table implementation, 256 blocks x 28 rounds: %d mismatches %s\n", bb, bb ? "FAIL" : "OK");
        fail |= bb != 0;
    }
    printf("RC:");
    uint8_t l = 1;
    for (int i = 0; i < 28; i++) { printf(" %02X", l); l = dip_lfsr_step(l); }
    printf("\n");
    return fail;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: dltool test|diff|dl|lin [options]\n");
        return 1;
    }
    cfg_t c;
    memset(&c, 0, sizeof c);
    c.mode = M_INDEP; c.nkeys = 1; c.threads = 4; c.seed = 1; c.log2n = 20;
    int rot[4] = {1, 4, 7, 11};
    if (!strcmp(argv[1], "enc")) {
        dip_init();
        int bits; char kh[64], ph[32];
        while (scanf("%d %63s %31s", &bits, kh, ph) == 3) {
            dip_key K; memset(&K, 0, sizeof K); K.words = bits / 16;
            size_t L = strlen(kh);
            for (int w = 0; w < K.words; w++) {
                char buf[5] = {0};
                long st = (long)L - 4 * (w + 1);
                for (int q = 0; q < 4; q++) buf[q] = st + q >= 0 ? kh[st + q] : '0';
                K.k[w] = (uint16_t)strtoul(buf, 0, 16);
            }
            uint64_t rk[28];
            dip_expand(K, 28, rk);
            printf("%016" PRIX64 "\n", dip_encrypt_rk(strtoull(ph, 0, 16), rk, 28));
        }
        return 0;
    }
    if (!strcmp(argv[1], "test")) c.cmd = 0;
    else if (!strcmp(argv[1], "diff")) c.cmd = 1;
    else if (!strcmp(argv[1], "dl")) c.cmd = 2;
    else if (!strcmp(argv[1], "lin")) c.cmd = 3;
    else { fprintf(stderr, "unknown command\n"); return 1; }
    for (int i = 2; i < argc; i++) {
        const char *o = argv[i], *v = i + 1 < argc ? argv[i + 1] : "";
        if (!strcmp(o, "-r")) c.rounds = atoi(v);
        else if (!strcmp(o, "-i")) c.din = strtoull(v, 0, 16);
        else if (!strcmp(o, "-o")) {
            char *s = strdup(v), *tok = strtok(s, ",");
            while (tok && c.nout < MAXOUT) { c.dout[c.nout++] = strtoull(tok, 0, 16); tok = strtok(0, ","); }
            free(s);
        } else if (!strcmp(o, "-a")) c.amask = strtoull(v, 0, 16);
        else if (!strcmp(o, "-b")) c.bmask = strtoull(v, 0, 16);
        else if (!strcmp(o, "-n")) c.log2n = atoi(v);
        else if (!strcmp(o, "-k")) c.mode = parse_mode(v);
        else if (!strcmp(o, "-K")) c.nkeys = atoi(v);
        else if (!strcmp(o, "-s")) c.seed = strtoull(v, 0, 10);
        else if (!strcmp(o, "-t")) c.threads = atoi(v);
        else if (!strcmp(o, "-R")) sscanf(v, "%d,%d,%d,%d", &rot[0], &rot[1], &rot[2], &rot[3]);
        else { fprintf(stderr, "unknown option %s\n", o); return 1; }
        i++;
    }
    dip_init_rot(rot[0], rot[1], rot[2], rot[3]);
    bs_init();
    if (c.cmd == 0) return run_test();
    if (c.nkeys < 1 || c.nkeys > MAXKEYS) { fprintf(stderr, "bad -K\n"); return 1; }
    if (c.mode == M_INDEP) c.nkeys = 1;
    int T = c.threads;
    if (c.mode != M_INDEP && T > c.nkeys) T = c.nkeys;
    job_t *jobs = calloc((size_t)T, sizeof(job_t));
    pthread_t *th = calloc((size_t)T, sizeof(pthread_t));
    uint64_t *keycnt = calloc((size_t)c.nkeys, sizeof(uint64_t));
    double *keycorr = calloc((size_t)c.nkeys, sizeof(double));
    uint64_t *keyset = calloc((size_t)c.nkeys, sizeof(uint64_t));
    uint64_t total = 1ULL << c.log2n;
    for (int t = 0; t < T; t++) {
        jobs[t].c = &c; jobs[t].tid = t;
        jobs[t].pairs = total / (uint64_t)T + ((uint64_t)t < total % (uint64_t)T);
        jobs[t].key0 = c.nkeys * t / T; jobs[t].key1 = c.nkeys * (t + 1) / T;
        jobs[t].keycnt = keycnt; jobs[t].keycorr = keycorr; jobs[t].keyset = keyset;
        pthread_create(&th[t], NULL, worker, &jobs[t]);
    }
    uint64_t cnt[MAXOUT] = {0}, bit[64] = {0}, setc = 0;
    for (int t = 0; t < T; t++) {
        pthread_join(th[t], NULL);
        for (int o = 0; o < MAXOUT; o++) cnt[o] += jobs[t].cnt[o];
        for (int b = 0; b < 64; b++) bit[b] += jobs[t].bitcnt[b];
        setc += jobs[t].setcnt;
    }
    uint64_t used = 0;
    if (c.mode == M_INDEP)
        for (int t = 0; t < T; t++) used += jobs[t].pairs / BS_BLOCKS * BS_BLOCKS;
    else
        used = (total / (uint64_t)c.nkeys / BS_BLOCKS) * BS_BLOCKS * (uint64_t)c.nkeys;
    printf("cmd=%s rounds=%d mode=%s nkeys=%d seed=%" PRIu64 " rot=%d,%d,%d,%d\n", argv[1], c.rounds,
           mode_name(c.mode), c.nkeys, c.seed, rot[0], rot[1], rot[2], rot[3]);
    if (c.cmd == 1) {
        printf("din=%016" PRIx64 " pairs=%" PRIu64 " log2pairs=%.3f\n", c.din, used, log2((double)used));
        for (int o = 0; o < c.nout; o++)
            printf("dout=%016" PRIx64 " hits=%" PRIu64 " log2p=%.4f\n", c.dout[o], cnt[o],
                   cnt[o] ? log2((double)cnt[o] / (double)used) : -INFINITY);
        printf("set hits=%" PRIu64 " log2p=%.4f\n", setc, setc ? log2((double)setc / (double)used) : -INFINITY);
        if (c.mode != M_INDEP) {
            printf("perkey_hits_dout0=");
            for (int k = 0; k < c.nkeys; k++) printf("%s%" PRIu64, k ? "," : "", keycnt[k]);
            printf("\nperkey_hits_set=");
            for (int k = 0; k < c.nkeys; k++) printf("%s%" PRIu64, k ? "," : "", keyset[k]);
            printf("\n");
        }
    } else if (c.cmd == 2) {
        printf("din=%016" PRIx64 " pairs=%" PRIu64 "\n", c.din, used);
        for (int b = 0; b < 64; b++) {
            double corr = 1.0 - 2.0 * (double)bit[b] / (double)used;
            printf("bit=%d ones=%" PRIu64 " corr=%.8f log2abs=%.4f\n", b, bit[b], corr,
                   corr != 0 ? log2(fabs(corr)) : -INFINITY);
        }
    } else {
        double s2 = 0, s1 = 0;
        for (int k = 0; k < c.nkeys; k++) { s2 += keycorr[k] * keycorr[k]; s1 += keycorr[k]; }
        double N = (double)(total / BS_BLOCKS * BS_BLOCKS), elp = s2 / c.nkeys - 1.0 / N;
        printf("amask=%016" PRIx64 " bmask=%016" PRIx64 " samples_per_key=%" PRIu64 "\n", c.amask, c.bmask, total / BS_BLOCKS * BS_BLOCKS);
        printf("mean_corr=%.8f mean_sq=%.10e elp_debiased=%.10e log2sqrt_elp=%.4f\n", s1 / c.nkeys,
               s2 / c.nkeys, elp, elp > 0 ? 0.5 * log2(elp) : -INFINITY);
        printf("perkey_corr=");
        for (int k = 0; k < c.nkeys; k++) printf("%s%.6f", k ? "," : "", keycorr[k]);
        printf("\n");
    }
    return 0;
}
