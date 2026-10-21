#!/usr/bin/env python3
"""Step 15 (table): Dipper-32 model against the full codebook.

Inputs: results/raw/15_toy_book.txt and results/data/toy/hits.bin (toymeasure, 64 keys of
the Dipper-32/48 key schedule, all 2^32 plaintexts), results/raw/15_toy_indep.txt
(independent round keys), results/data/toy/model_d*.txt (toygraph), and the test of
steps/15_toy_model.py.

Table A: Pr[0x8000 -> 0x8000] after d rounds: model lower bound (exact, <= 3 nibbles),
independent round keys (measured), and the exact per-key values over the codebook
(pooled mean, quantiles over the 64 keys, index of dispersion D of the per-key counts).

Table B: the distinguisher on R = d + 1 rounds with structures on nibble 3 and the 16
targets of the model. N chosen plaintexts are taken as N/16 structures, each structure
of the codebook kept independently with probability q = (N/16)/2^28. For one key the
number of right pairs is then sum_i c_i Bernoulli(q) over the structures i with c_i > 0
(exact distribution from the codebook counts). The test accepts when the count is >= tau,
tau from the Poisson false-alarm rule of step 06 (alpha = 0.01). Reported: model N and
predicted success; the success at that N for each key (mean, minimum, share of keys
>= 0.9); and the smallest N at which the mean success over the keys is >= 0.9.
"""
import importlib.util, math, os, struct, sys
from collections import defaultdict

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("m", os.path.join(root, "steps", "15_toy_model.py"))
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
RAW = os.path.join(root, "results", "raw")
DATA = os.path.join(root, "results", "data", "toy")


def binom_pmf_trunc(m, q, cap):
    """P[Bin(m,q) = j] for j = 0..cap-1 and the tail mass at index cap."""
    out = [0.0] * (cap + 1)
    if m == 0: out[0] = 1.0; return out
    lp = math.log(q) if q > 0 else -math.inf
    lq = math.log1p(-q)
    acc = 0.0
    for j in range(min(cap, m + 1)):
        lc = math.lgamma(m + 1) - math.lgamma(j + 1) - math.lgamma(m - j + 1)
        out[j] = math.exp(lc + j * lp + (m - j) * lq)
        acc += out[j]
    out[cap] = max(0.0, 1.0 - acc)
    return out


def tail(groups, q, tau):
    """P[sum_c c * Bin(m_c, q) >= tau] for groups {c: m_c}."""
    if tau <= 0: return 1.0
    dist = [1.0] + [0.0] * tau  # index tau = ">= tau"
    for c, m in groups.items():
        pm = binom_pmf_trunc(m, q, tau)  # pm[j] for j < tau, pm[tau] = P[>= tau]
        new = [0.0] * (tau + 1)
        for s, ps in enumerate(dist):
            if ps == 0: continue
            if s == tau: new[tau] += ps; continue
            for j, pj in enumerate(pm):
                if pj == 0: continue
                t = s + c * j if j < tau else tau
                new[min(t, tau)] += ps * pj
        dist = new
    return dist[tau]


def parse_book():
    c8, sets = {}, defaultdict(dict)
    for line in open(os.path.join(RAW, "15_toy_book.txt")):
        f = line.split()
        if len(f) > 2 and f[0] == "key" and f[2] == "c8:":
            c8[int(f[1])] = [int(x) for x in f[3:]]
        elif len(f) > 4 and f[0] == "key" and f[2] == "set":
            d = int(f[4].split("=")[1]); sets[int(f[1])][d] = int(f[5].split("=")[1])
    return c8, sets


def parse_hits():
    h = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))  # key -> d -> count -> m
    raw = open(os.path.join(DATA, "hits.bin"), "rb").read()
    for off in range(0, len(raw), 10):
        k, d, s, a, v8 = struct.unpack_from("<BBIHH", raw, off)
        h[k][d][a] += 1
    return h


def main():
    c8, sets = parse_book()
    keys = sorted(c8)
    nk = len(keys)
    indep = {}
    for line in open(os.path.join(RAW, "15_toy_indep.txt")):
        if line.startswith("r="):
            f = dict(kv.split("=") for kv in line.split())
            indep[int(f["r"])] = (int(f["hits"]), float(f["log2p"]))
    print(f"Table A: Pr[0x8000 -> 0x8000], d rounds; codebook: {nk} keys x 2^31 pairs each (exact per key)")
    print(f"{'d':>3} {'model':>8} {'indep':>8} {'(hits)':>9} {'pooled':>8} {'min':>8} {'10%':>8} {'median':>8} {'90%':>8} {'max':>8} {'D':>7}")
    D = len(c8[keys[0]])
    for d in range(1, D + 1):
        rows = [l.split() for l in open(os.path.join(DATA, f"model_d{d}.txt"))]
        mod = next((float(r[2]) for r in rows if int(r[0], 16) == 0x8000), 0.0)
        cnt = sorted(c8[k][d - 1] for k in keys)
        tot = sum(cnt)
        mean = tot / nk
        var = sum((c - mean) ** 2 for c in cnt) / (nk - 1)
        disp = var / mean if mean > 0 else float("nan")
        lg = lambda c: math.log2(c / 2 ** 31) if c > 0 else float("-inf")
        q = lambda p: cnt[min(nk - 1, int(p * nk))]
        ih, il = indep.get(d, (0, float("nan")))
        print(f"{d:>3} {math.log2(mod):8.3f} {il:8.3f} {ih:>9} {lg(tot / nk):8.3f} {lg(cnt[0]):8.3f} "
              f"{lg(q(0.1)):8.3f} {lg(q(0.5)):8.3f} {lg(q(0.9)):8.3f} {lg(cnt[-1]):8.3f} {disp:7.2f}")
    hits = parse_hits()
    print(f"\nTable B: distinguisher on R = d + 1 rounds, structures on nibble 3, 16 targets; {nk} keys")
    print(f"{'d':>3} {'R':>3} {'model log2N':>11} {'tau':>4} {'pred succ':>9} {'lam1 model':>10} {'mean count':>10} "
          f"{'mean succ':>9} {'min succ':>8} {'keys>=0.9':>9} {'measured log2N':>14}")
    tdir = {}
    for line in open(os.path.join(DATA, "targets.txt")):
        if line.startswith("#"): continue
        f = line.split(); tdir[int(f[0])] = [int(x, 16) for x in f[1:]]
    for d in sorted(tdir):
        rows = M.load(d)
        sig = sum(r[1] for r in sorted(rows, key=lambda r: -r[1])[:16]) * 8 / 16
        ppp = 120 / 16
        r = M.min_log2N(ppp, sig / ppp, 16)
        if r is None:
            print(f"{d:>3} {d + 1:>3} model needs more than the codebook"); continue
        lN, tau, l0, l1, s = r

        def success(lN_):
            q = 2.0 ** lN_ / 16 / 2 ** 28
            t, _, _, _ = M.best_test(ppp * 2.0 ** lN_, sig / ppp, 16)
            return [tail(hits[k][d], q, t) for k in keys], t
        succ, _ = success(lN)
        mean_count = sum(sets[k][d] for k in keys) / nk * (2.0 ** lN / 2 ** 32)
        meas = None
        for t16 in range(int(16 * (lN - 4)), 16 * 32 + 1):
            ss, _ = success(t16 / 16)
            if sum(ss) / nk >= 0.9: meas = t16 / 16; break
        print(f"{d:>3} {d + 1:>3} {lN:11.3f} {tau:>4} {s:9.3f} {l1:10.3f} {mean_count:10.3f} {sum(succ) / nk:9.3f} "
              f"{min(succ):8.3f} {sum(x >= 0.9 for x in succ):>6}/{nk} {meas if meas is not None else float('nan'):14.3f}")


if __name__ == "__main__":
    main()
