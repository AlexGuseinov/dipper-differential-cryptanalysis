#!/usr/bin/env python3
"""Step 17 (model): expected survivors of the key-recovery attack (step 14).

Model of step 06 (independent round keys; first round exact, then the 64-element cluster):
p(d) = average over the 15 nibble-5 input differences of Pr[difference after d rounds is
0x800000 or 0x900000], with d = R-2. A structure gives 120 pairs, so a run with 2^s
structures has mean signal lam1 = 120 * 2^s * p(d). A wrong pair passes the filter with
probability 31 * 2^-64 (equal z except in the 5 top bits of word C, and z0 != z1), so
the mean noise is lam0 = 120 * 2^s * 31 * 2^-64.

Prints lam1 against the survivors of all finished 14_krec runs, then lam1 and lam0 for
the full cipher. Inferred beyond 13 rounds, as every number of step 06.
"""
import glob, importlib.util, math, os, re, sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("dc", os.path.join(here, "06_data_complexity.py"))
dc = importlib.util.module_from_spec(spec)
argv, sys.argv = sys.argv, ["x", "--none"]
spec.loader.exec_module(dc)
sys.argv = argv
root = dc.root
TARGETS = (0x800000, 0x900000)


def p_of_d(dmax):
    """p(d) for d = 1..dmax (d = 1: exact first round restricted to the cluster)."""
    nodes, edges = dc.load_cluster(os.path.join(root, "results", "data", "cluster_800000.txt"))
    n = len(nodes)
    idx = {d: i for i, d in nodes.items()}
    t = [idx[x] for x in TARGETS]
    tot = [0.0] * (dmax + 1)
    for x in range(1, 16):
        w = [0.0] * n
        for o, p in dc.one_round(x << 20).items():
            if o in idx: w[idx[o]] += p
        for d in range(1, dmax + 1):
            tot[d] += sum(w[i] for i in t)
            if d < dmax: w = dc.propagate(w, edges, n)
    return [v / 15 for v in tot]


def observed():
    runs = {}
    for path in glob.glob(os.path.join(root, "results", "raw", "14_krec_r*_s*.txt")):
        t = open(path).read()
        if "# seconds:" not in t: continue
        R = int(re.search(r"rounds=(\d+)", t).group(1))
        s = int(re.search(r"log2structures=(\d+)", t).group(1))
        n = int(re.search(r"survivors=(\d+)", t).group(1))
        runs.setdefault((R, s), []).append(n)
    return runs


def main():
    p = p_of_d(27)
    print("key-recovery survivors: model mean against the 14_krec runs")
    print(f"{'R':>3} {'log2 p(R-2)':>11} {'log2 CP':>7} {'keys':>4} {'mean surv':>9} {'model':>8} {'ratio':>6}")
    tot_o = tot_m = 0.0
    for (R, s), v in sorted(observed().items()):
        m = 120 * 2.0 ** s * p[R - 2]
        o = sum(v) / len(v)
        tot_o += sum(v); tot_m += m * len(v)
        print(f"{R:>3} {math.log2(p[R - 2]):11.3f} {s + 4:>7} {len(v):>4} {o:9.2f} {m:8.2f} {o / m:6.3f}")
    print(f"all runs: {tot_o:.0f} survivors, model {tot_m:.1f}, ratio {tot_o / tot_m:.3f}")
    print("\nfull cipher, R = 28: d = 26 rounds, 16-plaintext structures")
    print(f"log2 p(26) = {math.log2(p[26]):.3f}")
    print(f"{'log2 CP':>7} {'lam1':>8} {'lam0':>8}")
    for c in (56, 57, 58, 59, 60, 61):
        pairs = 120 * 2.0 ** (c - 4)
        print(f"{c:>7} {pairs * p[26]:8.3f} {pairs * 31 * 2.0 ** -64:8.3f}")


if __name__ == "__main__":
    main()
