#!/usr/bin/env python3
"""Step 06c: the data complexity of step 06 with and without the background term.

Step 06 uses lam1 = pairs * (p_S + |S| 2^-64). But p_S is a lower bound on the probability of
reaching S, from the trails in the cluster, so adding the random background assumes that
all other trails together behave like it. That assumption does not follow from independent
round keys. This script reports, for each row of step 06,
  background: lam1 = pairs * (p_S + |S| 2^-64), the value of step 06;
  signal only: lam1 = pairs * p_S. The mean count of the cipher is at least this, so this N
  is an upper bound within the independent-round-key model.
The threshold tau comes from the random permutation in both cases. It then tests the
assumption on Dipper-32 with the full-codebook runs of step 15 (15_toy_table).
"""
import importlib.util, math, os, re, sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("dc", os.path.join(here, "06_data_complexity.py"))
dc = importlib.util.module_from_spec(spec)
argv, sys.argv = sys.argv, ["x", "--none"]
spec.loader.exec_module(dc)
sys.argv = argv
sf = dc.pois_sf


def solve(q, pS, size, background):
    for t in range(16 * 40, 16 * 64 + 1):
        lN = t / 16
        pairs = q * 2.0 ** lN
        l0 = pairs * size * 2.0 ** -64
        tau = 0
        while sf(tau, l0) > dc.ALPHA: tau += 1
        lam = pairs * pS + (l0 if background else 0.0)
        if sf(tau, lam) >= dc.SUCCESS: return lN, tau, l0, lam, sf(tau, lam)
    return None


def row(name, q, pS, size):
    a, b = solve(q, pS, size, True), solve(q, pS, size, False)
    print(f"{name:8} {size:>3} {math.log2(pS):9.3f} | background log2N {a[0]:7.3f} tau {a[1]:>3} lam0 {a[2]:7.3f} "
          f"lam1 {a[3]:7.3f} succ {a[4]:.3f} | signal only log2N {b[0]:7.3f} tau {b[1]:>3} lam0 {b[2]:7.3f} "
          f"lam1 {b[3]:7.3f} succ {b[4]:.3f}")
    pairs = q * 2.0 ** a[0]
    print(f"         at the background N: signal-only mean {pairs * pS:.3f}, background mean {a[2]:.3f}, "
          f"success with the signal alone {sf(a[1], pairs * pS):.3f}")


def main():
    nodes, edges = dc.load_cluster(os.path.join(dc.root, "results", "data", "cluster_800000.txt"))
    n = len(nodes)
    idx = {d: i for i, d in nodes.items()}
    v = [0.0] * n; v[idx[0x800000]] = 1.0
    for r in range(dc.ROUNDS): v = dc.propagate(v, edges, n)
    final = {}
    for x in range(1, 16):
        w = [0.0] * n
        for o, p in dc.one_round(x << 20).items():
            if o in idx: w[idx[o]] += p
        for r in range(dc.ROUNDS - 1): w = dc.propagate(w, edges, n)
        final[x] = w
    tot = [sum(final[x][i] for x in final) for i in range(n)]
    order_a = sorted(range(n), key=lambda i: -v[i])
    order = sorted(range(n), key=lambda i: -tot[i])
    print(f"success >= {dc.SUCCESS}, false alarm <= {dc.ALPHA}; rows of step 06")
    for k in (1, 16):
        row("pairs", 0.5, sum(v[i] for i in order_a[:k]), k)
        row("nibble", 7.5, sum(tot[i] for i in order[:k]) * (8 / 16) / 7.5, k)
    print("\nDipper-32, step 15 (measured over each key's codebook, 64 keys):")
    print(f"{'R':>3} {'log2N':>7} {'tau':>4} {'lam0':>7} {'lam1':>7} {'mean count':>10} {'succ bg':>8} {'succ signal':>11} {'measured':>8}")
    text = open(os.path.join(dc.root, "results", "raw", "15_toy_table.txt")).read().split("Table B")[1]
    for m in re.finditer(r"^\s*(\d+)\s+(\d+)\s+([\d.]+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", text, re.M):
        d, R, lN, tau, pred, l1, mean, ms = m.groups()
        lN, tau, l1 = float(lN), int(tau), float(l1)
        l0 = 7.5 * 2.0 ** lN * 16 * 2.0 ** -32
        print(f"{R:>3} {lN:7.3f} {tau:>4} {l0:7.3f} {l1:7.3f} {float(mean):10.3f} {sf(tau, l1):8.3f} {sf(tau, l1 - l0):11.3f} {float(ms):8.3f}")


if __name__ == "__main__":
    main()
