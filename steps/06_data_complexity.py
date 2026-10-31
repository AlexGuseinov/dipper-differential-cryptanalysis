#!/usr/bin/env python3
"""Step 06: data complexity of the full-round distinguisher.

Model (independent round keys, Markov): the probability that a pair with input
difference d has difference delta after 27 rounds is bounded below by the cluster
matrix of `onestep cluster` (results/data/cluster_800000.txt).  Round 28 has no key
after it, so the attacker peels it off the ciphertexts (keyless inverse) and observes
the 27-round output difference directly.

Settings compared (N = number of chosen plaintexts, N <= 2^64):
  single : pairs (x, x ^ 0x800000), N/2 pairs, target set S
  nibble : structures of 16 plaintexts that take all values in nibble 5 (bits 20..23);
           120 pairs per structure, 8 pairs for each nonzero nibble difference
For a target set S the count of pairs whose 27-round output difference is in S is
Poisson with mean lam1 = pairs * (p_S + |S| 2^-64) for Dipper and lam0 = pairs * |S| 2^-64
for a random permutation.  The test accepts "Dipper" when count >= tau.  For each N
we pick tau as the smallest value with false-alarm <= ALPHA and report the success rate.
All 27-round probabilities are lower bounds under the independent-round-key model,
i.e. "inferred" numbers. Adding the background to them assumes that all other trails
behave like it; step 06c gives the values without that assumption, which are the upper
bounds within the model.
"""
import math, os, subprocess, sys

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROUNDS = 27
ALPHA = 0.01
SUCCESS = 0.90

def load_cluster(path):
    nodes, edges = {}, []
    for line in open(path):
        if line.startswith("n "):
            _, i, d = line.split(); nodes[int(i)] = int(d, 16)
        elif line.startswith("e "):
            _, a, b, p = line.split(); edges.append((int(a), int(b), float(p)))
    return nodes, edges

def propagate(v, edges, n):
    w = [0.0] * n
    for a, b, p in edges:
        if v[a]: w[b] += v[a] * p
    return w

def one_round(din):
    out = subprocess.run([os.path.join(root, "bin", "onestep"), "round", "-i", f"{din:x}", "-m", "-12"],
                         capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        if line.startswith("dout="):
            f = dict(kv.split("=") for kv in line.split())
            res[int(f["dout"], 16)] = int(f["numL"]) * int(f["numR"]) / 2.0 ** 64
    return res

def pois_sf(k, lam):
    """P[X >= k] for X ~ Poisson(lam), computed in log space."""
    if k <= 0: return 1.0
    # sum_{j<k} e^-lam lam^j / j!
    s, term = 0.0, math.exp(-lam)
    for j in range(k):
        s += term; term *= lam / (j + 1)
    return max(0.0, 1.0 - s)

def best_test(pairs, pS, size):
    lam0 = pairs * size * 2.0 ** -64
    lam1 = lam0 + pairs * pS
    tau = 0
    while pois_sf(tau, lam0) > ALPHA: tau += 1
    return tau, lam0, lam1, pois_sf(tau, lam1)

def min_log2N(pairs_per_pt, pS, size):
    """smallest N (in 1/16 steps of log2) with success >= SUCCESS at false alarm <= ALPHA."""
    for t in range(16 * 40, 16 * 64 + 1):
        lN = t / 16
        tau, l0, l1, s = best_test(pairs_per_pt * 2.0 ** lN, pS, size)
        if s >= SUCCESS: return lN, tau, l0, l1, s
    return None

def main():
    nodes, edges = load_cluster(os.path.join(root, "results", "data", "cluster_800000.txt"))
    n = len(nodes)
    idx = {d: i for i, d in nodes.items()}
    a = idx[0x800000]
    # forward vectors from a: v_r = Pr[a -> node, r rounds] lower bounds
    v = [0.0] * n; v[a] = 1.0
    V = {0: v}
    for r in range(1, ROUNDS + 1):
        v = propagate(v, edges, n); V[r] = v
    print(f"cluster: {n} differences, {len(edges)} transitions")
    for r in (10, 13, 20, 26, 27):
        print(f"  r={r}: Pr[a->a] >= 2^{math.log2(V[r][a]):.3f}, Pr[a->cluster] >= 2^{math.log2(sum(V[r])):.3f}")
    # input side: one-round distributions of the nibble-5 differences, entering the cluster
    print("\ninput differences in nibble 5 (exact first round, then 26 cluster rounds):")
    starts = {}
    for x in range(1, 16):
        d = x << 20
        dist = one_round(d)
        w = [0.0] * n
        for o, p in dist.items():
            if o in idx: w[idx[o]] += p
        starts[d] = w
        print(f"  {d:#09x}: mass into cluster after 1 round 2^{math.log2(sum(w)) if sum(w) else float('-inf'):.3f}")
    # 27-round vectors for each start: round 1 exact, rounds 2..27 cluster
    final = {}
    for d, w in starts.items():
        for r in range(ROUNDS - 1): w = propagate(w, edges, n)
        final[d] = w
    # choose target sets: top-k outputs by probability for the summed nibble setting
    tot = [sum(final[d][i] for d in final) for i in range(n)]
    order = sorted(range(n), key=lambda i: -tot[i])
    order_a = sorted(range(n), key=lambda i: -final[0x800000][i])
    print(f"\nsuccess >= {SUCCESS}, false alarm <= {ALPHA}; N in chosen plaintexts (log2), codebook 2^64")
    print(f"{'setting':8} {'|S|':>4} {'log2 p_S':>9} {'log2 N':>7} {'tau':>4} {'lam0':>8} {'lam1':>8} {'succ':>6}")
    rows = []
    for k in (1, 2, 4, 8, 16, 32, n):
        if k > n: continue
        S = order_a[:k]
        pS = sum(final[0x800000][i] for i in S)
        r = min_log2N(0.5, pS, k)
        rows.append(("single", k, pS, r))
        S = order[:k]
        # per plaintext: 120/16 pairs per structure-plaintext, signal 8/16 * sum_d p_d(S)
        sig = sum(sum(final[d][i] for i in S) for d in final) * 8 / 16  # per plaintext
        pairs_pp = 120 / 16
        r2 = min_log2N(pairs_pp, sig / pairs_pp, k)
        rows.append(("nibble", k, sig / pairs_pp, r2))
    for name, k, pS, r in rows:
        if r is None:
            print(f"{name:8} {k:>4} {math.log2(pS):9.3f}   >64")
        else:
            lN, tau, l0, l1, s = r
            print(f"{name:8} {k:>4} {math.log2(pS):9.3f} {lN:7.3f} {tau:>4} {l0:8.3f} {l1:8.3f} {s:6.3f}")
    print("\nlargest outputs for the single setting:")
    for i in order_a[:8]:
        print(f"  {nodes[i]:#018x}: 2^{math.log2(final[0x800000][i]):.3f}")

if __name__ == "__main__" and len(sys.argv) == 1:
    main()

def validation():
    """Model predictions for reduced rounds, to be compared with dltool measurements."""
    nodes, edges = load_cluster(os.path.join(root, "results", "data", "cluster_800000.txt"))
    n = len(nodes); idx = {d: i for i, d in nodes.items()}
    T4 = [0x800000, 0x900000, 0x2000800000, 0x2000900000]
    for d in (0x800000, 0x100000, 0xc00000, 0xe00000):
        dist = one_round(d)
        w = [0.0] * n
        for o, p in dist.items():
            if o in idx: w[idx[o]] += p
        for r in range(2, 13):
            w = propagate(w, edges, n)
            if r in (8, 10, 12):
                pa = w[idx[0x800000]]; p4 = sum(w[idx[t]] for t in T4)
                L = lambda x: f"{math.log2(x):.3f}" if x > 0 else "-inf"
                print(f"model {d:#09x} r={r}: Pr[->0x800000] >= 2^{L(pa)}, Pr[->4-set] >= 2^{L(p4)}")

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--validate":
    validation()
