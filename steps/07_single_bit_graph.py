#!/usr/bin/env python3
"""Step 07: exact single-bit transition graph of one Dipper round.

Nodes are the 64 one-bit differences.  An edge u -> v carries the exact probability
that a one-bit input difference u gives the one-bit output difference v after one round
(uniform state / uniform round key).  Only bits of the added words A and C have outgoing
edges: a difference in B or D passes through the S-box into the passed word and also
changes the sum, so the output has at least two active bits.  Edges from word A depend
only on the rotation of A, edges from word C only on the rotation of C.

  07_single_bit_graph.py [--rot a,b,c,d]    print the graph and its best cycles
  07_single_bit_graph.py --check            check the claims used in the paper
Exact counts come from `bin/onestep bitgraph` (cached in results/data/bitgraph.txt).
"""
import argparse, math, os, subprocess, sys
from fractions import Fraction

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cache = os.path.join(root, "results", "data", "bitgraph.txt")

def P(i): return 4*(i//16) + 16*((3*((i % 16)//4) + (i % 4)) % 4) + (i % 4)

def load():
    if not os.path.exists(cache):
        out = subprocess.run([os.path.join(root, "bin", "onestep"), "bitgraph"], capture_output=True,
                             text=True, check=True).stdout
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        open(cache, "w").write(out)
    T = {}
    for line in open(cache):
        f = dict(kv.split("=") for kv in line.split())
        T[(int(f["r"]), int(f["i"]), int(f["j"]))] = Fraction(int(f["num"]), 2 ** int(f["log2den"]))
    return T

def graph(T, rot):
    ra, rc = rot[0] % 16, rot[2] % 16
    E = {}
    for (r, i, j), p in T.items():
        if r == ra: E[(48 + i, P(48 + j))] = p
        if r == rc: E[(16 + i, P(16 + j))] = p
    return E

def best_cycles(E, maxlen=4):
    """All simple cycles up to maxlen; returns list of (mean log2 p per round, cycle, prob)."""
    adj = {}
    for (u, v), p in E.items(): adj.setdefault(u, []).append((v, p))
    res = []
    def dfs(start, u, path, prob):
        for v, p in adj.get(u, []):
            if v == start:
                res.append((math.log2(prob * p) / len(path), tuple(path), prob * p))
            elif v > start and v not in path and len(path) < maxlen:
                dfs(start, v, path + [v], prob * p)
    for s in sorted(adj): dfs(s, s, [s], Fraction(1))
    res.sort(key=lambda x: -x[0])
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rot", default="1,4,7,11")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    T = load()
    if a.check:
        ok = True
        E = graph(T, (1, 4, 7, 11))
        p = E.get((23, 23), 0)
        print(f"(1,4,7,11): self-loop at bit 23 with probability {p}", "OK" if p == Fraction(1, 8) else "FAIL")
        ok &= p == Fraction(1, 8)
        loops = sorted((u, E[(u, v)]) for (u, v) in E if u == v)
        print(f"(1,4,7,11): all self-loops: {[(u, str(q)) for u, q in loops]}")
        E0 = graph(T, (0, 0, 0, 0))
        loops0 = [u for (u, v) in E0 if u == v]
        print(f"(0,0,0,0): self-loops: {loops0}", "OK (none)" if not loops0 else "FAIL")
        ok &= not loops0
        c0 = best_cycles(E0)
        print("(0,0,0,0): best cycle:", f"{c0[0][1]} prob {c0[0][2]} ({c0[0][0]:.3f} bits/round)" if c0 else "none")
        sys.exit(0 if ok else 1)
    rot = tuple(int(x) for x in a.rot.split(","))
    E = graph(T, rot)
    print(f"rotation {rot}: {len(E)} edges")
    for (u, v), p in sorted(E.items()):
        print(f"  {u:2d} -> {v:2d}  {p}  (2^{math.log2(p):.3f})")
    print("best cycles (up to length 4):")
    for m, c, p in best_cycles(E)[:10]:
        print(f"  {c}  prob {p}  {m:.3f} log2/round")

if __name__ == "__main__":
    main()
