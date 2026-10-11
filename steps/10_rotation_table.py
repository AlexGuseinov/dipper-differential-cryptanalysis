#!/usr/bin/env python3
"""Step 10: single-bit self-loops and best single-bit cycles for all rotation choices.

The single-bit graph depends only on the rotations of the added words (r_A, r_C);
r_B and r_D do not enter it.  For every (r_A, r_C) in 0..15 x 0..15 this prints the
self-loops, the best simple cycle of length <= 6, and its decay in bits per round.
Writes results/data/rotation_table.csv.
"""
import csv, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
g = importlib.import_module("07_single_bit_graph")

T = g.load()
rows = []
for ra in range(16):
    for rc in range(16):
        E = g.graph(T, (ra, 0, rc, 0))
        loops = sorted((u, E[(u, v)]) for (u, v) in E if u == v)
        cyc = g.best_cycles(E, maxlen=6)
        best = cyc[0] if cyc else None
        rows.append({
            "rA": ra, "rC": rc,
            "self_loops": " ".join(f"{u}:{p}" for u, p in loops) or "-",
            "max_self_loop_log2": f"{max(math.log2(p) for _, p in loops):.3f}" if loops else "-",
            "best_cycle": " ".join(map(str, best[1])) if best else "-",
            "best_cycle_prob": str(best[2]) if best else "-",
            "best_cycle_bits_per_round": f"{-best[0]:.3f}" if best else "-",
        })
out = os.path.join(g.root, "results", "data", "rotation_table.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
noloop = [r for r in rows if r["self_loops"] == "-"]
print(f"(rA,rC) pairs without a single-bit self-loop: {len(noloop)} of 256")
for r in rows:
    if (r["rA"], r["rC"]) == (1, 7): print("current (1,7):", r)
print("\nself-loop map (row rA, column rC; '.' = no self-loop, else -log2 of the best self-loop):")
print("rA\\rC " + " ".join(f"{c:>2}" for c in range(16)))
for ra in range(16):
    line = []
    for rc in range(16):
        r = rows[ra * 16 + rc]
        line.append(" ." if r["self_loops"] == "-" else f"{-float(r['max_self_loop_log2']):>2.0f}")
    print(f"{ra:>5} " + " ".join(line))
acyc = [r for r in noloop if r["best_cycle"] == "-"]
print(f"\npairs without self-loop and without any single-bit cycle of length <= 6: {len(acyc)}")
print("  " + " ".join(f"({r['rA']},{r['rC']})" for r in acyc))
cyc = sorted((r for r in noloop if r["best_cycle"] != "-"), key=lambda r: float(r["best_cycle_bits_per_round"]))
print(f"pairs without self-loop but with a longer single-bit cycle: {len(cyc)}; strongest:")
for r in cyc[:10]:
    print(f"  rA={r['rA']:2d} rC={r['rC']:2d} cycle {r['best_cycle']} prob {r['best_cycle_prob']} = {r['best_cycle_bits_per_round']} bits/round")
