#!/usr/bin/env python3
"""Print the 16-element target set of the nibble setting (step 06) as a dltool -o list."""
import importlib.util, os, sys
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("dc", os.path.join(here, "06_data_complexity.py"))
ARGS = sys.argv[1:]
dc = importlib.util.module_from_spec(spec); sys.argv = ["x", "--none"]; spec.loader.exec_module(dc)
nodes, edges = dc.load_cluster(os.path.join(dc.root, "results", "data", "cluster_800000.txt"))
n = len(nodes); idx = {d: i for i, d in nodes.items()}
final = {}
for x in range(1, 16):
    w = [0.0] * n
    for o, p in dc.one_round(x << 20).items():
        if o in idx: w[idx[o]] += p
    tot = {}
    w10 = list(w)
    for r in range(dc.ROUNDS - 1): w = dc.propagate(w, edges, n)
    for r in range(9): w10 = dc.propagate(w10, edges, n)
    final[x] = (w, w10)
tot = [sum(final[x][0][i] for x in final) for i in range(n)]
order = sorted(range(n), key=lambda i: -tot[i])[:16]
if ARGS and ARGS[0] == "--predict":
    import math
    for x in range(1, 16):
        p = sum(final[x][1][i] for i in order)
        print(f"model {x << 20:#09x} r=10 -> 16-set: " + (f"2^{math.log2(p):.3f}" if p else "0"))
else:
    print(",".join(f"{nodes[i]:x}" for i in order))
