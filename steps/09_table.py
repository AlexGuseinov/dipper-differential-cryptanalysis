#!/usr/bin/env python3
"""Summarise step 09: one-round correlations of the 4-round cycle (formula vs enumeration),
and measured sqrt(ELP) of the cycle on 1..5 rounds against the single-trail correlation."""
import math, os, re
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
raw = lambda n: open(os.path.join(root, "results", "raw", n + ".txt")).read()
cors = []
print("| edge | input mask | output mask | formula | enumeration |")
print("|---|---|---|---|---|")
for i in range(4):
    a = re.search(r"a=(\w+) b=(\w+) cor=(\S+)", raw(f"09_edge_{i}_corr")); e = re.search(r"cor=(\S+)", raw(f"09_edge_{i}_enum"))
    cors.append(float(a.group(3)))
    print(f"| {i} | 0x{int(a.group(1),16):x} | 0x{int(a.group(2),16):x} | {float(a.group(3)):+g} | {float(e.group(1)):+g} |")
print("\n| r | samples/key x keys | trail log2 abs c | measured log2 sqrt(ELP), rk | measured, ks96 |")
print("|---|---|---|---|---|")
for r in range(1, 6):
    tr = sum(math.log2(abs(cors[i % 4])) for i in range(r))
    v = {}
    for m in ("rk", "ks96"):
        t = raw(f"09_lin_r{r}_{m}")
        v[m] = float(re.search(r"log2sqrt_elp=(\S+)", t).group(1))
        n = int(re.search(r"samples_per_key=(\d+)", t).group(1))
    print(f"| {r} | 2^{math.log2(n):.0f} x 16 | {tr:.2f} | {v['rk']:.2f} | {v['ks96']:.2f} |")
