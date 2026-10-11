#!/usr/bin/env python3
"""Summarise step 04: measured log2 p with a 95% interval (normal approximation to the
Poisson count) against the bound of step 03 and the 4-set cluster model."""
import math, os, re
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def bound(r): return math.log2(1/8) + (r - 1) * math.log2(7/32)
print("| r | log2 pairs | hits | measured log2 p | 95% interval | bound | measured - bound |")
print("|---|---|---|---|---|---|---|")
for r in range(1, 14):
    t = open(os.path.join(root, "results", "raw", f"04_iter_r{r}.txt")).read()
    N = int(re.search(r"pairs=(\d+)", t).group(1))
    h = int(re.search(r"dout=0000000000800000 hits=(\d+)", t).group(1))
    lp = math.log2(h / N); d = 1.96 * math.sqrt(h)
    lo, hi = math.log2((h - d) / N), math.log2((h + d) / N)
    print(f"| {r} | {math.log2(N):.0f} | {h} | {lp:.3f} | [{lo:.3f}, {hi:.3f}] | {bound(r):.3f} | {lp - bound(r):+.3f} |")
