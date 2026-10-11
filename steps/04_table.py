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
# repeat of r = 13 with another seed, and the two runs combined
f = os.path.join(root, "results", "raw", "04_iter_r13b.txt")
if os.path.exists(f):
    tot_h = tot_N = 0
    for name in ("04_iter_r13.txt", "04_iter_r13b.txt"):
        t = open(os.path.join(root, "results", "raw", name)).read()
        tot_N += int(re.search(r"pairs=(\d+)", t).group(1))
        tot_h += int(re.search(r"dout=0000000000800000 hits=(\d+)", t).group(1))
    t = open(f).read(); N = int(re.search(r"pairs=(\d+)", t).group(1)); h = int(re.search(r"dout=0000000000800000 hits=(\d+)", t).group(1))
    for lab, h_, N_ in (("13 (repeat)", h, N), ("13 (both)", tot_h, tot_N)):
        lp = math.log2(h_ / N_); d = 1.96 * math.sqrt(h_)
        print(f"| {lab} | {math.log2(N_):.0f} | {h_} | {lp:.3f} | [{math.log2((h_-d)/N_):.3f}, {math.log2((h_+d)/N_):.3f}] | {bound(13):.3f} | {lp - bound(13):+.3f} |")
    e = tot_N * 2 ** bound(13)
    print(f"\nboth runs: expected under the bound {e:.1f}, observed {tot_h}, z = {(tot_h - e) / math.sqrt(e):.2f}")
