#!/usr/bin/env python3
"""Summarise step 05: fixed-key measurements of 0x800000 -> 0x800000.
For each run: pooled probability with 95% interval, and the spread of the per-key counts.
Index of dispersion D = var/mean of the per-key counts; D ~ 1 if every key has the same
probability (Poisson sampling only). The chi-square statistic (K-1)D is compared with its
K-1 degrees of freedom by a normal approximation (z)."""
import math, os, re
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def bound(r): return math.log2(1/8) + (r - 1) * math.log2(7/32)
print("| mode | r | keys | pairs/key | hits | log2 p (95% interval) | bound | D | z | per-key log2 p: min / 10% / median / 90% / max |")
print("|---|---|---|---|---|---|---|---|---|---|")
for mode in ("ks96", "ks128", "rk"):
    for r in (6, 8, 10):
        f = os.path.join(root, "results", "raw", f"05_{mode}_r{r}.txt")
        if not os.path.exists(f): continue
        t = open(f).read()
        if "perkey_hits_dout0" not in t: continue
        N = int(re.search(r"pairs=(\d+)", t).group(1))
        K = int(re.search(r"nkeys=(\d+)", t).group(1))
        ks = [int(x) for x in re.search(r"perkey_hits_dout0=([\d,]+)", t).group(1).split(",")]
        h = sum(ks); per = N // K
        lp = math.log2(h / N); d = 1.96 * math.sqrt(h)
        m = h / K; var = sum((x - m) ** 2 for x in ks) / (K - 1); D = var / m
        z = ((K - 1) * D - (K - 1)) / math.sqrt(2 * (K - 1))
        sk = sorted(ks); q = lambda f: sk[min(K - 1, int(f * K))]
        qs = [sk[0], q(0.1), q(0.5), q(0.9), sk[-1]]
        L = lambda x: f"{math.log2(x / per):.2f}" if x else "-inf"
        print(f"| {mode} | {r} | {K} | 2^{math.log2(per):.0f} | {h} | {lp:.3f} [{math.log2((h-d)/N):.3f}, {math.log2((h+d)/N):.3f}] | {bound(r):.3f} | {D:.2f} | {z:+.1f} | {' / '.join(L(x) for x in qs)} |")
