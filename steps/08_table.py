#!/usr/bin/env python3
"""Summarise step 08: largest single-bit differential-linear correlation per round
(input difference 0x800000, independent round keys), the per-round decay, and the data
needed by the DL test (1/c^2 pairs) against the differential (1/p pairs, p from step 03 bound,
which step 04 measured up to 13 rounds).  Note: the maximum over 64 output bits is biased
upward when the correlation is close to the noise level 1/sqrt(pairs)."""
import math, os, re
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def bound(r): return math.log2(1/8) + (r - 1) * math.log2(7/32)
rows = []
print("| r | log2 pairs | best bit | corr | log2 abs | noise log2 (1/sqrt N) | DL data log2(1/c^2) | diff data log2(1/p) | fewer pairs |")
print("|---|---|---|---|---|---|---|---|---|")
for r in range(1, 12):
    t = open(os.path.join(root, "results", "raw", f"08_dl_r{r}.txt")).read()
    N = int(re.search(r"pairs=(\d+)", t).group(1))
    best = max((abs(float(c)), int(b), float(c)) for b, c in re.findall(r"bit=(\d+) ones=\d+ corr=(-?[\d.]+)", t))
    la = math.log2(best[0]); rows.append((r, la))
    dl, df = -2 * la, -bound(r)
    print(f"| {r} | {math.log2(N):.0f} | {best[1]} | {best[2]:+.8f} | {la:.3f} | {-math.log2(N)/2:.2f} | {dl:.2f} | {df:.2f} | {'DL' if dl < df else 'differential'} |")
for a, b in ((7, 11), (8, 11), (6, 11)):
    la, lb = dict(rows)[a], dict(rows)[b]
    print(f"\ndecay of the best DL correlation, rounds {a}..{b}: {(la - lb) / (b - a):.3f} bits/round", end="")
print()
