#!/usr/bin/env python3
"""Step 04c table: the three 13-round runs (04_iter_r13, 04_iter_r13b, 04_iter_r13c) against
the lower bound (1/8)(7/32)^12 for 0x800000 -> 0x800000. Counts are read from results/raw.
For each run and for the pooled counts: hits, expected hits under the bound, z = (hits - E)/sqrt(E)
(Poisson), measured log2 p with an exact 95% Poisson interval, and the one-sided probability
of seeing at most this many hits if the probability equals the bound."""
import math, os, re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
LOG2_BOUND = math.log2(1 / 8) + 12 * math.log2(7 / 32)
RUNS = ["04_iter_r13", "04_iter_r13b", "04_iter_r13c"]


def read(name):
    txt = open(os.path.join(ROOT, "results", "raw", name + ".txt")).read()
    pairs = int(re.search(r"pairs=(\d+)", txt).group(1))
    hits = int(re.search(r"dout=0000000000800000 hits=(\d+)", txt).group(1))
    return pairs, hits


def poisson_cdf(k, lam):
    # P[X <= k], summed in log space
    terms = [i * math.log(lam) - lam - math.lgamma(i + 1) for i in range(k + 1)]
    m = max(terms)
    return math.exp(m) * sum(math.exp(t - m) for t in terms)


def poisson_ci(k, alpha=0.05):
    # exact (Garwood) interval by bisection on the Poisson CDF
    def solve(f, lo, hi):
        for _ in range(200):
            mid = (lo + hi) / 2
            if f(mid):
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2
    lo = 0.0 if k == 0 else solve(lambda l: 1 - poisson_cdf(k - 1, l) < alpha / 2, 0, k)
    hi = solve(lambda l: poisson_cdf(k, l) > alpha / 2, k, 3 * k + 50)
    return lo, hi


rows = []
for name in RUNS:
    try:
        rows.append((name,) + read(name))
    except FileNotFoundError:
        print(f"{name}: missing")
if len(rows) > 1:
    rows.append(("pooled", sum(r[1] for r in rows), sum(r[2] for r in rows)))

print(f"bound: log2 p = {LOG2_BOUND:.4f}")
print("| run | log2 pairs | hits | expected (bound) | z | measured log2 p | 95% interval | P[X <= hits] |")
print("|---|---|---|---|---|---|---|---|")
for name, pairs, hits in rows:
    e = pairs * 2.0 ** LOG2_BOUND
    lo, hi = poisson_ci(hits)
    print(f"| {name} | {math.log2(pairs):.3f} | {hits} | {e:.1f} | {(hits - e) / math.sqrt(e):+.2f} | "
          f"{math.log2(hits / pairs):.4f} | [{math.log2(lo / pairs):.4f}, {math.log2(hi / pairs):.4f}] | "
          f"{poisson_cdf(hits, e):.3g} |")
