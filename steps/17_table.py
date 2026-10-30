#!/usr/bin/env python3
"""Step 17 (table): summary of results/raw/17_check_s*.txt and 17_sim_c*_s*.txt.

Per data volume: mean signal and noise pairs, fraction of trials with the true guess
ranked first (no higher score), and the bits of the 32 determined (32 - log2(guesses with
score >= the true guess)): 10%, median, 90%. A trial with a higher wrong guess is still
counted, with its own number of kept guesses.
"""
import glob, os, re

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = re.compile(r"trial \d+ K=(\d+) W=(\d+) score=(\d+) higher=(\d+) same=(\d+) bits=([\d.-]+)")


def load(pattern):
    rows = []
    for path in glob.glob(os.path.join(root, "results", "raw", pattern)):
        rows += [tuple(float(x) for x in m) for m in T.findall(open(path).read())]
    return rows


def line(name, rows):
    n = len(rows)
    bits = sorted(r[5] for r in rows)
    first = sum(1 for r in rows if r[3] == 0)
    print(f"{name:>10} {n:>6} {sum(r[0] for r in rows) / n:7.2f} {sum(r[1] for r in rows) / n:7.2f} "
          f"{first / n:7.3f} {bits[n // 10]:7.2f} {bits[n // 2]:7.2f} {bits[9 * n // 10]:7.2f}")


def keys():
    v = []
    for path in glob.glob(os.path.join(root, "results", "raw", "17_keys_r14_s*.txt")):
        m = re.search(r"survivors=(\d+)", open(path).read())
        if m: v.append(int(m.group(1)))
    if not v: return
    n = len(v); mean = sum(v) / n
    var = sum((x - mean) ** 2 for x in v) / (n - 1)
    print(f"R=14, 2^30 CP, {n} fresh keys: mean survivors {mean:.2f} +- {(var / n) ** 0.5:.2f} (s.e.), "
          f"dispersion {var / mean:.2f}, min {min(v)}, max {max(v)}")


def main():
    keys()
    print(f"{'setting':>10} {'trials':>6} {'K':>7} {'W':>7} {'first':>7} {'b 10%':>7} {'b med':>7} {'b 90%':>7}")
    c = load("17_check_s*.txt")
    if c: line("check", c)
    for x in (56, 57, 58, 59, 60, 61):
        r = load(f"17_sim_c{x}_s*.txt")
        if r: line(f"2^{x} CP", r)


if __name__ == "__main__":
    main()
