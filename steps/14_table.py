#!/usr/bin/env python3
"""Step 14 (table): summary of the key-recovery runs results/raw/14_krec_r*_s*.txt.

Per run: rounds R, key schedule, seed, chosen plaintexts, survivors of the 59-bit filter,
right pairs (Delta x_{R-2} = 0x800000, checked with the known key), pairs with
Delta x_{R-2} = 0x900000, rank of the correct 32-bit guess (number of guesses with a
higher score), number of guesses kept (score >= correct), bits of the 32 determined,
seconds for data and ranking. Then per R: totals and the range over keys.
"""
import glob, math, os, re

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def parse(path):
    t = open(path).read()
    g = lambda pat: re.search(pat, t)
    r = {}
    r["R"] = int(g(r"rounds=(\d+)").group(1))
    r["ks"] = g(r"key=(ks\d+)").group(1)
    r["seed"] = int(g(r"seed=(\d+)").group(1))
    r["cp"] = float(g(r"chosen plaintexts=2\^([\d.]+)").group(1))
    r["dsec"] = float(g(r"data: .* seconds=([\d.]+)").group(1))
    m = g(r"survivors=(\d+) right\(Delta x_\{R-2\}=0x800000\)=(\d+) nibble5=(\d+) other=(\d+)")
    r["surv"], r["right"], r["nib5"], r["other"] = map(int, m.groups())
    m = g(r"\(0x100000\*v\):(.*)")
    byv = dict(kv.split(":") for kv in m.group(1).split())
    r["v9"] = int(byv["9"])
    m = g(r"guesses with a higher score=(\d+) with the same score=(\d+)")
    if m:
        r["above"], r["tie"] = int(m.group(1)), int(m.group(2))
        r["bits"] = 32 - math.log2(r["above"] + r["tie"])
        r["rsec"] = float(g(r"ranking: seconds=([\d.]+)").group(1))
        r["check"] = "ok" in g(r"ranking selfcheck.*").group(0)
    else:
        r["above"] = r["tie"] = None; r["bits"] = 0.0; r["rsec"] = 0.0; r["check"] = None
    return r


def main():
    paths = glob.glob(os.path.join(root, "results", "raw", "14_krec_r*_s*.txt"))
    runs = [parse(p) for p in paths if "# seconds:" in open(p).read()]  # finished runs only
    runs.sort(key=lambda r: (r["R"], r["seed"]))
    print("Key recovery, differential over rounds 1..R-2, 32 bits of P^-1(rk_R) guessed, round R peeled")
    print(f"{'R':>3} {'key':>6} {'seed':>4} {'log2 CP':>7} {'surv':>5} {'right':>5} {'0x9..':>5} {'other':>5} "
          f"{'higher':>6} {'kept':>9} {'bits':>6} {'data s':>8} {'rank s':>6}")
    for r in runs:
        kept = (r["above"] + r["tie"]) if r["tie"] is not None else None
        print(f"{r['R']:>3} {r['ks']:>6} {r['seed']:>4} {r['cp']:7.3f} {r['surv']:>5} {r['right']:>5} {r['v9']:>5} "
              f"{r['other']:>5} {r['above'] if r['above'] is not None else '-':>6} "
              f"{('2^%.3f' % math.log2(kept)) if kept else '-':>9} {r['bits']:6.2f} {r['dsec']:8.1f} {r['rsec']:6.1f}")
    print("\nper R: keys, correct guess ranked first (no higher score) in all runs, bits determined min/median/max, "
          "right pairs min/median/max")
    for R in sorted({r["R"] for r in runs}):
        rs = [r for r in runs if r["R"] == R]
        bits = sorted(r["bits"] for r in rs); right = sorted(r["right"] for r in rs)
        first = sum(1 for r in rs if r["above"] == 0)
        print(f"  R={R}: {len(rs)} keys, 2^{rs[0]['cp']:.0f} CP, first in {first}/{len(rs)}, "
              f"bits {bits[0]:.2f}/{bits[len(bits) // 2]:.2f}/{bits[-1]:.2f}, "
              f"right pairs {right[0]}/{right[len(right) // 2]}/{right[-1]}, survivors not in nibble 5: "
              f"{sum(r['other'] for r in rs)}")


if __name__ == "__main__":
    main()
