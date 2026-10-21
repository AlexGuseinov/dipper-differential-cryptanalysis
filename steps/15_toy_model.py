#!/usr/bin/env python3
"""Step 15 (model): data complexity of the Dipper-32 distinguisher, computed exactly as in
step 06 for Dipper but with the toy's codebook 2^32.

Input: results/data/toy/model_d<d>.txt from `toygraph` (top outputs after d rounds:
"dout sig p8000", sig = sum over the 15 nibble-3 input differences, p8000 from 0x8000;
lower bounds under independent round keys).
A distinguisher on R = d + 1 rounds peels the keyless last round and looks at the
d-round difference. For each d and setting it reports the smallest N (chosen plaintexts,
log2 in steps of 1/16) with success >= 0.90 at false alarm <= 0.01 (Poisson counts) and
writes the 16-element target sets to results/data/toy/targets.txt for `toymeasure`.

usage: 15_toy_model.py [DMIN DMAX [TMIN TMAX]]   (target sets are written for d = TMIN..TMAX)
"""
import math, os, sys

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR = os.path.join(root, "results", "data", "toy")
ALPHA, SUCCESS, NBITS = 0.01, 0.90, 32


def pois_sf(k, lam):
    if k <= 0: return 1.0
    s, term = 0.0, math.exp(-lam)
    for j in range(k):
        s += term; term *= lam / (j + 1)
    return max(0.0, 1.0 - s)


def best_test(pairs, pS, size):
    lam0 = pairs * size * 2.0 ** -NBITS
    lam1 = lam0 + pairs * pS
    tau = 0
    while pois_sf(tau, lam0) > ALPHA: tau += 1
    return tau, lam0, lam1, pois_sf(tau, lam1)


def min_log2N(pairs_per_pt, pS, size):
    for t in range(16 * 8, 16 * NBITS + 1):
        lN = t / 16
        tau, l0, l1, s = best_test(pairs_per_pt * 2.0 ** lN, pS, size)
        if s >= SUCCESS: return lN, tau, l0, l1, s
    return None


def load(d):
    rows = []
    for line in open(os.path.join(DIR, f"model_d{d}.txt")):
        o, sig, p8 = line.split()
        rows.append((int(o, 16), float(sig), float(p8)))
    return rows


def main():
    dmin, dmax = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (8, 16)
    tmin, tmax = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (dmin, dmax)
    print(f"Dipper-32 model (independent round keys, <= 3 active nibbles); success >= {SUCCESS}, "
          f"false alarm <= {ALPHA}; N = chosen plaintexts, codebook 2^{NBITS}")
    print(f"{'d':>3} {'R':>3} {'setting':8} {'|S|':>4} {'log2 p_S':>9} {'log2 N':>7} {'tau':>4} {'lam0':>8} {'lam1':>8} {'succ':>6}")
    targets = []
    for d in range(dmin, dmax + 1):
        rows = load(d)
        by_sig = sorted(rows, key=lambda r: -r[1])
        by_p8 = sorted(rows, key=lambda r: -r[2])
        for k in (1, 16):
            pS = sum(r[2] for r in by_p8[:k])
            res = [("single", k, pS, min_log2N(0.5, pS, k))]
            sig = sum(r[1] for r in by_sig[:k]) * 8 / 16  # per plaintext
            res.append(("nibble", k, sig / (120 / 16), min_log2N(120 / 16, sig / (120 / 16), k)))
            for name, kk, p, r in res:
                if r is None:
                    print(f"{d:>3} {d + 1:>3} {name:8} {kk:>4} {math.log2(p):9.3f}   >{NBITS}")
                else:
                    lN, tau, l0, l1, s = r
                    print(f"{d:>3} {d + 1:>3} {name:8} {kk:>4} {math.log2(p):9.3f} {lN:7.3f} {tau:>4} {l0:8.3f} {l1:8.3f} {s:6.3f}")
        if tmin <= d <= tmax: targets.append((d, [r[0] for r in by_sig[:16]]))
    with open(os.path.join(DIR, "targets.txt"), "w") as f:
        f.write("# d, then the 16 outputs with the largest structure signal (toygraph model)\n")
        for d, outs in targets:
            f.write(f"{d} " + " ".join(f"{o:08x}" for o in outs) + "\n")
    print(f"wrote {os.path.join(DIR, 'targets.txt')}")


if __name__ == "__main__":
    main()
