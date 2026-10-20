#!/usr/bin/env python3
"""Step 11: optimal differential and linear trails in the standard (independence) model, by CP-SAT.

Model of one round (any rotation constants):
  differential: each S-box is a DDT entry, weight -log2(DDT/16); each 16-bit addition follows
                Lipmaa-Moriai (validity and weight HW(~eq(a,b,c)) on bits 0..14); rotations and
                the bit permutation are wiring.
  linear:       each S-box is a LAT entry, weight -log2|cor|; each addition follows Wallen
                (z_15 = 0, z_i = z_{i+1} ^ u_{i+1} ^ v_{i+1} ^ w_{i+1}, u^v <= z, u^w <= z,
                weight HW(z)); the passed word B (or D) branches into the adder.
The model treats S-boxes and additions as independent.  For Dipper that assumption fails
(Section 3: the S-box output enters the addition with no key in between), so the weights are
those of the standard model, not exact probabilities.  Every exactly possible differential trail
is valid in the model (the model's support contains the exact support); exact probabilities of
the trails found are computed separately with `onestep round`.

S-box weights are stored in thousandths of a bit (log2(16/6) = 1.415).

  11_sat.py selftest                       model of the addition checked against brute force (n=5)
  11_sat.py best -m diff|lin -r R [-R a,b,c,d] [-t SECONDS] [--trail]
       optimal r-round trail weight for r = 1..R (each r solved separately, proven optimal
       unless the time limit is hit, which is reported)
"""
import argparse, itertools, math, sys, time
from ortools.sat.python import cp_model

SBOX = [0x1, 0xA, 0x4, 0xC, 0x6, 0xF, 0x3, 0x9, 0x2, 0xD, 0xB, 0x7, 0x5, 0x0, 0x8, 0xE]
SCALE = 1000

def perm_bit(i):
    return 4 * (i // 16) + 16 * ((3 * ((i % 16) // 4) + (i % 4)) % 4) + (i % 4)

def ddt():
    T = [[0] * 16 for _ in range(16)]
    for x in range(16):
        for d in range(16):
            T[d][SBOX[x] ^ SBOX[x ^ d]] += 1
    return T

def lat():
    T = [[0] * 16 for _ in range(16)]
    for a in range(16):
        for b in range(16):
            T[a][b] = sum(1 if bin((a & x) ^ (b & SBOX[x])).count("1") % 2 == 0 else -1 for x in range(16))
    return T  # correlation = T / 16

def wt(p):  # weight in thousandths of a bit for probability or |correlation| p
    return round(-math.log2(p) * SCALE)

def bits(v, n):
    return [(v >> i) & 1 for i in range(n)]

# ---- addition constraints (shared by the model and the self-test) ----
def lm_valid(a, b, c, n):
    """Lipmaa-Moriai: (a,b) -> c for n-bit addition possible; weight."""
    eq = lambda i: (a >> i & 1) == (b >> i & 1) == (c >> i & 1)
    if ((a ^ b ^ c) & 1) != 0: return None
    for i in range(1, n):
        if eq(i - 1) and (((a ^ b ^ c) >> i) ^ (b >> (i - 1))) & 1: return None
    return sum(0 if eq(i) else 1 for i in range(n - 1))

def wallen(u, v, w, n):
    """Wallen: output mask u, input masks v, w of n-bit addition; weight HW(z) or None."""
    z, zi = 0, 0
    for i in range(n - 2, -1, -1):
        zi ^= ((u ^ v ^ w) >> (i + 1)) & 1
        z |= zi << i
    if (u ^ v) & ~z & ((1 << n) - 1) or (u ^ w) & ~z & ((1 << n) - 1): return None
    return bin(z).count("1")

def selftest(n=5):
    M = 1 << n
    bad = 0
    for a in range(M):
        for b in range(M):
            cnt = [0] * M
            for x in range(M):
                for y in range(M):
                    cnt[((x + y) % M) ^ (((x ^ a) + (y ^ b)) % M)] += 1
            for c in range(M):
                w = lm_valid(a, b, c, n)
                if (cnt[c] == 0) != (w is None) or (w is not None and cnt[c] != M * M >> w):
                    bad += 1
    print(f"Lipmaa-Moriai, n={n}: {bad} mismatches out of {M**3}")
    bad2 = 0
    for u in range(M):
        for v in range(M):
            for w in range(M):
                s = sum((-1) ** (bin((u & ((x + y) % M)) ^ (v & x) ^ (w & y)).count("1") & 1)
                        for x in range(M) for y in range(M))
                c = abs(s) / (M * M)
                ww = wallen(u, v, w, n)
                if (c == 0) != (ww is None) or (ww is not None and c != 2.0 ** -ww):
                    bad2 += 1
    print(f"Wallen, n={n}: {bad2} mismatches out of {M**3}")
    return bad + bad2

# ---- CP-SAT model ----
def build(R, rots, mode):
    m = cp_model.CpModel()
    X = [[m.NewBoolVar(f"x{r}_{i}") for i in range(64)] for r in range(R + 1)]
    weights = []
    D, L = ddt(), lat()
    if mode == "diff":
        stab = [(*bits(d, 4), *bits(o, 4), wt(D[d][o] / 16)) for d in range(16) for o in range(16) if D[d][o]]
    else:
        stab = [(*bits(a, 4), *bits(b, 4), wt(abs(L[a][b]) / 16)) for a in range(16) for b in range(16) if L[a][b]]
    swts = sorted({t[-1] for t in stab})
    for r in range(R):
        x = X[r]
        s = [m.NewBoolVar(f"s{r}_{i}") for i in range(64)]
        for j in range(16):
            w = m.NewIntVarFromDomain(cp_model.Domain.FromValues(swts), f"ws{r}_{j}")
            m.AddAllowedAssignments(x[4 * j:4 * j + 4] + s[4 * j:4 * j + 4] + [w], stab)
            weights.append(w)
        # rotation: t[base + (k + rot) % 16] = s[base + k]
        t = [None] * 64
        for wd in range(4):
            base = 48 - 16 * wd
            for k in range(16):
                t[base + (k + rots[wd]) % 16] = s[base + k]
        u = [None] * 64  # state after ARX, before the permutation
        for (a, b) in ((48, 32), (16, 0)):
            if mode == "diff":
                al, be = t[a:a + 16], t[b:b + 16]
                ga = [m.NewBoolVar(f"g{r}_{a}_{k}") for k in range(16)]
                m.AddBoolXOr([al[0], be[0], ga[0], m.NewConstant(1)])  # a0^b0^c0 = 0
                e = []
                for k in range(16):
                    ek = m.NewBoolVar(f"e{r}_{a}_{k}")  # ek <=> al==be==ga
                    m.AddAllowedAssignments([al[k], be[k], ga[k], ek],
                                            [(p, q, z, int(p == q == z)) for p in (0, 1) for q in (0, 1) for z in (0, 1)])
                    e.append(ek)
                for k in range(1, 16):
                    # e_{k-1} => al_k ^ be_k ^ ga_k ^ be_{k-1} = 0
                    m.AddBoolXOr([al[k], be[k], ga[k], be[k - 1], m.NewConstant(1)]).OnlyEnforceIf(e[k - 1])
                for k in range(15):
                    weights.append(SCALE * (1 - e[k]))
                for k in range(16):
                    u[a + k] = ga[k]; u[b + k] = be[k]
            else:
                # masks: output mask on A2 (ua), on B out (ub); adder input masks: v on t_A, w on t_B
                ua = [m.NewBoolVar(f"ua{r}_{a}_{k}") for k in range(16)]
                ub = [m.NewBoolVar(f"ub{r}_{a}_{k}") for k in range(16)]
                v, tb = t[a:a + 16], t[b:b + 16]  # mask on t_A is v; mask on t_B is w ^ ub
                w = [m.NewBoolVar(f"w{r}_{a}_{k}") for k in range(16)]
                for k in range(16):
                    m.AddBoolXOr([tb[k], w[k], ub[k], m.NewConstant(1)])  # tb = w ^ ub
                z = [m.NewBoolVar(f"z{r}_{a}_{k}") for k in range(16)]
                m.Add(z[15] == 0)
                for k in range(15):
                    # z_k = z_{k+1} ^ ua_{k+1} ^ v_{k+1} ^ w_{k+1}
                    m.AddBoolXOr([z[k], z[k + 1], ua[k + 1], v[k + 1], w[k + 1], m.NewConstant(1)])
                for k in range(16):
                    # (ua ^ v) <= z and (ua ^ w) <= z
                    m.AddAllowedAssignments([ua[k], v[k], w[k], z[k]],
                                            [(p, q, s_, zz) for p in (0, 1) for q in (0, 1) for s_ in (0, 1) for zz in (0, 1)
                                             if zz or (p == q and p == s_)])
                for k in range(15):
                    weights.append(SCALE * z[k])
                for k in range(16):
                    u[a + k] = ua[k]; u[b + k] = ub[k]
        for i in range(64):
            m.Add(X[r + 1][perm_bit(i)] == u[i])
    m.AddBoolOr(X[0])
    obj = sum(weights)
    m.Minimize(obj)
    return m, X

def hexv(sol, xs):
    return sum(sol.Value(b) << i for i, b in enumerate(xs))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd")
    ap.add_argument("-m", default="diff")
    ap.add_argument("-r", type=int, default=4)
    ap.add_argument("-R", default="1,4,7,11")
    ap.add_argument("-t", type=float, default=3600)
    ap.add_argument("-w", type=int, default=4)
    ap.add_argument("--from", dest="r0", type=int, default=1)
    ap.add_argument("--trail", action="store_true")
    a = ap.parse_args()
    if a.cmd == "selftest":
        sys.exit(1 if selftest() else 0)
    rots = [int(v) for v in a.R.split(",")]
    print(f"mode={a.m} rotations={rots} standard model, weights in bits (S-box weights to 0.001)")
    for r in range(a.r0, a.r + 1):
        m, X = build(r, rots, a.m)
        sv = cp_model.CpSolver()
        sv.parameters.max_time_in_seconds = a.t
        sv.parameters.num_workers = a.w
        t0 = time.time()
        st = sv.Solve(m)
        dt = time.time() - t0
        name = sv.StatusName(st)
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            print(f"r={r} status={name} weight={sv.ObjectiveValue() / SCALE:.3f} "
                  f"lower_bound={sv.BestObjectiveBound() / SCALE:.3f} seconds={dt:.1f}", flush=True)
            if a.trail:
                print("  trail: " + " -> ".join(f"{hexv(sv, X[k]):016x}" for k in range(r + 1)), flush=True)
        else:
            print(f"r={r} status={name} seconds={dt:.1f}", flush=True)

if __name__ == "__main__":
    main()
