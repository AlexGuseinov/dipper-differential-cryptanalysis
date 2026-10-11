#!/usr/bin/env python3
"""Step 02: exact one-round transitions of 0x800000 and 0x900000.

Runs `onestep round` (exhaustive over each 32-bit half) and compares the exact
probability with the product model (S-box DDT x Lipmaa-Moriai xdp+), which treats
the S-box output and the addition as independent.
"""
import os, subprocess
from fractions import Fraction

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S = [1, 10, 4, 12, 6, 15, 3, 9, 2, 13, 11, 7, 5, 0, 8, 14]
ROT = {"C": 7, "D": 11}
M16 = 0xFFFF

def P(i): return 4*(i//16) + 16*((3*((i % 16)//4) + (i % 4)) % 4) + (i % 4)
def perm(s): return sum(((s >> i) & 1) << P(i) for i in range(64))
def rotl(x, r): return ((x << r) | (x >> (16 - r))) & M16 if r else x
def hw(x): return bin(x).count("1")

def xdp_add(a, b, g, n=16):
    """Lipmaa-Moriai: Pr_{x,y}[((x^a)+(y^b)) ^ (x+y) = g] for n-bit words."""
    m = (1 << n) - 1
    eq = lambda x, y, z: (~x ^ y) & (~x ^ z) & m
    if eq(a << 1, b << 1, g << 1) & (a ^ b ^ g ^ (b << 1)) & m:
        return Fraction(0)
    return Fraction(1, 2 ** hw(~eq(a, b, g) & (m >> 1)))

def exact(din):
    out = subprocess.run([os.path.join(root, "bin", "onestep"), "round", "-i", f"{din:x}", "-m", "-6"],
                         capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        if line.startswith("dout="):
            f = dict(kv.split("=") for kv in line.split())
            res[int(f["dout"], 16)] = Fraction(int(f["numL"]) * int(f["numR"]), 2 ** 64)
    return res

def product_model(din, dout):
    """Sum over S-box output differences of DDT * xdp+, for differences active only in word C."""
    nib = (din >> 16) & M16
    assert din >> 32 == 0 and din & M16 == 0
    tot = Fraction(0)
    act = [j for j in range(4) if (nib >> (4*j)) & 15]
    assert len(act) == 1
    j = act[0]; a = (nib >> (4*j)) & 15
    for b in range(16):
        ddt = sum(1 for x in range(16) if S[x] ^ S[x ^ a] == b)
        if not ddt: continue
        alpha = rotl(b << (4*j), ROT["C"])
        for g in range(1, 1 << 16):
            if perm(g << 16) == dout:
                tot += Fraction(ddt, 16) * xdp_add(alpha, 0, g)
    return tot

for din in (0x800000, 0x900000):
    ex = exact(din)
    for dout in (0x800000, 0x900000):
        e = ex.get(dout, Fraction(0))
        pm = product_model(din, dout)
        print(f"{din:#x} -> {dout:#x}: exact {e} (2^{float(__import__('math').log2(e)):.4f}), "
              f"product model {pm}")
    best = sorted(ex.items(), key=lambda kv: -kv[1])[:4]
    print("   top outputs:", ", ".join(f"{d:#018x}:{p}" for d, p in best))

# Per S-box output difference: product model vs exact (enumeration over the 16-bit
# word C and a carry DP over the uniform operand D1).
def add_count(u, u2, g):
    st = {(0, 0): 1}
    for i in range(16):
        ns = {}
        ui, u2i, gi = (u >> i) & 1, (u2 >> i) & 1, (g >> i) & 1
        for (c, c2), w in st.items():
            for y in (0, 1):
                if (ui ^ y ^ c) ^ (u2i ^ y ^ c2) != gi: continue
                k = ((ui & y) | (ui & c) | (y & c), (u2i & y) | (u2i & c2) | (y & c2))
                ns[k] = ns.get(k, 0) + w
        st = ns
    return sum(st.values())

def sub16(x): return sum(S[(x >> (4*j)) & 15] << (4*j) for j in range(4))

print("\nper-characteristic breakdown (S-box output difference b in nibble 1 of C):")
for din in (0x800000, 0x900000):
    nib = (din >> 16) & M16; a = (nib >> 4) & 15
    for dout in (0x800000, 0x900000):
        gs = [g for g in range(1, 1 << 16) if perm(g << 16) == dout]
        assert len(gs) == 1; g = gs[0]
        for b in range(16):
            ddt = sum(1 for x in range(16) if S[x] ^ S[x ^ a] == b)
            if not ddt: continue
            pm = Fraction(ddt, 16) * xdp_add(rotl(b << 4, 7), 0, g)
            if pm == 0: continue
            num = 0
            for C in range(1 << 16):
                if ((sub16(C) ^ sub16(C ^ nib)) >> 4) & 15 != b: continue
                num += add_count(rotl(sub16(C), 7), rotl(sub16(C ^ nib), 7), g)
            print(f"  {din:#x}->{dout:#x} b={b:#x}: product {pm}, exact {Fraction(num, 2**32)}")
