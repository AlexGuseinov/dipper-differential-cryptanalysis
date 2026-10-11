#!/usr/bin/env python3
"""Step 03: lower bound for the iterative differential under independent round keys.

From step 02 (exact): 0x800000 and 0x900000 both go to 0x800000 with probability 1/8
and to 0x900000 with probability 3/32.  With independent uniform round keys the round
transitions form a Markov chain, so summing only the paths that stay in {0x800000, 0x900000}
gives
    Pr[0x800000 -> 0x800000, r rounds] >= (1/8 + 3/32)^(r-1) * 1/8 = (1/8)(7/32)^(r-1),
    Pr[0x800000 -> {0x800000, 0x900000}, r rounds] >= (7/32)^r.
"""
from fractions import Fraction
import math

a, b = Fraction(1, 8), Fraction(3, 32)
print(" r  Pr[a->a] >= (log2)   Pr[a->{a,b}] >= (log2)")
for r in range(1, 29):
    p = a * (a + b) ** (r - 1)
    q = (a + b) ** r
    print(f"{r:2d}  {math.log2(p):10.4f}           {math.log2(q):10.4f}")
print(f"per-round decay: log2(32/7) = {math.log2(Fraction(32, 7)):.4f} bits")
