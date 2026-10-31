#!/usr/bin/env python3
"""Step 01b: which 64-bit window of the 128-bit key register is the round key.

Equation (7) of the published paper takes K[127:64] for Dipper-64/128; specification v1.1 and
dipper_reference.py take K[63:0]. Encrypts the two Dipper-64/128 test vectors of the published
Table 3 with each window. Uses the reference as step 01 does ($DIPPER_REF).
"""
import importlib.util, os, sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
ref = os.environ.get("DIPPER_REF", os.path.join(root, "..", "Dipper-lightweight-cipher", "dipper_reference.py"))
if not os.path.exists(ref):
    print(f"reference not found at {ref}; set DIPPER_REF"); sys.exit(1)
spec = importlib.util.spec_from_file_location("dipper_reference", ref)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

TV = [("TV3", 0x0123456789ABCDEF, 0x000102030405060708090A0B0C0D0E0F, 0x4B46284387969060),
      ("TV4", 0xFFFFFFFFFFFFFFFF, (1 << 128) - 1, 0x6ABB4518063706B0)]
low = mod.extract_round_key
high = lambda K, *rest: (K >> 64) & 0xFFFFFFFFFFFFFFFF
for name, p, k, c in TV:
    mod.extract_round_key = low
    a = mod.encrypt(p, k, key_size=128)
    mod.extract_round_key = high
    b = mod.encrypt(p, k, key_size=128)
    print(f"{name}: published {c:016X}, K[63:0] {a:016X} {'match' if a == c else 'NO MATCH'}, "
          f"K[127:64] {b:016X} {'match' if b == c else 'NO MATCH'}")
