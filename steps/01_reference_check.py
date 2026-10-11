#!/usr/bin/env python3
"""Compare the C engine with the published Python reference implementation.

Uses dipper_reference.py from the Dipper-lightweight-cipher repository
(default: ../Dipper-lightweight-cipher/dipper_reference.py, or $DIPPER_REF).
Checks 2000 random (key, plaintext) pairs for each key size, full 28 rounds.
"""
import importlib.util, os, random, subprocess, sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
ref = os.environ.get("DIPPER_REF", os.path.join(root, "..", "Dipper-lightweight-cipher", "dipper_reference.py"))
if not os.path.exists(ref):
    print(f"reference not found at {ref}; set DIPPER_REF"); sys.exit(1)
spec = importlib.util.spec_from_file_location("dipper_reference", ref)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

rnd = random.Random(2026)
lines, exp = [], []
for bits in (96, 128):
    for _ in range(2000):
        k, p = rnd.getrandbits(bits), rnd.getrandbits(64)
        lines.append(f"{bits} {k:0{bits//4}x} {p:016x}")
        exp.append(mod.encrypt(p, k, key_size=bits))
out = subprocess.run([os.path.join(root, "bin", "dltool"), "enc"], input="\n".join(lines) + "\n",
                     capture_output=True, text=True, check=True).stdout.split()
bad = sum(int(o, 16) != e for o, e in zip(out, exp))
print(f"C engine vs Python reference: {len(exp)} encryptions, {bad} mismatches {'OK' if bad == 0 else 'FAIL'}")
sys.exit(1 if bad else 0)
