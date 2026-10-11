#!/bin/sh
# Step 05: the iterative differential 0x800000 -> 0x800000 under fixed keys.
#   ks96 / ks128 : real key schedules, NKEYS random master keys
#   rk           : independent uniform round keys, fixed per key (same layout, no schedule)
# Per-key hit counts are printed so that the key dependence can be examined.
cd "$(dirname "$0")/.."
T=${THREADS:-4}
OUT=800000,900000,2000800000,2000900000
for mode in ks96 ks128 rk; do
  steps/run.sh 05_${mode}_r6  ./bin/dltool diff -r 6  -i 800000 -o $OUT -n 30 -k $mode -K 256 -s 56 -t $T
  steps/run.sh 05_${mode}_r8  ./bin/dltool diff -r 8  -i 800000 -o $OUT -n 33 -k $mode -K 256 -s 58 -t $T
  steps/run.sh 05_${mode}_r10 ./bin/dltool diff -r 10 -i 800000 -o $OUT -n 35 -k $mode -K 64  -s 510 -t $T
done
