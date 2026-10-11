#!/bin/sh
# Step 06b: measured 10-round probabilities for all 15 input differences in nibble 5,
# to the 16-element target set used in step 06 (independent round keys).
cd "$(dirname "$0")/.."
T=${THREADS:-4}
S=$(python3 steps/06b_targets.py)
for x in 1 2 3 4 5 6 7 8 9 a b c d e f; do
  steps/run.sh 06b_nib_${x}_r10 ./bin/dltool diff -r 10 -i ${x}00000 -o $S -n 32 -k indep -s 6$x -t $T
done
