#!/bin/sh
# Step 04: measure Pr[0x800000 -> 0x800000] (and nearby outputs) for 1..13 rounds
# with independent uniform round keys drawn fresh for every pair.
# Targets: 0x800000, 0x900000, 0x2000800000, 0x2000900000.
cd "$(dirname "$0")/.."
T=${THREADS:-4}
OUT=800000,900000,2000800000,2000900000
for spec in 1:28 2:28 3:28 4:30 5:30 6:30 7:30 8:30 9:32 10:32 11:34 12:36 13:37; do
  r=${spec%%:*}; n=${spec##*:}
  steps/run.sh 04_iter_r$r ./bin/dltool diff -r $r -i 800000 -o $OUT -n $n -k indep -s 4$r -t $T
done
