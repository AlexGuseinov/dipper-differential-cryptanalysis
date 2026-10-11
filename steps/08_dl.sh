#!/bin/sh
# Step 08: single-bit differential-linear correlations for input difference 0x800000,
# independent round keys per pair; all 64 single-bit output masks at once.
cd "$(dirname "$0")/.."
T=${THREADS:-4}
for spec in 1:26 2:26 3:28 4:28 5:30 6:30 7:32 8:34 9:35 10:36 11:37; do
  r=${spec%%:*}; n=${spec##*:}
  steps/run.sh 08_dl_r$r ./bin/dltool dl -r $r -i 800000 -n $n -k indep -s 8$r -t $T
done
