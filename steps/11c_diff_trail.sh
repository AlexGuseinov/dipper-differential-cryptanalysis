#!/bin/sh
# Step 11c: exact one-round probabilities (onestep round) along the optimal 7-round differential
# trail of 11_sat_diff (standard-model weight 27)
cd "$(dirname "$0")/.."
set -- 0000000000e00000 0000000000800000 0000000000800000 0000000000900000 0000000000800000 0000000000900000 0040000000000000 0000800000002000
while [ $# -ge 2 ]; do
  p=$(./bin/onestep round -i $1 -m -20 | awk -v o=$2 '$1 == "dout=" o {split($5,a,"="); print a[2]}')
  echo "$1 -> $2 log2p=${p:-none}"
  shift
done | awk '{print} {split($4,a,"="); s+=a[2]} END {printf "sum of log2 p over 7 rounds = %.4f\n", s}'
