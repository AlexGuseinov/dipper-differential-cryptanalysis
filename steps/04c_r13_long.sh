#!/bin/sh
# Step 04c: 13 rounds of the iterative differential 0x800000 -> 0x800000 with 2^40 pairs,
# independent uniform round keys per pair (same model as step 04), with bin/longdiff:
# all cores, 4096 chunks of 2^28 pairs, finished chunks kept in results/state/04_iter_r13c.state.
# If the run is interrupted, run this script again: it continues from the state file and the
# final counts are the same as for an uninterrupted run. Progress: results/state/04_iter_r13c.log
cd "$(dirname "$0")/.."
T=${THREADS:-$(getconf _NPROCESSORS_ONLN)}
mkdir -p results/state
KEEP=""
command -v caffeinate >/dev/null 2>&1 && KEEP="caffeinate -i"   # macOS: no idle sleep while it runs
exec $KEEP steps/run.sh 04_iter_r13c ./bin/longdiff -r 13 -i 800000 -o 800000,900000,2000800000,2000900000 \
  -n 40 -c 28 -s 4131 -t "$T" -f results/state/04_iter_r13c.state -l results/state/04_iter_r13c.log
