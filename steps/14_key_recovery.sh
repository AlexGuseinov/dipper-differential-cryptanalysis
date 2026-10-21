#!/bin/sh
# Step 14: key recovery on R-round Dipper with the real key schedules (one random master
# key per seed; seeds 1-4: Dipper-64/96, seeds 5-8: Dipper-64/128).
# Differential 0x800000 -> {0x800000, 0x900000} over rounds 1..R-2 from structures of 16
# plaintexts on nibble 5, 32 last-round key bits guessed, round R peeled without a key.
# usage: steps/14_key_recovery.sh R LOG2STRUCT [SEEDS]
set -e
cd "$(dirname "$0")/.."
R=$1; N=$2; SEEDS=${3:-"1 2 3 4 5 6 7 8"}
T=${THREADS:-$(getconf _NPROCESSORS_ONLN)}
for s in $SEEDS; do
    if [ "$s" -le 4 ]; then k=ks96; else k=ks128; fi
    steps/run.sh 14_krec_r${R}_s${s} ./bin/krec -r "$R" -k $k -s "$s" -n "$N" -t "$T" -T 8,9
done
