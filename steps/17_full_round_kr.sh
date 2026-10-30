#!/bin/sh
# Step 17: the key recovery of step 14 on the full 28 rounds (inferred).
# `pool` runs krec on 14 rounds to collect real pairs with Delta x_{R-2} in {0x800000, 0x900000}
# for one key per seed. `sim` ranks all 2^32 guesses on such pairs mixed with synthetic noise,
# with the means of 17_model for 28 rounds at LOG2CP chosen plaintexts.
# `keys` repeats the 14-round attack for 256 fresh master keys (seeds 101-356) to test the mean
# number of survivors against 17_model: the runs of step 14 share 8 master keys across all R.
# usage: steps/17_full_round_kr.sh pool | keys | sim | check
set -e
cd "$(dirname "$0")/.."
T=${THREADS:-$(getconf _NPROCESSORS_ONLN)}
SEEDS="1 2 3 4 5 6 7 8"
TRIALS=${TRIALS:-40}
D=results/data/krsim
key() { grep -o 'kc=[0-9a-f]* kd=[0-9a-f]*' "results/raw/17_pool_s$1.txt" | sed 's/k[cd]=//g'; }
lam() { awk -v c="$1" -v col="$2" '$1 == c && NF == 3 { print $col }' results/raw/17_model.txt; }
case $1 in
pool)
    mkdir -p $D
    for s in $SEEDS; do
        if [ "$s" -le 4 ]; then k=ks96; else k=ks128; fi
        steps/run.sh 17_pool_s$s ./bin/krec -r 14 -k $k -s "$s" -n 30 -t "$T" -T 8,9 -o $D/pool_s$s.txt
    done ;;
keys)
    s=101
    while [ $s -le 356 ]; do
        if [ "$s" -le 228 ]; then k=ks96; else k=ks128; fi
        steps/run.sh 17_keys_r14_s$s ./bin/krec -r 14 -k $k -s "$s" -n 26 -t "$T" -T 8,9
        s=$((s + 1))
    done ;;
sim)
    for c in 56 57 58 59 60 61; do
        for s in $SEEDS; do
            set -- $(key "$s")
            steps/run.sh 17_sim_c${c}_s$s ./bin/krsim -f $D/pool_s$s.txt -c "$1" -d "$2" -l "$(lam "$c" 2)" -w "$(lam "$c" 3)" \
                -n "$TRIALS" -s "$s" -t "$T"
        done
    done ;;
check)
    # no noise, the signal of the 20-round runs: must reproduce the measured bits of step 14
    for s in $SEEDS; do
        set -- $(key "$s")
        steps/run.sh 17_check_s$s ./bin/krsim -f $D/pool_s$s.txt -c "$1" -d "$2" -l 26 -w 0 -n "$TRIALS" -s "$s" -t "$T"
    done ;;
*) echo "usage: $0 pool|keys|sim|check" >&2; exit 1 ;;
esac
