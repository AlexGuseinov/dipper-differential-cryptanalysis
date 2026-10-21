#!/bin/sh
# Step 15: Dipper-32 (src/toy32.h), a 32-bit variant with Dipper's iterative differential
# (self-loop 0x8000 -> 0x8000 with probability 2^-3), to check the full-round data
# complexity of step 06 against the full codebook.
#   model : toysearch (choice of the variant), toygraph (Markov model), 15_toy_model.py
#   indep : Markov process measured with independent round keys
#   book  : exact counts over all 2^32 plaintexts for 64 keys of the Dipper-32/48 key schedule
#   table : 15_toy_table.py (per-key probabilities, success of the distinguisher)
# usage: steps/15_toy.sh model|indep|book|table
set -e
cd "$(dirname "$0")/.."
T=${THREADS:-$(getconf _NPROCESSORS_ONLN)}
mkdir -p results/data/toy
case "$1" in
model)
    steps/run.sh 15_toysearch ./bin/toysearch 20
    steps/run.sh 15_toygraph ./bin/toygraph -k 3 -m -20 -d 16 -o results/data/toy/model
    steps/run.sh 15_toy_model python3 steps/15_toy_model.py 8 16 10 13 ;;
indep)
    steps/run.sh 15_toy_indep ./bin/toymeasure indep -d 16 -l 36 -s 1 -t "$T" ;;
book)
    steps/run.sh 15_toy_book ./bin/toymeasure -s 1 -n 64 -d 16 -S results/data/toy/targets.txt -t "$T" \
        -o results/data/toy/hits.bin ;;
table)
    steps/run.sh 15_toy_table python3 steps/15_toy_table.py ;;
*) echo "usage: $0 model|indep|book|table"; exit 1 ;;
esac
