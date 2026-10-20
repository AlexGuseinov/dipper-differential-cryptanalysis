#!/bin/sh
# Step 13: the fix of Section 7 evaluated for r_C in {3, 9, 11} (r_A=1, r_B=4, r_D=11 kept),
# with the original r_C=7 for comparison.
#  13_lowgraph_k1_rC* : exact graph of differences with <= 1 active nibble (K=1)
#                       best iterative trail and cluster decay
#  13_sat_diff_rC*    : optimal differential trail weight, standard model, r = 1..8
#  13_sat_lin_rC*     : optimal linear trail weight, standard model, r = 1..8
#  13_lin_search_rC*  : exact iterative linear approximations among masks with <= 2 active nibbles
cd "$(dirname "$0")/.."
T=${THREADS:-4}
for rc in 7 3 9 11; do
  mkdir -p results/state/lowgraph_k1_rC$rc
  steps/run.sh 13_lowgraph_k1_rC$rc ./bin/lowgraph -k 1 -m -16 -R 1,4,$rc,11 -t $T -c results/state/lowgraph_k1_rC$rc
done
for rc in 3 9 11; do
  steps/run.sh 13_sat_diff_rC$rc python3 steps/11_sat.py best -m diff -r 8 -R 1,4,$rc,11 --trail -t 3600 -w $T
  steps/run.sh 13_sat_lin_rC$rc python3 steps/11_sat.py best -m lin -r 8 -R 1,4,$rc,11 --trail -t 3600 -w $T
done
for rc in 3 9 11; do
  steps/run.sh 13_lin_search_rC$rc ./bin/lintool search -k 2 -m -5 -L 4 -R 1,4,$rc,11 &
done
wait
# exact graph of differences with <= 2 active nibbles for r_C = 3 and 9 (resumable: rerun to continue)
for rc in 3 9; do
  mkdir -p results/state/lowgraph_k2_rC$rc
  steps/run.sh 13_lowgraph_k2_rC$rc ./bin/lowgraph -k 2 -m -16 -R 1,4,$rc,11 -t $T -c results/state/lowgraph_k2_rC$rc
done
