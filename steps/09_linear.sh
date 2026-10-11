#!/bin/sh
# Step 09: linear side.
#  09_lin_search_k2 : best iterative approximations among masks with <= 2 active nibbles
#  09_edge_*        : exact one-round correlations of the best cycle, by the formula (corr)
#                     and by exhaustive enumeration of both halves (enum)
#  09_lin_r*        : measured correlations of the cycle trail on 1..5 rounds, 16 keys
#                     with fixed independent round keys (rk) and the 64/96 key schedule
cd "$(dirname "$0")/.."
T=${THREADS:-4}
M0=0000000000000013; M1=00000000000b0000; M2=0000001900000000; M3=000000000c000600
[ -f results/raw/09_lin_search_k2.txt ] || steps/run.sh 09_lin_search_k2 ./bin/lintool search -k 2 -m -5 -L 4
i=0
for pair in $M0:$M1 $M1:$M2 $M2:$M3 $M3:$M0; do
  a=${pair%%:*}; b=${pair##*:}
  steps/run.sh 09_edge_${i}_corr ./bin/lintool corr -a $a -b $b
  steps/run.sh 09_edge_${i}_enum ./bin/lintool enum -a $a -b $b
  i=$((i+1))
done
set -- $M1 $M2 $M3 $M0 $M1
for spec in 1:20 2:22 3:26 4:30 5:32; do
  r=${spec%%:*}; n=${spec##*:}
  eval b=\${$r}
  for mode in rk ks96; do
    steps/run.sh 09_lin_r${r}_$mode ./bin/dltool lin -r $r -a $M0 -b $b -n $n -K 16 -k $mode -s 9$r -t $T
  done
done
