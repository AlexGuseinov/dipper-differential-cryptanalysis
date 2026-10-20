#!/bin/sh
# Step 11b: exact one-round correlations (lintool corr) along the optimal 10-round linear trail of 11_sat_lin
cd "$(dirname "$0")/.."
set -- 0000a01000000000 000000000a000000 0000000000000011 0000000600000000 0000000004000000 0000000000220000 000000bb00000000 0000000000000600 0000000000000002 0000000000090000 0002000000000044
while [ $# -ge 2 ]; do ./bin/lintool corr -a $1 -b $2; shift; done | awk '{print} {split($4,a,"="); s+=a[2]} END {printf "sum of log2|cor| over 10 rounds = %.4f\n", s}'
