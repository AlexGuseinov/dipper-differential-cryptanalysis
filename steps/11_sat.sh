#!/bin/sh
# Step 11: optimal trails in the standard (independence) model, CP-SAT (needs: pip install ortools).
#  11_sat_selftest : the addition constraints against brute force on 5-bit words
#  11_sat_diff     : optimal differential trail weight, r = 1..10
#  11_sat_lin      : optimal linear trail weight, r = 1..10
cd "$(dirname "$0")/.."
W=${THREADS:-4}
steps/run.sh 11_sat_selftest python3 steps/11_sat.py selftest
steps/run.sh 11_sat_diff python3 steps/11_sat.py best -m diff -r 10 --trail -t 3600 -w $W
steps/run.sh 11_sat_lin python3 steps/11_sat.py best -m lin -r 10 --trail -t 3600 -w $W
