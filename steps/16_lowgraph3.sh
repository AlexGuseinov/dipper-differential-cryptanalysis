#!/bin/sh
# Step 16: exact one-round graph of all differences with <= 3 active nibbles.
# 16_hfcheck: the conditioned half distribution (src/halffast.h) equals the 2^32 brute force
#             of lowgraph on 100 random half differences with <= 3 active nibbles.
# 16_lowgraph3_k2: the new tool reproduces step 12 (K = 2) exactly.
# 16_lowgraph3_k3: K = 3 (resumable through the cache directory).
set -e
cd "$(dirname "$0")/.."
T=${THREADS:-$(getconf _NPROCESSORS_ONLN)}
steps/run.sh 16_hfcheck ./bin/hfcheck 100 2026 3
steps/run.sh 16_lowgraph3_k2 ./bin/lowgraph3 -k 2 -m -16 -t "$T"
mkdir -p results/state/lowgraph3_k3_m16
steps/run.sh 16_lowgraph3_k3 ./bin/lowgraph3 -k 3 -m -16 -t "$T" -c results/state/lowgraph3_k3_m16
