#!/bin/sh
# usage: steps/run.sh NAME COMMAND...
# Runs COMMAND from the repository root, saves its output to results/raw/NAME.txt
# and records the command, host and date in results/COMMANDS.md.
set -e
cd "$(dirname "$0")/.."
name=$1; shift
out=results/raw/$name.txt
start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
{ echo "# command: $*"; echo "# host: $(uname -srm), $(getconf _NPROCESSORS_ONLN 2>/dev/null || echo ?) cpus"; echo "# start: $start"; } > "$out"
t0=$(date +%s)
"$@" >> "$out" 2>&1
t1=$(date +%s)
echo "# seconds: $((t1 - t0))" >> "$out"
printf -- '- `%s` -> `%s` (%s, %ss)\n' "$*" "$out" "$start" "$((t1 - t0))" >> results/COMMANDS.md
