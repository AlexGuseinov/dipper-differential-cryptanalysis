# dipper-differential-cryptanalysis

Differential, differential-linear and linear analysis of the Dipper block cipher
(spec v1.1, https://github.com/AlexGuseinov/Dipper-lightweight-cipher): code and data for the
paper "Differential Cryptanalysis of Full-Round Dipper" (A. Huseynli, Y. Imamverdiyev, J. Alizadeh).

## Build and test

    make                # portable build
    make NATIVE=1       # adds -march=native (faster bitsliced engine); used for the measurements
    make test           # spec test vectors, fast vs bit-level round, bitsliced vs table engine,
                        # and 4000 encryptions against the published Python reference

On macOS, if the linker rejects the default SDK ("tapi error: malformed file"), pick an
installed one for the build only, e.g. `SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk make NATIVE=1`.

`make test` expects the cipher repository next to this one
(`../Dipper-lightweight-cipher/dipper_reference.py`) or `DIPPER_REF=<path>`.

## Tools (src/)

| binary | purpose |
|---|---|
| `bin/dltool` | measurements: `diff` (differential counts), `dl` (single-bit differential-linear correlations), `lin` (fixed-key linear correlations), `test`, `enc` |
| `bin/onestep` | exact one-round analysis by exhaustive enumeration of each 32-bit half: `round` (output distribution), `bitgraph` (single-bit graph data), `cluster` (Markov-chain lower bound over a set of differences) |
| `bin/longdiff` | long differential runs with independent round keys (as `dltool diff -k indep`): chunked, all cores, checkpoint/resume (state file), progress log; `longdiff bench -r R` measures throughput |
| `bin/lintool` | exact one-round linear correlations (formula and enumeration) and a search for iterative linear approximations |
| `bin/lowgraph` | exact one-round transition graph of all differences with at most K active nibbles: best iterative trail (minimum mean cycle) and cluster decay (spectral radius); `-c DIR` caches the half enumerations so a run can be resumed |

Key modes for `dltool`: `indep` (fresh uniform round keys for every pair: the Markov model),
`rk` (uniform round keys fixed per key), `ks96`, `ks128` (real key schedules, random master keys).

## Steps (steps/)

| step | what |
|---|---|
| 01 | C engine vs Python reference |
| 02 | exact one-round transitions of 0x800000 and 0x900000; product model |
| 03 | the lower bound (1/8)(7/32)^(r-1) |
| 04 | measurements of the iterative differential, 1..13 rounds, independent round keys |
| 04c | 13 rounds with 2^40 pairs (`longdiff`, resumable) and the table of all 13-round runs against the bound |
| 05 | fixed-key measurements: Dipper-64/96, Dipper-64/128, fixed independent round keys |
| 06 | data complexity of the full-round distinguisher |
| 07 | exact single-bit transition graph (`--check` verifies the claims) |
| 08 | single-bit differential-linear correlations |
| 09 | linear side: exact one-round correlations, iterative approximation, measurements |
| 10 | single-bit self-loops for all rotation constants |
| 11 | optimal differential and linear trails in the standard model by CP-SAT (`11_sat.py`, needs `pip install ortools`); 11b, 11c: exact correlations and probabilities along the optimal 10-round linear and 7-round differential trails |
| 12 | exact graph of all differences with at most two active nibbles (`lowgraph -k 2`) |
| 13 | the rotation fix evaluated for r_C in {3, 9, 11}: exact graphs, CP-SAT trails, iterative linear approximations |

`steps/run.sh NAME CMD...` runs a command, stores its output in `results/raw/NAME.txt`
and appends the command to `results/COMMANDS.md`. Results are summarised in
`results/RESULTS.md`; every number there is labelled exact, measured or inferred.
