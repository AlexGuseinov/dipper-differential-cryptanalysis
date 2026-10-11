# Results

Labels: **exact** = computed by exhaustive enumeration or exact arithmetic;
**measured** = sampled on the cipher (pairs, hits, 95% interval from the Poisson count);
**inferred** = beyond the measured rounds, under the independent-round-key (Markov) model.
Raw outputs: `results/raw/`; commands: `results/COMMANDS.md`.
Differences are 64-bit hex, bit 0 = least significant; 0x800000 = bit 23 (word C, bit 7).

## 0. Implementation checks (`00_make_test`)
Spec test vectors TV1–TV4 reproduce; the table-based round equals the bit-level reference round
on 10^6 random inputs; the bitsliced engine equals the table engine on 256 blocks x 28 rounds;
4000 full encryptions agree with `dipper_reference.py` from the cipher repository.

## 1. Exact one-round transitions (`02_one_round_exact`, exact)
Enumeration of all 2^32 inputs of the (C,D) half (the (A,B) half has zero difference):

| transition | exact | product model (DDT x xdp+, summed over S-box outputs) |
|---|---|---|
| 0x800000 -> 0x800000 | 1/8 | 7/64 |
| 0x800000 -> 0x900000 | 3/32 | 7/64 |
| 0x900000 -> 0x800000 | 1/8 | 7/64 |
| 0x900000 -> 0x900000 | 3/32 | 11/128 |

For 0x800000 -> 0x800000 only the S-box output difference 0x3 contributes. Its product-model
probability is 1/4 x 1/4 = 1/16, the exact probability is 1/8: the S-box output enters the
modular addition with no key in between, so the S-box value and the carry are correlated.
The S-box output differences 0x7 and 0xF have product-model probabilities 1/32 and 1/64 but are
impossible (exact 0). Other one-round outputs of 0x800000 with probability 1/8: 0x2000800000.

## 2. Lower bound (`03_bound`, exact arithmetic, valid under independent round keys)
Pr[0x800000 -> 0x800000, r rounds] >= (1/8)(7/32)^(r-1); decay log2(32/7) = 2.1926 bits/round.
r = 10: 2^-22.734; r = 13: 2^-29.312; r = 26: 2^-57.816; r = 27: 2^-60.009 (inferred).
Pr[0x800000 -> {0x800000, 0x900000}] >= (7/32)^r; r = 27: 2^-59.201 (inferred).

## 3. Measurements, independent round keys per pair (`04_iter_r*`, measured)

| r | log2 pairs | hits | measured log2 p | 95% interval | bound | measured - bound |
|---|---|---|---|---|---|---|
| 1 | 28 | 33554577 | -3.000 | [-3.000, -3.000] | -3.000 | +0.000 |
| 2 | 28 | 7341453 | -5.192 | [-5.193, -5.191] | -5.193 | +0.000 |
| 3 | 28 | 1605225 | -7.386 | [-7.388, -7.383] | -7.385 | -0.000 |
| 4 | 30 | 1403715 | -9.579 | [-9.582, -9.577] | -9.578 | -0.001 |
| 5 | 30 | 306818 | -11.773 | [-11.778, -11.768] | -11.771 | -0.002 |
| 6 | 30 | 66926 | -13.970 | [-13.981, -13.959] | -13.963 | -0.006 |
| 7 | 30 | 14990 | -16.128 | [-16.152, -16.105] | -16.156 | +0.028 |
| 8 | 30 | 3136 | -18.385 | [-18.437, -18.336] | -18.349 | -0.037 |
| 9 | 32 | 2826 | -20.535 | [-20.590, -20.483] | -20.541 | +0.006 |
| 10 | 32 | 657 | -22.640 | [-22.755, -22.534] | -22.734 | +0.094 |
| 11 | 34 | 527 | -24.958 | [-25.087, -24.840] | -24.926 | -0.032 |
| 12 | 36 | 454 | -27.173 | [-27.313, -27.046] | -27.119 | -0.054 |
| 13 | 37 | 179 | -29.516 | [-29.745, -29.319] | -29.312 | -0.204 |

The bound is attained within sampling error for r <= 12: paths leaving {0x800000, 0x900000}
add nothing measurable to 0x800000 -> 0x800000. At r = 13 the measurement is 0.20 bits below the
bound and its 95% interval just misses it (179 hits against 206 expected). With fresh keys in every
round the measured process is exactly the Markov chain, so the expectation cannot be below the bound;
this is read as a sampling deviation (about 1.9 sigma). A repeat with another seed is queued (`04_iter_r13b`).

The published paper (Table 20) reports 2^-11.98, 2^-14.76, 2^-18.35 at 5, 6, 7 rounds for input
0x00E00000 and extrapolates about 3.2 bits per round, reaching 2^-64 "around r = 21".
The decay measured here is 2.19 bits per round over 13 rounds.

## 4. Full cipher (`06_cluster`, `06_data`, `06_validate_model`, `06b_*`; inferred)
Round 28 is followed by no key, so the attacker inverts it on the ciphertexts without the key; a
27-round differential therefore gives a 28-round distinguisher.

Cluster: 64 differences reached from 0x800000 by exact transitions of probability >= 2^-7
(`onestep cluster`). It does not raise Pr[0x800000 -> 0x800000] above the two-state bound, but adds
other output differences. Model predictions agree with measurements at reduced rounds:
4-set {0x800000, 0x900000, 0x2000800000, 0x2000900000}: model 2^-16.721 / 2^-25.492 at r = 8 / 12,
measured 2^-16.708 / 2^-25.523. Nibble-5 input differences to the 16-element target set at r = 10
(model vs measured, 2^32 pairs each): all 15 agree within sampling error, e.g. 0x600000 2^-19.990
vs 2^-19.954, 0x300000 2^-22.798 vs 2^-22.886; the three differences with model 0 (0x1, 0x7, 0xC)
gave 0 hits.

Data complexity, success >= 0.90 at false-alarm <= 0.01 (Poisson counts; N chosen plaintexts):

| setting | target set size | log2 p per pair (27 r) | log2 N |
|---|---|---|---|
| pairs (x, x^0x800000) | 1 | -60.009 | 63.375 |
| pairs (x, x^0x800000) | 16 | -57.265 | 61.062 |
| 16-plaintext structures on nibble 5 | 1 | -60.865 | 60.562 |
| 16-plaintext structures on nibble 5 | 16 | -58.122 | 58.312 |

All four are below the codebook (2^64). Best: 2^58.3 chosen plaintexts (structures, 16 targets,
threshold 7 right pairs, expected 2.3 for a random permutation and 10.9 for Dipper).
These numbers are inferred (Markov model beyond 13 rounds) and are upper bounds on the data within
that model, because the cluster probabilities are lower bounds.

## 5. Single-bit graph (`07_single_bit_graph --check`, exact)
For (1,4,7,11) the only single-bit self-loop is bit 23 -> bit 23, probability 1/8.
With all rotations 0 there is no self-loop (best single-bit cycle: (20,22,52,30), probability 2^-16,
4 bits per round).

## 6. Rotation constants (`10_rotation_table`, exact)
The single-bit graph depends only on the rotations of the added words, r_A and r_C.
A self-loop exists exactly when r_A or r_C is in {2, 5, 6, 7, 12}; its probability is 1/8 for 7 and
1/16 otherwise. 121 of the 256 pairs (r_A, r_C) have no self-loop; 90 of them have no single-bit cycle
of length <= 6. Full table: `results/data/rotation_table.csv`.
Limitation: this is the single-bit graph only; multi-bit iterative differentials are not excluded by it.

## 7. Linear side (`09_*`)
Exact one-round correlations by formula (S-box LAT and exact carry-matrix correlation of the
addition), cross-checked by exhaustive enumeration (`09_edge_*`).
Best iterative approximation among masks with <= 2 active nibbles (exact, search threshold 2^-5 per
round): the 4-round cycle 0x13 -> 0xb0000 -> 0x1900000000 -> 0xc000600 -> 0x13 with correlations
-2^-3, 2^-2, -2^-4, 2^-2: |c| = 2^-11 per 4 rounds, 2.75 bits per round.
Ten rounds of this cycle give 2^-27.5, weaker than the published optimal 10-round trail (weight 21,
Table 21), so it does not contradict it. Reduced-round measurements: pending.

## 8. Fixed keys (`05_*`) and differential-linear (`08_*`)
Pending.
