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

Round key window (`01b_key_window`): Equation (7) of the published paper takes K[127:64] for Dipper-64/128,
specification v1.1 and `dipper_reference.py` take K[63:0]. The published test vectors TV3 and TV4 are reproduced
by K[63:0] only. This repository analyses the cipher defined by the test vectors. Results for independent round
keys do not depend on the window; measurements with `ks128` do.

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

Why 1/8 (`02b_carry_check`, exact over all 2^32 inputs, 0 mismatches): with u, v the inputs of the addition
C + D, the output difference of 0x800000 is 0x800000 exactly when (i) the S-box output difference is 0x3
and (ii) v_11 != c_11, the carry into bit 11; count 2^29. The S-box outputs with output difference 0x3 are
1, 2, 5, 6, whose two low bits always differ. Given (ii) the carry into bit 12 equals u_11, so the second carry
condition u_12 != c_12 holds automatically, which the product model charges a factor 1/2 for.
The S-box outputs for 0x7 (9, a, d, e) and 0xF (0, 4, b, f) contradict the carry conditions, so hits = 0.

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
| 13 (repeat, seed 4130) | 37 | 192 | -29.415 | [-29.635, -29.224] | -29.312 | -0.103 |
| 13 (both runs) | 38 | 371 | -29.465 | [-29.620, -29.325] | -29.312 | -0.153 |
| 13 (2^40 run, seed 4131, `04_iter_r13c`) | 40 | 1669 | -29.295 | [-29.365, -29.226] | -29.312 | +0.017 |
| 13 (all three runs) | 40.322 | 2040 | -29.328 | [-29.391, -29.265] | -29.312 | -0.016 |

The bound is attained within sampling error for r <= 12: paths leaving {0x800000, 0x900000}
add nothing measurable to 0x800000 -> 0x800000. At r = 13 the first two runs were below the bound; together
371 hits against 412.5 expected (z = -2.04, 0.15 bits). The 2^40-pair run (`04_iter_r13c`, step 04c)
gives 1669 hits against 1650.0 expected (z = +0.47); all three runs together give 2040 against 2062.5
(z = -0.50), and the bound lies inside the 95% interval (`04c_table`). With fresh keys in every round the measured
process is exactly the Markov chain, so its expectation cannot be below the bound, and the r = 2 value
(7.3 million hits) matches the exact two-round value to 0.0002 bits, which checks the engine.
The r = 13 gap of the first two runs was a sampling deviation: the 2^40-pair run attains the bound.

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
These numbers are inferred (Markov model beyond 13 rounds). They use lam1 = pairs * (p_S + |S| 2^-64), so they
also assume that all trails outside the cluster behave like the random background; that does not follow from
independent round keys. See section 16 for the values without this assumption, which are the upper bounds on the
data within the model. Best without it: 2^59.06 chosen plaintexts.

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
Table 21), so it does not contradict it.

Formula vs enumeration and reduced-round measurements (`09_table`; 16 keys; sqrt(ELP) debiased):

| edge | input mask | output mask | formula | enumeration |
|---|---|---|---|---|
| 0 | 0x13 | 0xb0000 | -0.125 | -0.125 |
| 1 | 0xb0000 | 0x1900000000 | +0.25 | +0.25 |
| 2 | 0x1900000000 | 0xc000600 | -0.0625 | -0.0625 |
| 3 | 0xc000600 | 0x13 | +0.25 | +0.25 |

| r | samples/key x keys | trail log2 abs c (exact) | measured log2 sqrt(ELP), fixed indep. keys | measured, Dipper-64/96 |
|---|---|---|---|---|
| 1 | 2^20 x 16 | -3.00 | -2.99 | -2.99 |
| 2 | 2^22 x 16 | -5.00 | -4.99 | -5.01 |
| 3 | 2^26 x 16 | -9.00 | -8.51 | -8.56 |
| 4 | 2^30 x 16 | -11.00 | -10.23 | -10.51 |
| 5 | 2^32 x 16 | -14.00 | -13.63 | -13.32 |

From 3 rounds the measured hull is 0.4 to 0.8 bits stronger than the single trail.

## 8. Fixed keys (`05_*`, `05_table`; measured)
Per-key counts for 256 keys (6, 8 rounds) and 64 keys (10 rounds). D = index of dispersion of the
per-key counts (1 if all keys had the same probability); z = normal approximation of the chi-square test.

| mode | r | keys | pairs/key | hits | log2 p (95% interval) | bound | D | z | per-key log2 p: min / 10% / median / 90% / max |
|---|---|---|---|---|---|---|---|---|---|
| ks96 | 6 | 256 | 2^22 | 67681 | -13.954 [-13.964, -13.943] | -13.963 | 20.85 | +224.1 | -14.80 / -14.44 / -14.04 / -13.51 / -12.90 |
| ks96 | 8 | 256 | 2^25 | 24625 | -18.412 [-18.430, -18.394] | -18.349 | 10.12 | +102.9 | -19.33 / -18.98 / -18.51 / -17.91 / -16.97 |
| ks96 | 10 | 64 | 2^29 | 4874 | -22.749 [-22.790, -22.709] | -22.734 | 12.62 | +65.2 | -24.00 / -23.51 / -22.87 / -22.03 / -21.57 |
| ks128 | 6 | 256 | 2^22 | 67604 | -13.955 [-13.966, -13.944] | -13.963 | 18.73 | +200.2 | -14.69 / -14.42 / -14.06 / -13.53 / -12.91 |
| ks128 | 8 | 256 | 2^25 | 25268 | -18.375 [-18.393, -18.357] | -18.349 | 10.32 | +105.2 | -19.39 / -18.96 / -18.48 / -17.82 / -17.03 |
| ks128 | 10 | 64 | 2^29 | 4799 | -22.771 [-22.813, -22.731] | -22.734 | 4.43 | +19.2 | -24.19 / -23.36 / -22.73 / -22.40 / -22.18 |
| rk | 6 | 256 | 2^22 | 66782 | -13.973 [-13.984, -13.962] | -13.963 | 17.57 | +187.1 | -14.76 / -14.45 / -14.08 / -13.54 / -12.92 |
| rk | 8 | 256 | 2^25 | 25324 | -18.372 [-18.390, -18.354] | -18.349 | 12.00 | +124.2 | -19.45 / -19.00 / -18.49 / -17.83 / -16.97 |
| rk | 10 | 64 | 2^29 | 5450 | -22.588 [-22.627, -22.550] | -22.734 | 14.60 | +76.3 | -23.75 / -23.39 / -22.73 / -21.89 / -21.39 |

Findings:
- Averaged over keys, both key schedules give the independent-round-key value (within 0.06 bits).
- For a fixed key the probability depends strongly on the key (D from 4 to 21). The spread is the same
  for fixed independent round keys (rk), so it comes from fixing the keys, not from the key schedule.
- The median key is between 0.04 bits above and 0.15 bits below the pooled value and the 10% quantile 0.5 to 0.8 bits below,
  at 6 to 10 rounds. How the spread behaves at 27 rounds is not measured; the data complexity of
  section 4 is for the average key.
- At 10 rounds each key has only about 75 right pairs, so part of the spread there is sampling noise.

## 9. Differential-linear (`08_*`, `08_table`; measured, independent round keys)
Input difference 0x800000; all 64 single-bit output masks; largest |correlation| per round:

| r | log2 pairs | best bit | corr | log2 abs | noise log2 (1/sqrt N) | DL data log2(1/c^2) | diff data log2(1/p) | fewer pairs |
|---|---|---|---|---|---|---|---|---|
| 1 | 26 | 63 | +1.00000000 | 0.000 | -13.00 | -0.00 | 3.00 | DL |
| 2 | 26 | 63 | +1.00000000 | 0.000 | -13.00 | -0.00 | 5.19 | DL |
| 3 | 28 | 43 | +0.91795719 | -0.124 | -14.00 | 0.25 | 7.39 | DL |
| 4 | 28 | 41 | +0.38305487 | -1.384 | -14.00 | 2.77 | 9.58 | DL |
| 5 | 30 | 11 | +0.09931290 | -3.332 | -15.00 | 6.66 | 11.77 | DL |
| 6 | 30 | 11 | +0.02192372 | -5.511 | -15.00 | 11.02 | 13.96 | DL |
| 7 | 32 | 50 | +0.00493563 | -7.663 | -16.00 | 15.33 | 16.16 | DL |
| 8 | 34 | 50 | +0.00108895 | -9.843 | -17.00 | 19.69 | 18.35 | differential |
| 9 | 35 | 50 | +0.00023237 | -12.071 | -17.50 | 24.14 | 20.54 | differential |
| 10 | 36 | 50 | +0.00005338 | -14.193 | -18.00 | 28.39 | 22.73 | differential |
| 11 | 37 | 11 | +0.00001255 | -16.282 | -18.50 | 32.56 | 24.93 | differential |

Decay of the best correlation: 2.155 bits/round over rounds 7..11 (2.146 over 8..11).
With the crude comparison 1/c^2 pairs (DL) against 1/p pairs (differential), the differential needs
fewer pairs from 8 rounds on; at 7 rounds the DL test still needs fewer (2^15.3 vs 2^16.2).
At 11 rounds the best correlation is 2.2 bits above the noise level, and the maximum over 64 bits is
biased upward near the noise level.

## 10. Search for better trails (`11_sat_*`, `11b_lin_trail`, `12_lowgraph_k2`)

- Standard (independence) model, CP-SAT, each r proven optimal (exact within the model; `11_sat_selftest`:
  the Lipmaa-Moriai and Wallen constraints match brute force on 5-bit words):
  - differential, r = 1..10: 2, 5, 9, 13.415, 19, 23, 27, 31, 35, 39. From r = 5 every optimal trail is
    0x800000 -> 0x800000/0x900000 (weight 4 per round in the model, 3 exactly).
  - linear, r = 1..10: 1, 2, 3, 6, 8, 11, 13, 16, 19, 21 (the published 10-round optimum is 21).
    Exact one-round correlations along the 10-round trail multiply to 2^-21.1943 (`11b_lin_trail`, exact).
  - the published single-trail weights (Table 16 of the Dipper paper) are 2, 5, 9, 13.42, 19, 23, 27.42 for
    r = 1..7; ours agree up to r = 6 and give 27 at r = 7. Exact one-round probabilities along our optimal
    7-round trail multiply to 2^-22.5081 (`11c_diff_trail`, exact), 4.5 bits above its model weight.
- The standard model is not an upper bound on exact differential probabilities (the S-box/addition
  dependence can raise them); its support contains every exactly possible transition.
- Exact graph of all 27240 differences with <= 2 active nibbles, all one-round transitions >= 2^-16
  (171032 edges; `12_lowgraph_k2`, exact for this class): minimum mean cycle = self-loop of 0x800000,
  3.0000 bits/round; spectral radius 2^-2.1926 per round (= the {0x800000, 0x900000} cluster).
  Differences with more than two active nibbles are not covered.

## 11. Evaluation of new rotation constants (`13_*`; r_A=1, r_B=4, r_D=11)

| r_C | 1-nibble graph cycle | 2-nibble cycle | 2-nibble cluster | SAT diff r=8 | SAT lin r=8 | best iterative linear (<= 2 nibbles) |
|---|---|---|---|---|---|---|
| 7 | 3.000 | 3.000 | 2.1926 | 31 | 16 | 2.75 |
| 3 | none | 5.1383 | 4.5130 | 36.415 | 17 | 3.75 |
| 9 | none | 5.0000 | 4.6365 | 37.415 | 16 | 3.75 |
| 11 | none | not computed | not computed | 35 | 16 | 2.75 |

Cycles and clusters in bits per round (exact for the class); SAT weights exact within the standard model.

## 12. Differences with three active nibbles (`16_*`, exact for the class)
`lowgraph3` computes each half distribution by conditioning on the values of the active nibbles
(`src/halffast.h`; at most 16^3 cases times 2^16) instead of enumerating 2^32 values. Checks: equal
to the 2^32 brute force on 100 random half differences with <= 3 active nibbles (`16_hfcheck`), and
`lowgraph3 -k 2` reproduces step 12 exactly (27240 nodes, 171032 edges, 9256 on cycles; `16_lowgraph3_k2`).
All 1917240 differences with <= 3 active nibbles, all one-round transitions >= 2^-16 (89173524 edges;
`16_lowgraph3_k3`):
- minimum mean cycle: the self-loop of 0x800000, 3.0000 bits/round; spectral radius 2^-2.1926 per round;
- best self-loop with exactly 3 active nibbles: 2^-10 (0x2400020000000000); best cycle among them
  9.0538 bits/round; their cluster decays by 6.9069 bits/round;
- lower bound from 0x800000 to 0x800000 in this graph: r = 13: 2^-29.3117, r = 27: 2^-60.0088,
  equal to the two-state bound (1/8)(7/32)^(r-1). The third nibble adds nothing measurable.
Not covered: four or more active nibbles, transitions below 2^-16.

## 13. Dipper-32: the full-round data complexity checked on the codebook (`15_*`)
Dipper-32 (`src/toy32.h`): 8-bit words, the Dipper S-box, rotations (1,4,7,3) = Dipper's (1,4,7,11)
mod 8, a GIFT-like 32-bit permutation, modular additions mod 2^8, a scaled key schedule
(Dipper-32/48). Of the 8 permutation variants with these rotations, the one used (pc=1, ps=2) keeps
Dipper's self-loop: 0x8000 -> 0x8000 with exact probability 2^-3 (`15_toysearch`). Model as in step 06
(`toygraph -k 3`, exact one-round transitions, <= 3 nibbles): cycle 3.0000, cluster 2.2382 bits/round.

Pr[0x8000 -> 0x8000] (`15_toy_table`, Table A): model (exact lower bound), independent round keys
(2^36 pairs, measured) and the exact value over all 2^31 pairs for 64 keys (pooled; quantiles over keys):

| d | model | indep | pooled 64 keys | 10% | median | 90% |
|---|---|---|---|---|---|---|
| 4 | -9.524 | -9.524 | -9.515 | -9.588 | -9.499 | -9.461 |
| 8 | -18.481 | -18.480 | -18.466 | -18.585 | -18.455 | -18.367 |
| 10 | -22.957 | -22.951 | -22.935 | -23.052 | -22.939 | -22.815 |
| 11 | -25.195 | -25.173 | -25.153 | -25.445 | -25.167 | -24.934 |
| 12 | -27.434 | -27.473 | -27.338 | -28.193 | -27.300 | -26.913 |

Distinguisher on R = d + 1 rounds (last round peeled), structures on nibble 3, 16 targets, success
>= 0.90 at false alarm <= 0.01 (Table B). "Measured" = exact over each key's codebook, structures
taken independently with probability N/2^32:

| R | model log2 N | tau | model success | measured mean success (64 keys) | min over keys | keys >= 0.9 |
|---|---|---|---|---|---|---|
| 11 | 20.125 | 2 | 0.901 | 0.902 | 0.869 | 41/64 |
| 12 | 22.812 | 3 | 0.911 | 0.912 | 0.886 | 55/64 |
| 13 | 25.438 | 5 | 0.910 | 0.910 | 0.861 | 49/64 |
| 14 | 29.125 | 27 | 0.934 | 0.941 | 0.817 | 53/64 |

R = 13 is the analogue of full Dipper (data 2^6.6 below the codebook; Dipper: 2^58.3, 2^5.7 below).
At every R the key-averaged success equals the model prediction within 0.01; single keys range down
to 0.82-0.89. The step-06 method therefore predicts the data of the full-round distinguisher of the
toy cipher correctly. For Dipper itself the 2^58.3 remains inferred.

## 14. Key recovery on 14 to 20 rounds (`14_krec_*`, `14_table`; measured)
Attack on R rounds (`krec`): structures of 16 chosen plaintexts on nibble 5; the differential
0x800000 -> {0x800000, 0x900000} over rounds 1..R-2; round R is peeled without a key; the survivors of a
key-independent 59-bit filter (Delta z = 0 outside the top 5 bits of the C2 word, z = P^-1 of the peeled
ciphertext) are used to rank all 2^32 guesses of the C2 and D words of P^-1(rk_R): a guess scores one
for each survivor whose recomputed Delta x_{R-2} is 0x800000 or 0x900000. The ranking by tables is
checked against the plain formula on 4097 keys in every run. One random master key per seed (seeds 1-4:
Dipper-64/96, 5-8: Dipper-64/128), real key schedule.

| R | data (CP) | keys | correct guess ranked first | bits of the 32 determined, min / median / max | right pairs, min / median / max |
|---|---|---|---|---|---|
| 14 | 2^30 | 8 | 8/8 | 15.48 / 18.00 / 21.42 | 13 / 21 / 31 |
| 15 | 2^32 | 8 | 8/8 | 13.42 / 17.75 / 19.23 | 13 / 22 / 37 |
| 16 | 2^34 | 8 | 8/8 | 14.75 / 16.77 / 18.42 | 11 / 19 / 21 |
| 17 | 2^36 | 8 | 8/8 | 13.57 / 15.84 / 18.91 | 6 / 16 / 23 |
| 18 | 2^38 | 8 | 8/8 | 12.91 / 15.80 / 20.00 | 7 / 13 / 17 |
| 19 | 2^40 | 4 | 4/4 | 12.98 / 15.54 / 18.96 | 7 / 12 / 13 |
| 20 | 2^42 | 2 | 2/2 | 15.42 / 18.22 (2 keys) | 6 / 17 (2 keys) |

- No survivor came from a pair outside nibble 5 in any run, and every survivor had Delta x_{R-2} in
  {0x800000, 0x900000}. So the filter passed no wrong pairs, and the correct guess always gets the
  maximum score. The score depends on the guess only through the top 5 bits of the C2 word and, per
  pair, through the top 5 bits of the D word plus one borrow from the low 11 bits; the guesses tied
  with the correct one are those that agree on these for every survivor.
- Time: 2^42 encryptions of 20 rounds took 2.5 to 4.0 hours on 8 cores (Contabo); the ranking takes seconds.
- Full key (inferred): the remaining kappa - b key bits by exhaustive search, b = bits determined; e.g.
  Dipper-64/96 at 20 rounds: about 2^(96-15.4) = 2^80.6 trial encryptions, data 2^42 CP. Only 32 bits of
  one round key are reachable this way (the iterative differential lives in the (C, D) half).

## 15. Key recovery on the full 28 rounds (`17_*`; inferred)
The attack of section 14 at R = 28 needs more data than can be run, so it is estimated in two parts.

Signal (`17_model`): model of section 4 with p(d) = average over the 15 nibble-5 input differences of
Pr[-> {0x800000, 0x900000}, d rounds]; log2 p(26) = -57.865. Test of the model at R = 14, 2^30 CP
(`17_keys_r14_s*`, 256 fresh master keys, measured): mean survivors 53.49 +- 1.17 s.e. against 53.40, dispersion
6.59. The 46 runs of section 14 give fewer (1344 against 1822.3), but they use the same 8 master keys for every R,
so they sample only 8 keys of a spread distribution.

Noise: a wrong pair passes the filter with probability 31 * 2^-64, so wrong pairs appear at this data.
`krsim` draws K ~ Poisson(signal) real right pairs (`17_pool_s*`, collected on 14 rounds, 490 to 1058 pairs per key)
and W ~ Poisson(noise) synthetic wrong pairs (z0 uniform, z1 = z0 ^ (delta << 27), delta in 1..31) and
ranks all 2^32 guesses. 8 keys, 40 trials per key and data volume (`17_sim_c*_s*`, `17_table`).
Check without noise at the signal of the 20-round runs (`17_check_s*`, K = 26): true guess first in 320/320, bits
14.54 / 16.41 / 18.43 (10% / median / 90%); measured at 20 rounds: 15.42, 18.22.

| data (CP) | K | W | true guess first | bits determined, 10% / median / 90% |
|---|---|---|---|---|
| 2^56 | 2.08 | 0.91 | 0.650 | 1.63 / 5.83 / 10.41 |
| 2^57 | 4.09 | 1.94 | 0.697 | 3.42 / 10.00 / 12.70 |
| 2^58 | 8.16 | 3.43 | 0.844 | 9.79 / 12.61 / 14.91 |
| 2^59 | 16.51 | 7.49 | 0.863 | 12.68 / 14.79 / 16.74 |
| 2^60 | 32.57 | 14.48 | 0.884 | 14.83 / 16.88 / 18.88 |
| 2^61 | 65.99 | 28.68 | 0.891 | 16.96 / 18.74 / 20.96 |

"Bits determined" here and in section 14 is the reduction of the 2^32 candidates, 32 - log2(guesses kept), and does
not mean that particular key bits are known. At 2^59 CP (2^5 below the codebook) the median reduction is 14.8 bits;
testing the remaining candidates costs about 2^(96-14.8) = 2^81.2 trial encryptions for Dipper-64/96 and 2^113.2 for Dipper-64/128.
These numbers are inferred: rounds 15 to 26 rest on the independent-round-key model, the wrong pairs are synthetic,
and the attack was not run on 28 rounds.

## 16. Data complexity with and without the background (`06c_background`; inferred, Dipper-32 measured)
Step 06 adds the random background |S| 2^-64 to p_S. Since p_S is a lower bound on the probability of reaching S, the
sum is justified only if the trails outside the cluster behave like the background. The mean count is at least
pairs * p_S, so the N computed from p_S alone is an upper bound within the independent-round-key model.

| setting | target set size | log2 N with background | log2 N, signal only |
|---|---|---|---|
| pairs (x, x^0x800000) | 1 | 63.375 | 63.438 |
| pairs (x, x^0x800000) | 16 | 61.062 | 61.312 |
| 16-plaintext structures on nibble 5 | 1 | 60.562 | 60.750 |
| 16-plaintext structures on nibble 5 | 16 | 58.312 | 59.062 |

Best setting, signal only: threshold 10, expected 3.9 for a random permutation and at least 14.4 for Dipper,
success 0.908. At 2^58.312 the signal alone gives mean 8.561 and success 0.750 instead of 0.917.

Test of the assumption on Dipper-32 (full codebook, 64 keys, section 13):

| R | lam0 | lam1 with background | measured mean count | success with background | signal only | measured |
|---|---|---|---|---|---|---|
| 11 | 0.032 | 3.900 | 3.919 | 0.901 | 0.898 | 0.902 |
| 12 | 0.206 | 5.487 | 5.514 | 0.911 | 0.897 | 0.912 |
| 13 | 1.270 | 8.175 | 8.179 | 0.910 | 0.818 | 0.910 |
| 14 | 16.358 | 35.213 | 35.699 | 0.934 | 0.045 | 0.941 |

The assumption holds on Dipper-32. This supports 2^58.3 for Dipper but does not prove it.
