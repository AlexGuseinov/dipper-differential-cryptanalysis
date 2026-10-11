# Commands

Every number in RESULTS.md and in the paper comes from one of the commands below.
Each line: command, raw output file, start time (UTC), wall time.
Build used for the measurements: `make NATIVE=1` unless stated otherwise.

- `python3 steps/02_one_round_exact.py` -> `results/raw/02_one_round_exact.txt` (2026-10-11T10:51:39Z, 64s)
- `make test NATIVE=1` -> `results/raw/00_make_test.txt` (2026-10-11T10:52:51Z, 5s)
- `./bin/dltool diff -r 1 -i 800000 -o 800000,900000,2000800000,2000900000 -n 28 -k indep -s 41 -t 4` -> `results/raw/04_iter_r1.txt` (2026-10-11T10:53:06Z, 0s)
- `./bin/dltool diff -r 2 -i 800000 -o 800000,900000,2000800000,2000900000 -n 28 -k indep -s 42 -t 4` -> `results/raw/04_iter_r2.txt` (2026-10-11T10:53:06Z, 1s)
- `./bin/dltool diff -r 3 -i 800000 -o 800000,900000,2000800000,2000900000 -n 28 -k indep -s 43 -t 4` -> `results/raw/04_iter_r3.txt` (2026-10-11T10:53:07Z, 1s)
- `./bin/dltool diff -r 4 -i 800000 -o 800000,900000,2000800000,2000900000 -n 30 -k indep -s 44 -t 4` -> `results/raw/04_iter_r4.txt` (2026-10-11T10:53:08Z, 3s)
- `./bin/dltool diff -r 5 -i 800000 -o 800000,900000,2000800000,2000900000 -n 30 -k indep -s 45 -t 4` -> `results/raw/04_iter_r5.txt` (2026-10-11T10:53:11Z, 4s)
- `./bin/dltool diff -r 6 -i 800000 -o 800000,900000,2000800000,2000900000 -n 30 -k indep -s 46 -t 4` -> `results/raw/04_iter_r6.txt` (2026-10-11T10:53:15Z, 5s)
- `./bin/dltool diff -r 7 -i 800000 -o 800000,900000,2000800000,2000900000 -n 30 -k indep -s 47 -t 4` -> `results/raw/04_iter_r7.txt` (2026-10-11T10:53:20Z, 6s)
- `./bin/dltool diff -r 8 -i 800000 -o 800000,900000,2000800000,2000900000 -n 30 -k indep -s 48 -t 4` -> `results/raw/04_iter_r8.txt` (2026-10-11T10:53:26Z, 6s)
- `./bin/dltool diff -r 9 -i 800000 -o 800000,900000,2000800000,2000900000 -n 32 -k indep -s 49 -t 4` -> `results/raw/04_iter_r9.txt` (2026-10-11T10:53:32Z, 27s)
- `./bin/dltool diff -r 10 -i 800000 -o 800000,900000,2000800000,2000900000 -n 32 -k indep -s 410 -t 4` -> `results/raw/04_iter_r10.txt` (2026-10-11T10:53:59Z, 30s)
- `./bin/dltool diff -r 11 -i 800000 -o 800000,900000,2000800000,2000900000 -n 34 -k indep -s 411 -t 4` -> `results/raw/04_iter_r11.txt` (2026-10-11T10:54:29Z, 145s)
- `python3 steps/03_bound.py` -> `results/raw/03_bound.txt` (2026-10-11T11:03:42Z, 1s)
- `./bin/dltool diff -r 12 -i 800000 -o 800000,900000,2000800000,2000900000 -n 36 -k indep -s 412 -t 4` -> `results/raw/04_iter_r12.txt` (2026-10-11T10:56:54Z, 574s)
- `./bin/lintool search -k 2 -m -5 -L 4` -> `results/raw/09_lin_search_k2.txt` (2026-10-11T11:01:15Z, 664s)
- `./bin/onestep cluster -i 800000 -m -7 -N 64 -r 27 -f results/data/cluster_800000.txt` -> `results/raw/06_cluster.txt` (2026-10-11T11:15:38Z, 308s)
- `python3 steps/06_data_complexity.py` -> `results/raw/06_data.txt` (2026-10-11T11:20:57Z, 119s)
- `python3 steps/06_data_complexity.py --validate` -> `results/raw/06_validate_model.txt` (2026-10-11T11:23:35Z, 32s)
- `./bin/dltool diff -r 13 -i 800000 -o 800000,900000,2000800000,2000900000 -n 37 -k indep -s 413 -t 4` -> `results/raw/04_iter_r13.txt` (2026-10-11T11:06:28Z, 1417s)
