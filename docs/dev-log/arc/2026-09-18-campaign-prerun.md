# Campaign pre-run (arc B): with/without GNN vs Rphylopars and BACE

2026-09-18, Totoro, branch feat/gnn-off (1458f34), `script/campaign_gnn_off_cell.R`.
DGP `bm_mixed`: 4 BM continuous + 1 binary + 1 three-level categorical on `ape::rcoal(n)`, 30% user-level MCAR, same mask across arms.
Cells: n in {100, 1000} x seeds {1, 2}; all seven arms; GNN-on at epochs = 2000; BACE nitt 50,000 / burnin 10,000 / thin 25, OVR on, 2 chains, 5 final draws.
4 threads per cell (`OMP_NUM_THREADS=4`, `OPENBLAS_NUM_THREADS=1`), 4 cells in parallel, cap 40 min per cell. No arm errored; no cell hit the cap.
**Two seeds: a pre-run to ground the estimate and prove the invocation, not evidence about the methods.**

## Wall time per arm (seconds, mean of 2 seeds)

| arm | n = 100 | n = 1000 |
|---|---:|---:|
| bace | 152.6 | 1381.0 |
| gnn_on | 136.2 | 523.2 |
| gnn_off | 0.6 | 135.9 |
| gnn_off_pure | 0.2 | 2.2 |
| rphylopars | 0.4 | 1.9 |
| floor | 0.0 | 0.0 |
| gnn_on_full | 0.0 | 0.0 |

`gnn_on_full` is derived from the `gnn_on` fit (predict with `baseline_override`), so its wall is 0. The 136 s of `gnn_off` at n = 1000 is the gate calibration (`calibrate_gates()`, cv folds over the simplex grid), not the baseline: `gnn_off_pure` skips it and takes 2 s.

## Mean z-RMSE over the 4 continuous traits (mean of 2 seeds; sd across seeds in parentheses)

| arm | n = 100 | n = 1000 |
|---|---:|---:|
| rphylopars | 0.292 (0.011) | 0.090 (0.012) |
| gnn_off_pure | 0.295 (0.026) | 0.090 (0.011) |
| gnn_off | 0.296 (0.026) | 0.090 (0.011) |
| gnn_on_full | 0.307 (0.034) | 0.090 (0.011) |
| gnn_on | 0.320 (0.038) | 0.099 (0.012) |
| bace | 0.354 (0.031) | 0.101 (0.007) |
| floor | 1.101 (0.004) | 0.994 (0.031) |

## Mean accuracy over the 2 discrete traits, and production-interval coverage (nominal 0.95; pigauto arms only)

| arm | acc n=100 | acc n=1000 | coverage n=100 | coverage n=1000 |
|---|---:|---:|---:|---:|
| rphylopars |  |  |  |  |
| gnn_off_pure | 0.742 | 0.952 | 0.883 | 0.954 |
| gnn_off | 0.742 | 0.953 | 0.883 | 0.954 |
| gnn_on_full | 0.742 | 0.953 | 0.838 | 0.951 |
| gnn_on | 0.775 | 0.950 | 0.829 | 0.941 |
| bace | 0.767 | 0.946 |  |  |
| floor | 0.333 | 0.393 |  |  |

Baseline dispatch recorded by `fit$baseline$path` for the GNN-off arms: continuous and binary columns via `threshold_joint`, the categorical trait via `ovr_categorical` (both n).

## What the pre-run says (regime: BM-correct mixed-type DGP, 2 seeds)

- Under a BM-correct DGP the GNN-off arm and raw Rphylopars are within 0.005 z-RMSE of each other at both n, as the design predicts (the continuous columns come from the same joint BM fit; the threshold-joint liability step costs nothing visible here).
- GNN-on is slightly worse than GNN-off on the continuous traits (0.320 vs 0.296 at n = 100; 0.099 vs 0.090 at n = 1000); `gnn_on_full` closes about half of that gap at n = 100 and all of it at n = 1000, which is the held-out-cell tax the 2026-08-16 note suspected. Two seeds; the campaign decides.
- BACE at these settings trails every pigauto arm on the continuous traits and matches them on the discrete ones; it is 9x slower than GNN-on at n = 1000.
- Coverage of the GNN-off production interval: 0.88 at n = 100 (36 masked cells per trait), 0.95 at n = 1000.

## Re-stated estimate for the full campaign (D-139)

Per-cell wall with arms run sequentially, 4 threads: n = 100 about 290 s (measured), n = 1000 about 2,050 s (measured), n = 300 about 800 s (interpolated from the 2026-08-16 GNN-on 362 s plus BACE scaling).
Locked design: 3 DGPs x 3 n x 20 seeds = 180 cells, plus AVONET300 x 20 seeds = 20 cells at n = 300. Cell-time about 60 x 290 + 80 x 800 + 60 x 2,050 = 204,000 s = 57 h.
On Totoro at 36 concurrent cells (144 threads, under the 150-core cap): about 1.6 h wall; call it 2 to 3 h with the OU and BACE DGPs, which the pre-run did not time. BACE dominates: 23 min of each n = 1000 cell.
Options that shrink it without a strawman: BACE only at n <= 300 (removes 60 x 1,380 s, wall about 1 h); or 10 seeds for BACE at n = 1000.

This exceeds 30 minutes, so the full submit waits for Shinichi (plan: MUST STOP). Command, once approved, runs from `~/gnn-off` on Totoro with the same runner; results land in `~/gnn-off/results_full/`.

Raw cell rds: `script/campaign_gnn_off_prerun/` (4 files).
