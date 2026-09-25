# Evidence: G3 recovery and posterior calibration (2026-09-24)

Code: commit 7a0f470 (`git archive`), run on Totoro under `/home/snakagaw/pigauto_mi_posterior/7a0f47030f/`.
DGP: `script/mi_gls/gate_recovery.R` (n = 1000, K = 2, 25% MCAR per trait, `ape::rtree`).
Sampler defaults: 4 chains, 1,000 burn-in plus 5,000 sweeps each, thin 20.

| File | What it holds |
|---|---|
| `g3_official.log` | The G3 gate as written (settings A to D, seeds 1 to 3). 12 of 12 fits converged. The gate FAILED on B seed 2 and C seed 2; REML gives the same values on those datasets. |
| `g3_seeds4to10.log` | The same gate script with seeds 4 to 10: 28 of 28 converged. |
| `g3_calibration.R`, `g3_calib.log`, `g3_calib_BCD_101_150.rds` | Calibration check: settings B, C and D, seeds 101 to 150 (150 fits, all converged). Records 95% credible-interval coverage of the true lambda_k and phylogenetic correlation, and REML lambda. |

Calibration (`g3_calib.log`, 50 fits per setting; coverage MCSE about 0.03):

| Setting | lambda_1 cover | lambda_2 cover | rho_P cover | SD of posterior means vs mean posterior SD (lambda_1; lambda_2) |
|---|---|---|---|---|
| B (0.5, 0.5; 0.7) | 0.92 | 0.98 | 0.98 | 0.065 vs 0.057; 0.055 vs 0.058 |
| C (0.3, 0.9; 0) | 0.94 | 0.92 | 0.92 | 0.068 vs 0.063; 0.024 vs 0.022 |
| D (0.05, 0.95; 0.7) | 0.98 | 0.92 | 0.98 | 0.032 vs 0.033; 0.014 vs 0.014 |

Reading: the posterior is calibrated for lambda and the phylogenetic correlation at n = 1000. G3's fixed 0.1 tolerance is about 1.7 posterior SD at lambda = 0.5, so some honest fits fall outside it.

Timing on Totoro at n = 1000, K = 2: 206 to 220 s per fit at the defaults, 8.6 to 9.2 ms per sweep.

Incident: the first launch of the calibration run exported no thread caps. Setting them inside R is too late for BLAS, and the load average reached about 18,000 for about 10 minutes. The run was killed and relaunched with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` exported at launch, giving 76 processes and 77 threads.

## G3 gate evidence at the campaign commit (2026-09-24)

`g3_calibration.rds`, `g3_calibration_run.log`: `script/mi_gls/gate_calibration.R run` at commit 69670d4
on Totoro. Settings A to D, seeds 101 to 150 (200 fits, all converged), 10 cores, 7,566 s wall under
full machine load. `Rscript script/mi_gls/gate_calibration.R check` prints `RECOVERY_OK`.

- A (lambda 1, 1): bias -0.004 and -0.003; rho_P coverage 0.96.
- B to D reproduce the 7a0f470 numbers above exactly (same seeds), which confirms the sampler did not
  change between the two commits.
- D's rho_P (bias -0.14) is not gated, because lambda_1 = 0.05 leaves almost no phylogenetic
  variance in trait 1 to identify the phylogenetic correlation.

## Round-2 merge and byte-identity (2026-09-24)

- `summarise_round2.log`: a re-run of `script/mi_gls/03_summarise_v2.R` over the 69670d4 and 9597e18
  campaign folders on Totoro, with its console lines (SOURCE: 4,800 files, then 3,258 replacing 58;
  `SETTINGS 0 of 8000`; `MIXED_CODE_SHA`). Its three outputs have the same md5 as the committed
  `sim_summary.csv`, `sim_summary.md` and `cell_coverage.csv`.
- `smoke_round2/`: the round-2 smoke run and the byte-identity checks, including the 14 converged
  cells the round-2 campaign re-ran by accident (`smoke_round2/README.md`).
