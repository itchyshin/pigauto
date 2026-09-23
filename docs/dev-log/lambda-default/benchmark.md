# Core simulation benchmark under the lambda default

18 core cells of the four-arm simulation (types_mixed, BM, lambda {0.3, 0.7, 1}, rho {0, 0.5}, MCAR 30%,
n {100, 300, 1000}), pigauto 0.11.0.9000 from `~/R/lib-lambda` on Totoro, campaign scripts unmodified.
Committed yardstick: `script/campaign_sim_results/summary.csv` on arc/imputation-sim (#184). Comparison
script: `compare_benchmark.R` (mean z-RMSE over c1, c2, cnt, prp and rho).

## Wave 1: no-GNN arms (gnn_off, gnn_off_rphylopars), 200 seeds per cell

Launched 05:52 MDT 2026-09-23, 140 parallel, done 07:38. 3,595 of 3,600 jobs; the 5 failures are one
error in the campaign's F1 helper (`if (tp + fp == 0)`), outside the package.

In-house pigauto without the GNN (the default path):

| lambda | n | committed | lambda default | change | gap to freq_lambda closed |
|---|---|---|---|---|---|
| 0.3 | 100 | 1.019 | 0.998 | -0.021 | 17% |
| 0.3 | 300 | 1.012 | 0.977 | -0.034 | 28% |
| 0.3 | 1000 | 0.990 | 0.965 | -0.025 | 20% |
| 0.7 | 100 | 0.916 | 0.848 | -0.068 | 59% |
| 0.7 | 300 | 0.888 | 0.806 | -0.082 | 75% |
| 0.7 | 1000 | 0.873 | 0.787 | -0.086 | 81% |
| 1 | 100 | 0.491 | 0.466 | -0.025 | |
| 1 | 300 | 0.426 | 0.393 | -0.033 | |
| 1 | 1000 | 0.363 | 0.342 | -0.021 | |

Every cell improves. Conformal coverage for c1 and c2 at n >= 300 is 0.957 to 0.960 (committed 0.950
to 0.968). No seed exceeded z-RMSE 1.3.

Against the pre-registered gate G12:

- Gap closed at lambda 0.3, n 100 is 17%, below the 50% the gate asked for. At lambda 0.7 it is 59 to
  81%. Weak signal at lambda 0.3 is where the per-column estimate's downward bias and small n bite most.
- The lambda = 1 cells moved by 0.021 to 0.033, beyond the 0.01 limit, in the improving direction. The
  guard was written to catch regressions; none occurred.
- Coverage (>= 0.94) and completeness pass.

Rphylopars comparator (`joint_solver = "rphylopars"`, now `model = "lambda"`): medians improve (for
example 0.824 at lambda 0.7, n 100), but a few seeds explode (1 of 200 at lambda 0.7, n 100; 6 to 8 of
about 198 at lambda 1, n 1000) with z-RMSE up to 10^5, so its means are not interpretable. Fixed after
the run: `fit_joint_solver()` now falls back to the in-house solver when phylopars returns tip
predictions more than 10 times beyond the observed range (test in `test-joint-lambda.R`). The in-house
default path was never affected. This arm is 26 times slower under lambda (45 s against 1.7 s per fit).

## Wave 2: gnn_on, 100 seeds per cell

Launched 07:38, 30 parallel x 4 torch threads. Capped at 100 seeds because the full 200 would take
about 8 hours at the 150-core limit. Results appended here when it finishes.
