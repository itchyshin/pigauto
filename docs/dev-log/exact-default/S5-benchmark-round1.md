# S5 round 1: exact default against current main

Totoro 2026-09-25. Main = origin/main eccb298 (lambda estimated, full REML) in ~/R/lib-main; exact =
feat/exact-default after S3 (exact default, GLS-mean centring) in ~/R/lib-exact3. Simulation: 18 core cells,
pigauto GNN off, 200 seeds, unmodified campaign scripts. Real data: `script/bench_lambda_datasets.R`,
13 cases x 5 masks x 4 lambda modes; the table uses lambda_mode = "estimate".

## Simulation (mean over c1, c2, cnt, prp and rho)

| true lambda | n | z-RMSE main | exact | change | discrete accuracy change |
|---|---|---|---|---|---|
| 0.3 | 100 | 0.998 | 0.961 | -3.6% | +0.013 |
| 0.3 | 300 | 0.978 | 0.931 | -4.8% | +0.019 |
| 0.3 | 1000 | 0.965 | 0.899 | -6.8% | +0.028 |
| 0.7 | 100 | 0.847 | 0.786 | -7.2% | +0.023 |
| 0.7 | 300 | 0.805 | 0.732 | -9.1% | +0.037 |
| 0.7 | 1000 | 0.787 | 0.708 | -10.0% | +0.046 |
| 1 | 100 | 0.465 | 0.465 | 0.0% | +0.001 |
| 1 | 300 | 0.393 | 0.402 | +2.0% | -0.001 |
| 1 | 1000 | 0.341 | 0.350 | +2.6% | -0.002 |

Coverage (c1, c2) rises in every cell (for example 0.958 to 0.962 at lambda 0.3, n 300). Binary and ordinal
gain most (binary 0.755 to 0.796, ordinal 0.490 to 0.553 at lambda 0.7); categorical is unchanged.

## Real data (mean z-RMSE, lambda estimated)

| dataset | species | main | exact | change |
|---|---|---|---|---|
| AVONET 300 (bundled) | 300 | 0.555 | 0.275 | -50.5% |
| AVONET | 300 | 0.627 | 0.402 | -35.8% |
| AVONET | 2000 | 0.419 | 0.281 | -33.0% |
| PanTHERIA | 300 | 0.496 | 0.402 | -18.9% |
| PanTHERIA | 2000 | 0.358 | 0.297 | -17.1% |
| AmphiBIO | 2000 | 0.665 | 0.644 | -3.1% |
| BIEN | 2000 | 0.807 | 0.800 | -0.8% |
| BIEN | 300 | 1.051 | 1.043 | -0.8% |
| LepTraits | 2000 | 0.901 | 0.895 | -0.6% |
| GlobTherm | 300 | 0.665 | 0.663 | -0.4% |
| AmphiBIO | 300 | 0.875 | 0.879 | +0.5% |
| GlobTherm | 1969 | 0.622 | 0.625 | +0.5% |
| LepTraits | 300 | 0.935 | 0.940 | +0.5% |

Discrete accuracy mostly rises (AVONET 300 migration 0.813 to 0.817, primary lifestyle 0.753 to 0.770,
PanTHERIA 300 terrestriality 0.850 to 0.871); AmphiBIO 300 habitat falls 0.794 to 0.784 and PanTHERIA 2000
terrestriality 0.936 to 0.934.

## Gate

The pre-registered gate (every cell and dataset not worse than main) fails at lambda = 1, n >= 300 (+2.0% and
+2.6%), on three real-data cases (+0.5% each) and on AmphiBIO 300 habitat accuracy (-0.010). Likely cause: the
Kronecker solve uses one shared lambda, and the count and proportion traits pull it below 1, over-shrinking
continuous traits whose true lambda is 1. Decision (Shinichi, 2026-09-25): add a per-trait choice between the
exact and per-column predictions on the validation split, then re-run the benchmark.
