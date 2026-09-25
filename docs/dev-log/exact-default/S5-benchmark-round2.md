# S5 round 2: the "auto" default against current main

Same design as round 1 (`S5-benchmark-round1.md`): 18 core cells, pigauto GNN off, 200 seeds, and 13
real-data cases; main = origin/main eccb298. The branch build now defaults to predict_method = "auto": both
the exact cross-trait route and the per-column route are fitted, and each trait takes the one with lower
validation loss (squared error on the latent scale; log-loss for binary and categorical).

## Simulation (mean z-RMSE over c1, c2, cnt, prp and rho)

| true lambda | n | main | exact (round 1) | auto | auto vs main | discrete accuracy change |
|---|---|---|---|---|---|---|
| 0.3 | 100 | 0.998 | 0.961 | 0.967 | -3.1% | +0.011 |
| 0.3 | 300 | 0.978 | 0.931 | 0.933 | -4.6% | +0.016 |
| 0.3 | 1000 | 0.965 | 0.899 | 0.900 | -6.8% | +0.027 |
| 0.7 | 100 | 0.847 | 0.786 | 0.795 | -6.1% | +0.019 |
| 0.7 | 300 | 0.805 | 0.732 | 0.737 | -8.5% | +0.033 |
| 0.7 | 1000 | 0.787 | 0.708 | 0.710 | -9.8% | +0.044 |
| 1 | 100 | 0.465 | 0.465 | 0.459 | -1.2% | 0.000 |
| 1 | 300 | 0.393 | 0.402 | 0.392 | -0.3% | -0.001 |
| 1 | 1000 | 0.341 | 0.350 | 0.341 | -0.1% | -0.002 |

Every cell is at or below main in z-RMSE; the round-1 losses at lambda = 1 are gone. Coverage is NOT uniformly
at or above main (an earlier version of this report said it was; that was a misreading, corrected after
review). Mean coverage for c1 and c2:

| true lambda | n 100 main / auto | n 300 main / auto | n 1000 main / auto |
|---|---|---|---|
| 0.3 | 0.877 / 0.880 | 0.958 / 0.962 | 0.956 / 0.964 |
| 0.7 | 0.861 / 0.869 | 0.960 / 0.965 | 0.957 / 0.966 |
| 1 | 0.900 / 0.894 | 0.969 / 0.966 | 0.968 / 0.968 |

At lambda = 1 coverage falls slightly (per trait, 27 of 72 trait-cells are lower; for example c1 0.889 to 0.883
and count 0.876 to 0.858 at n = 100). Round 1 (exact everywhere) did not lower coverage, which points to the
route choice and the conformal scores using the same validation cells; the fix (split the validation cells)
is in progress and this benchmark will be re-run. The discrete change at lambda = 1, n = 1000 comes from ordinal accuracy, about -0.006
(roughly 3 unpaired standard errors); binary and categorical are unchanged there. Ordinal gains about +0.06
at lambda 0.7. Likely cause: the route choice scores ordinal traits by squared error rather than class
accuracy. Follow-up agreed with Shinichi: score ordinal by class accuracy in the next change.

## Real data (mean z-RMSE, lambda estimated)

| dataset | species | main | auto | change |
|---|---|---|---|---|
| AVONET 300 (bundled) | 300 | 0.555 | 0.275 | -50.5% |
| AVONET | 300 | 0.627 | 0.402 | -35.8% |
| AVONET | 2000 | 0.419 | 0.281 | -33.0% |
| PanTHERIA | 300 | 0.496 | 0.410 | -17.4% |
| PanTHERIA | 2000 | 0.358 | 0.298 | -16.8% |
| AmphiBIO | 2000 | 0.665 | 0.648 | -2.6% |
| LepTraits | 2000 | 0.901 | 0.893 | -0.9% |
| BIEN | 2000 | 0.807 | 0.801 | -0.7% |
| BIEN | 300 | 1.051 | 1.047 | -0.4% |
| AmphiBIO | 300 | 0.875 | 0.873 | -0.2% |
| GlobTherm | 300 | 0.665 | 0.664 | -0.2% |
| LepTraits | 300 | 0.935 | 0.934 | -0.1% |
| GlobTherm | 1969 | 0.622 | 0.622 | +0.1% |

Discrete accuracy is equal to or above main in every dataset (AmphiBIO 300 habitat 0.794, as on main).

Decision (Shinichi, 2026-09-25): ship "auto" with the one ordinal cell documented; then, after review found the
coverage drop, split the validation cells between route choice and conformal calibration and re-benchmark.
Checks on this code: suite FAIL 0 / PASS 2663; R CMD check 0 errors, 0 warnings, 1 note (dev version).
