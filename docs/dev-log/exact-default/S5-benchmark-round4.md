# S5 round 4: "auto" with the size-aware validation split (the version proposed for release)

Same design as rounds 1 to 3: 18 core cells, pigauto GNN off, 200 seeds, 13 real-data cases; main =
origin/main eccb298. Under "auto", a trait's validation cells are split between route choice and
calibration only when each half keeps at least 19 cells (the fewest for which a 95% split-conformal
interval can reach nominal coverage); below that, all cells do both (as in round 2). Round 3 split every
trait with at least 10 cells and lost up to 0.032 coverage at n <= 300, because the calibration half became
too small.

## Simulation, change against main

| true lambda | z-RMSE, n 100 / 300 / 1000 | coverage (c1, c2) | discrete accuracy |
|---|---|---|---|
| 0.3 | -3.1% / -4.6% / -6.6% | +0.003 / +0.004 / +0.009 | +0.009 / +0.014 / +0.022 |
| 0.7 | -6.1% / -8.5% / -9.4% | +0.008 / +0.005 / +0.010 | +0.012 / +0.025 / +0.032 |
| 1 | -1.2% / -0.3% / +0.4% | -0.007 / -0.003 / +0.002 | -0.002 / -0.002 / -0.001 |

Per trait, 10 of 36 trait-cells have coverage below main; the largest falls are the count (-0.017) and proportion (-0.013) traits at lambda = 1, n = 100,
where validation cells are too few to split and the same cells choose the route and calibrate.

Coverage and accuracy columns are means over the continuous traits (c1, c2) and the discrete traits (bin,
cat3, ord). Per scenario, the largest accuracy falls are for cat3 (at most -0.009, unpaired |z| < 1.5, so
within Monte Carlo error). The ordinal loss that auto showed in round 2 at lambda 1, n 1000 (-0.006) is
-0.001 here; the planned ordinal scoring fix (S8) is still worth doing but no longer carries a measured loss.

## Real data (mean z-RMSE, lambda estimated)

| dataset | species | main | auto | change |
|---|---|---|---|---|
| AVONET 300 (bundled) | 300 | 0.555 | 0.275 | -50.5% |
| AVONET | 300 | 0.627 | 0.402 | -35.8% |
| AVONET | 2000 | 0.419 | 0.281 | -33.0% |
| PanTHERIA | 300 | 0.496 | 0.410 | -17.4% |
| PanTHERIA | 2000 | 0.358 | 0.298 | -16.8% |
| AmphiBIO | 2000 | 0.665 | 0.649 | -2.5% |
| BIEN | 2000 | 0.807 | 0.799 | -0.9% |
| BIEN | 300 | 1.051 | 1.047 | -0.4% |
| LepTraits | 2000 | 0.901 | 0.897 | -0.4% |
| AmphiBIO | 300 | 0.875 | 0.873 | -0.2% |
| GlobTherm | 300 | 0.665 | 0.664 | -0.2% |
| LepTraits | 300 | 0.935 | 0.934 | -0.1% |
| GlobTherm | 1969 | 0.622 | 0.627 | +0.8% |

Round 3 (split every trait with at least 10 validation cells) is kept for the record in
`core_lambda_auto3_off_agg_summary.csv` and `lambda_datasets_auto3_*`; its coverage fell by up to 0.032 at
n <= 300, which is why this version splits only at 38 cells or more.

Checks on this code: suite FAIL 0 / PASS 2724 (`S5d-suite.log`); R CMD check in `S5d-check.log`.
