# Conformal-coverage re-measurement -- results

Run: 2026-05-21 10:28:14.  Source: `dev/simulation_results_coverage_remeasure/results.rds`.
Design: follow-up to `script/sim_bench/mechanism_sweep_methodology.md`, with the same mechanism families, 120 replicates, n=500, 30% missing.

`cov_conformal` = empirical coverage of pigauto's actual split-conformal interval (PRIMARY). `cov_drawband` = the mixed-type sweep's old metric (2.5/97.5 band of 20 draws) -- difference from cov_conformal is the measurement artifact. `cov_bm` = BM-kriging analytic interval (reference). Nominal target = 0.95. Values are mean +/- Monte Carlo SE.

## Coverage by mechanism

| scenario | mechanism | cov_conformal | cov_drawband | cov_bm | conf_width |
|---|---|---|---|---|---|
| bm_strong | MCAR | 0.937 +/- 0.015 | 0.872 +/- 0.020 | 0.931 +/- 0.007 | 1.60 |
| bm_strong | phylo_MAR | 0.952 +/- 0.010 | 0.890 +/- 0.010 | 0.949 +/- 0.005 | 1.63 |
| bm_strong | trait_MAR | 0.934 +/- 0.018 | 0.853 +/- 0.024 | 0.924 +/- 0.008 | 1.60 |
| bm_strong | trait_MNAR | 0.878 +/- 0.022 | 0.799 +/- 0.021 | 0.928 +/- 0.008 | 1.50 |
| weak_signal | MCAR | 0.966 +/- 0.007 | 0.898 +/- 0.014 | 0.918 +/- 0.010 | 4.23 |
| weak_signal | phylo_MAR | 0.948 +/- 0.012 | 0.876 +/- 0.017 | 0.954 +/- 0.009 | 4.06 |
| weak_signal | trait_MAR | 0.951 +/- 0.008 | 0.893 +/- 0.013 | 0.933 +/- 0.010 | 4.15 |
| weak_signal | trait_MNAR | 0.815 +/- 0.020 | 0.690 +/- 0.022 | 0.890 +/- 0.014 | 3.64 |

## Reading

- **MCAR** is the reference: split conformal's exchangeability assumption holds, so `cov_conformal` should sit near 0.95.
- A `cov_conformal` materially below 0.95 under MAR/MNAR is genuine under-coverage (exchangeability between observed calibration cells and missing target cells is broken).
- `cov_conformal` minus `cov_drawband` is the artifact: how much the mixed-type sweep's 20-draw-band metric understated true coverage.
