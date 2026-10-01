# Pre-run 1: does exact (cross-trait) prediction close the lambda 0.3 gap?

Run on Totoro 2026-09-24, pigauto main at 2a59dc1 (lambda estimated by default, before the REML fix).
Driver: `exact_prerun_cell.R` (this folder), which reuses the imputation-sim lane's `make_cell()`,
`run_arms()` and `score_arm()` unchanged. Cells: types_mixed, BM, rho 0, MCAR 30%, driver on, fixed
thresholds; lambda {0.3, 0.7} x n {100, 300, 1000}; seeds 1 to 20. Arms: pigauto GNN off (default
per-column prediction), the same with `predict_method = "exact"`, and the frequentist lambda reference.
120 of 120 jobs, 0 errors. Mean z-RMSE over c1, c2, cnt, prp:

| true lambda | n | frequentist lambda | default | exact | exact vs default | gap closed |
|---|---|---|---|---|---|---|
| 0.3 | 100 | 0.945 | 0.980 | 0.973 | -0.8% | 21% |
| 0.3 | 300 | 0.914 | 0.976 | 0.948 | -2.9% | 45% |
| 0.3 | 1000 | 0.893 | 0.947 | 0.910 | -3.9% | 69% |
| 0.7 | 100 | 0.858 | 0.875 | 0.844 | -3.5% | beats reference |
| 0.7 | 300 | 0.788 | 0.802 | 0.764 | -4.7% | beats reference |
| 0.7 | 1000 | 0.789 | 0.778 | 0.732 | -5.9% | beats reference |

Discrete accuracy: default 0.457 / 0.604, exact 0.466 / 0.631 at lambda 0.3 / 0.7. Coverage (c1, c2):
default 0.931 / 0.926, exact 0.938 / 0.931. Median time per job: 13 s for default plus the frequentist
arm, 6.7 s for the exact fit.

Reading: the lambda 0.3 gap is mostly cross-trait information, not lambda bias (the REML fix, PR #191,
cut the bias but moves predictions by well under 1%). Twenty seeds per cell is a pre-run, not evidence
to report; the full benchmark follows the centring fix and the default flip.
