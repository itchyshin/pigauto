# BACE multi-obs predictor-phylo-signal diagnostic

Generated: 2026-04-26 11:03:08

Sweep: predictor_phylo_signal in {0.0, 0.3, 0.6} x beta in {0, 0.5} x 2 reps.
Held constant: n_species=100, multi_obs_ratio=4, response_phylo_signal=0.4,
response_type=gaussian, n_predictors=3 (gaus, gaus, bin), miss_frac=0.30.

## Mean RMSE per cell (lower is better)

| beta | pred_phylo | column_mean | pigauto_no_cov | pigauto_cov |
|------|-----------|-------------|----------------|-------------|
| 0.0 | 0.0 | 1.7104 | 1.1090 | 1.1215 |
| 0.0 | 0.3 | 1.6237 | 1.1949 | 1.2417 |
| 0.0 | 0.6 | 1.7255 | 1.2582 | 1.2433 |
| 0.5 | 0.0 | 1.6344 | 1.2494 | 1.2753 |
| 0.5 | 0.3 | 1.7273 | 1.3340 | 1.3300 |
| 0.5 | 0.6 | 1.6667 | 1.2767 | 1.2739 |

## Diagnostic ratios

Ratios are pigauto / column_mean. < 1 means pigauto helps.

| beta | pred_phylo | no_cov / col_mean | cov / col_mean |
|------|-----------|-------------------|---------------|
| 0.0 | 0.0 | 0.648 | 0.656 |
| 0.0 | 0.3 | 0.736 | 0.765 |
| 0.0 | 0.6 | 0.729 | 0.721 |
| 0.5 | 0.0 | 0.764 | 0.780 |
| 0.5 | 0.3 | 0.772 | 0.770 |
| 0.5 | 0.6 | 0.766 | 0.764 |
