# Ordinal route choice: two alternative scorers tested, neither adopted

Under `predict_method = "auto"` (PR #192) the route choice scores ordinal traits by squared error on the
latent z-scale. Two alternatives were tested against that merged code (bf7ca7a):

- A, class error rate: decode to the class (ordinal traits are coded 0 to K-1 before z-scoring), score the share
  of wrong classes, and break ties by squared error.
- B, log-probability of the true class: under N(mu, se^2) on the z-scale, with class boundaries halfway
  between the class codes; ties broken by squared error.

Benchmark: 18 core cells (mixed types, Brownian motion with Pagel's lambda, 30% MCAR), GNN off, 200 seeds,
Totoro, `~/R/lib-exact10`, with a temporary switch selecting the scorer (branch
`feat/ordinal-route-logloss`, not merged). Summaries: `core_lambda_ordA_off_agg_summary.csv`,
`core_lambda_ordB_off_agg_summary.csv`; merged code: `../exact-default/core_lambda_auto5_off_agg_summary.csv`.

Ordinal accuracy, change against the merged code (mean over cross-trait correlations 0 and 0.5):

| true lambda | species | merged | A | B |
|---|---|---|---|---|
| 0.3 | 100 | 0.324 | -0.005 | -0.013 |
| 0.3 | 300 | 0.328 | -0.004 | -0.020 |
| 0.3 | 1000 | 0.330 | -0.004 | -0.030 |
| 0.7 | 100 | 0.509 | 0.000 | -0.013 |
| 0.7 | 300 | 0.529 | -0.003 | -0.038 |
| 0.7 | 1000 | 0.552 | -0.001 | -0.052 |
| 1 | 100, 300, 1000 | 0.891 to 0.960 | 0.000 | 0.000 |

A lowers accuracy in 10 of 18 scenarios and raises it in 1; B lowers it in 12 and raises it in none. Macro F1
moves the same way. No other trait changed. Squared error stays.

An earlier run of A (branch `feat/ordinal-route-accuracy`) decoded classes as 1 to K and so scored every
lowest-class prediction as wrong; its result is superseded by the run above.

Separate finding, not addressed here: the ordinal label-propagation candidate in `fit_baseline()`
(`R/fit_baseline.R`, Phase F, commit a3b89e6) also decodes classes as 1 to K, which drops every lowest-class
observation from that candidate. It is chosen only when it wins on validation error, so the effect is most
likely that it rarely competes; this has not been measured.
