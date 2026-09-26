# Ordinal route choice by class error rate: tested, not shipped

Branch `feat/ordinal-route-accuracy` (commit b182e2a) changed the "auto" route choice to score ordinal traits
by class error rate, with latent squared error breaking ties. Main at bf7ca7a (PR #192) scores them by latent
squared error.

Benchmark: 18 core cells, GNN off, 200 seeds, Totoro, `~/R/lib-exact9`. Summary:
`core_lambda_ord1_off_agg_summary.csv`, compared with `../exact-default/core_lambda_auto5_off_agg_summary.csv`
(the merged code).

Ordinal accuracy, change against the merged code (mean over two cross-trait correlations):

| true lambda | n 100 | n 300 | n 1000 |
|---|---|---|---|
| 0.3 | -0.005 | -0.006 | -0.006 |
| 0.7 | -0.001 | -0.010 | -0.007 |
| 1 | 0.000 | 0.000 | 0.000 |

Accuracy is lower in 13 of 18 scenarios and higher in one (+0.001); the largest single-scenario fall is -0.014
(unpaired |z| up to 1.4). Macro F1 moves the same way. No other trait changed. The merged code's ordinal
accuracy is already 0.011 to 0.044 above the previous default at lambda 0.3 and 0.7 and within 0.001 at
lambda 1.

Decision: not shipped. A likely reason (not tested) is that a 0-1 loss on half of a trait's validation rows is
a noisier route selector than squared error, which keeps how far each prediction is from the truth. The NEWS
"Known limitation" line about ordinal scoring on main is now out of date.
