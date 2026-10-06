# Discrete-trait lambda prototype: screen against BACE (2026-10-05)

Requested by Shinichi ("prototype λ for discrete traits"; aim: catch up with BACE on categorical accuracy).
Branch `feat/discrete-lambda`. Prototype switch, default off: `options(pigauto.discrete_lambda = "estimate")` adds
the binary and ordinal liability columns to the joint solver's `lambda_cols` in `fit_joint_threshold_baseline()`;
categorical traits get it through the one-vs-rest fits. Before this, every discrete trait was imputed at lambda = 1
(the "section 7 B iii cut" of D-278).

## Design

- Datasets: the four-arm study's core slice (branch `arc/imputation-sim`, `script/campaign_sim_cell.R`,
  `types_mixed`, BM, MCAR 30%, fixed thresholds, always-observed driver), n = 100 and 300, lambda = 0.3, 0.7, 1,
  rho = 0 and 0.5, seeds 1 to 100: the datasets BACE ran on. BACE's stored per-seed results
  (Totoro `~/pigauto_sim/results/core_bace/`) are compared dataset by dataset; 1,197 of 1,200 paired (BACE files
  missing for 3). The re-run mode floor matched the stored floor exactly on the smoke dataset.
- pigauto: the prototype (30eeedb) in a private Totoro library, `impute(..., gnn = FALSE)` as in the study, with the
  safety machinery on (`gnn_off`, the default) and off (`gnn_off_pure`: `safety_floor = FALSE,
  phylo_signal_gate = FALSE`), each with the switch off and on. Totoro, 80 cores, 2 x 1,200 jobs, about 4 min each.
- Screen 2 also tried giving the ordinal Brownian-motion candidate an estimated lambda (a4718d8); it lowered ordinal
  accuracy by 0.009 to 0.018 and was reverted (8540543). The tables below use screen 2's other columns, which are
  unaffected; `dlam_screen1.txt` has the ordinal figures without that change.

## Result (`dlam_screen2.txt`)

Discrete accuracy pooled over binary, ordinal and three-class traits and rho:

| n | lambda | mode floor | current default (gate on, lambda 1) | gate off, lambda estimated | BACE |
|---|---|---|---|---|---|
| 100 | 0.3 | 0.459 | 0.462 | 0.526 | 0.526 |
| 100 | 0.7 | 0.545 | 0.617 | 0.664 | 0.640 |
| 100 | 1 | 0.626 | 0.900 | 0.923 | 0.768 |
| 300 | 0.3 | 0.469 | 0.476 | 0.546 | 0.566 |
| 300 | 0.7 | 0.560 | 0.647 | 0.685 | 0.687 |
| 300 | 1 | 0.641 | 0.952 | 0.953 | 0.792 |

(screen 1 figures for the gate-off, lambda-estimated column, i.e. without the reverted ordinal change.)

- Paired against BACE, the gap at lambda = 0.3 shrinks from -0.064 / -0.090 (n = 100 / 300) to 0.000 / -0.020.
- Binary and categorical: at or above BACE in every cell (n = 300, lambda = 0.3: 0.712 vs 0.710 and 0.548 vs 0.544).
- Ordinal is the remaining gap at low signal: 0.356 / 0.379 vs BACE 0.410 / 0.443 at lambda = 0.3.
- Both parts are needed. Lambda alone (gate on) moves little, because the phylo-signal gate replaces weak-signal
  discrete predictions with the mode; gate off alone (lambda 1) leaves categorical below the mode floor at
  n = 300, lambda = 0.3 (0.448 vs 0.450).
- Discrete Brier score: gate off, lambda estimated is the lowest of all arms in 5 of 6 cells, and below BACE in all 6.

Turning the gate off also changes continuous traits (zRMSE pooled over c1, c2, cnt, prp; lower is better):

| n | lambda | gate on | gate off | BACE |
|---|---|---|---|---|
| 100 | 0.3 | 0.971 | 0.892 | 1.013 |
| 300 | 0.3 | 0.930 | 0.857 | 0.921 |
| 300 | 1 | 0.389 | 0.386 | 0.738 |

Gate off is better in every cell, and conformal coverage is unchanged or slightly higher (n = 300: 0.963 to 0.966;
n = 100: 0.880 to 0.893). The September GNN-off campaign found the opposite on its low-signal `bace_dgp` (the
safety machinery helped), when the continuous baseline still ran at lambda = 1; since D-278 lambda is estimated, which
may make the gate redundant. That has not been tested on `bace_dgp`, OU or AVONET.

## Does NOT cover

One DGP family (the four-arm core: BM, MCAR, fixed thresholds); n = 1000; multi-observation data; the LP path
(categorical traits with no continuous trait alongside); real data; the gate change on other DGPs; posterior MI
(continuous-only).
