# Multiple imputation under phylogenetic GLS: results

2026-09-23. Branch `arc/mi-gls-attenuation`. Question: do pigauto's multiple-imputation
draws support valid inference in a phylogenetic GLS downstream analysis, and does a draw
from the proper Gaussian conditional distribution fix what they get wrong?

## Design

Two traits simulated as bivariate Brownian motion (correlation 0.7, so the true PGLS
slope of y on x is 0.70) on random trees, 16 regimes: lambda 1 or 0.5 (the tree is
Pagel-transformed for simulation only), n 300 or 1000, missingness MCAR or clade-structured
(MAR_phylo), 30% of cells missing in x only or in both traits. 120 replicates per regime,
m = 20 imputations, pooled by Rubin's rules (`pool_mi()`). Downstream: `gls(y ~ x,
corBrownian)` (PGLS) and `lm(y ~ x)` (OLS). Methods:

- `complete`: the analysis on the true values (the reference).
- `oracle`: exact conditional draws under the true model (lambda 1, x only).
- `draw_cond_em`: the prototype `draw_conditional_bm()`, which draws all missing cells
  jointly from their Gaussian conditional under a multivariate Brownian model with the
  trait covariance estimated by EM (`R/draws_conditional.R`).
- `draw_cond_inhouse`: the same draws with pigauto's current plug-in covariance.
- `single`: `impute()` point prediction; `mi_dropout`, `mi_conf_exact`,
  `mi_conf_per_column`: `multi_impute()` with MC-dropout and conformal draws (exact and
  default prediction routes). 500 epochs for all GNN fits.

GNN-based methods come from fir array 61137481 (1,920 cells, 9 reran with more time after
2 node faults and 7 timeouts); the draw methods, complete and oracle from fir array
61174611 on the identical per-cell datasets. Largest bias MCSE across cells: 0.016.
Scripts: `script/mi_gls/01_cell.R`, `03_combine.R`; per-cell data `cells.csv`; full
summary `summary.md` (bias, MCSE, SE ratio, coverage, coverage MCSE).

## Table 1. PGLS slope: bias / 95% CI coverage

Complete-data coverage is below 0.95 at lambda 0.5 because the analysis model assumes
lambda 1; compare each method with `complete`, not with 0.95.

| regime | lambda | n | mechanism | missing | complete | oracle | draw_cond_em | draw_cond_inhouse | single | mi_dropout | mi_conf_exact | mi_conf_per_column |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 300 | MCAR | x_only | -0.002 / 0.93 | -0.004 / 0.97 | -0.004 / 0.93 | -0.017 / 0.93 | -0.061 / 0.81 | -0.208 / 0.18 | -0.258 / 0.09 | -0.361 / 0.02 |
| 2 | 0.5 | 300 | MCAR | x_only | -0.006 / 0.78 | . | +0.042 / 0.68 | +0.021 / 0.81 | -0.038 / 0.77 | -0.074 / 0.67 | -0.237 / 0.17 | -0.311 / 0.07 |
| 3 | 1 | 1000 | MCAR | x_only | -0.002 / 0.91 | -0.000 / 0.95 | -0.001 / 0.90 | -0.010 / 0.92 | -0.057 / 0.51 | -0.225 / 0.05 | -0.242 / 0.00 | -0.350 / 0.00 |
| 4 | 0.5 | 1000 | MCAR | x_only | -0.002 / 0.76 | . | +0.056 / 0.39 | +0.034 / 0.66 | -0.008 / 0.77 | -0.030 / 0.66 | -0.196 / 0.03 | -0.274 / 0.00 |
| 5 | 1 | 300 | MAR_phylo | x_only | -0.001 / 0.93 | -0.001 / 0.94 | +0.008 / 0.91 | -0.009 / 0.93 | -0.041 / 0.80 | -0.223 / 0.23 | -0.265 / 0.10 | -0.359 / 0.03 |
| 6 | 0.5 | 300 | MAR_phylo | x_only | -0.008 / 0.78 | . | +0.049 / 0.57 | +0.024 / 0.80 | -0.021 / 0.82 | -0.061 / 0.74 | -0.250 / 0.17 | -0.323 / 0.05 |
| 7 | 1 | 1000 | MAR_phylo | x_only | +0.001 / 0.88 | +0.001 / 0.95 | +0.009 / 0.90 | -0.007 / 0.91 | -0.050 / 0.59 | -0.263 / 0.07 | -0.234 / 0.00 | -0.350 / 0.00 |
| 8 | 0.5 | 1000 | MAR_phylo | x_only | -0.009 / 0.74 | . | +0.052 / 0.43 | +0.029 / 0.70 | -0.010 / 0.78 | -0.037 / 0.62 | -0.209 / 0.05 | -0.289 / 0.01 |
| 9 | 1 | 300 | MCAR | both | -0.009 / 0.90 | . | -0.009 / 0.84 | -0.089 / 0.65 | -0.236 / 0.02 | -0.346 / 0.02 | -0.358 / 0.04 | -0.450 / 0.00 |
| 10 | 0.5 | 300 | MCAR | both | -0.005 / 0.83 | . | +0.046 / 0.72 | -0.112 / 0.43 | -0.254 / 0.04 | -0.278 / 0.04 | -0.380 / 0.07 | -0.440 / 0.03 |
| 11 | 1 | 1000 | MCAR | both | -0.001 / 0.93 | . | -0.002 / 0.89 | -0.082 / 0.21 | -0.235 / 0.00 | -0.366 / 0.00 | -0.346 / 0.00 | -0.441 / 0.00 |
| 12 | 0.5 | 1000 | MCAR | both | +0.002 / 0.76 | . | +0.053 / 0.47 | -0.101 / 0.28 | -0.244 / 0.00 | -0.258 / 0.00 | -0.352 / 0.01 | -0.418 / 0.00 |
| 13 | 1 | 300 | MAR_phylo | both | +0.001 / 0.97 | . | +0.008 / 0.89 | -0.086 / 0.65 | -0.226 / 0.03 | -0.360 / 0.03 | -0.376 / 0.03 | -0.459 / 0.02 |
| 14 | 0.5 | 300 | MAR_phylo | both | -0.007 / 0.68 | . | +0.035 / 0.64 | -0.105 / 0.48 | -0.245 / 0.10 | -0.275 / 0.06 | -0.391 / 0.02 | -0.448 / 0.00 |
| 15 | 1 | 1000 | MAR_phylo | both | +0.000 / 0.98 | . | +0.008 / 0.90 | -0.081 / 0.26 | -0.229 / 0.00 | -0.384 / 0.00 | -0.361 / 0.00 | -0.452 / 0.00 |
| 16 | 0.5 | 1000 | MAR_phylo | both | -0.006 / 0.78 | . | +0.039 / 0.49 | -0.089 / 0.33 | -0.228 / 0.01 | -0.247 / 0.01 | -0.355 / 0.03 | -0.418 / 0.02 |

## Table 2. OLS slope: bias

OLS ignores phylogenetic dependence, so its confidence intervals are invalid even on
complete data (coverage 0.28 to 0.68); only bias is informative here.

| regime | lambda | n | mechanism | missing | complete | oracle | draw_cond_em | draw_cond_inhouse | single | mi_dropout | mi_conf_exact | mi_conf_per_column |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 300 | MCAR | x_only | +0.004 | +0.003 | +0.002 | -0.001 | +0.006 | -0.028 | -0.043 | -0.081 |
| 2 | 0.5 | 300 | MCAR | x_only | -0.004 | . | -0.053 | -0.070 | -0.021 | -0.043 | -0.169 | -0.214 |
| 3 | 1 | 1000 | MCAR | x_only | +0.013 | +0.013 | +0.013 | +0.011 | +0.017 | -0.014 | -0.013 | -0.040 |
| 4 | 0.5 | 1000 | MCAR | x_only | -0.001 | . | -0.050 | -0.067 | -0.011 | -0.024 | -0.140 | -0.184 |
| 5 | 1 | 300 | MAR_phylo | x_only | -0.002 | +0.000 | +0.002 | -0.003 | +0.004 | -0.037 | -0.042 | -0.084 |
| 6 | 0.5 | 300 | MAR_phylo | x_only | -0.003 | . | -0.062 | -0.083 | -0.017 | -0.039 | -0.175 | -0.225 |
| 7 | 1 | 1000 | MAR_phylo | x_only | -0.006 | -0.005 | -0.004 | -0.007 | -0.002 | -0.040 | -0.021 | -0.054 |
| 8 | 0.5 | 1000 | MAR_phylo | x_only | +0.006 | . | -0.058 | -0.077 | -0.003 | -0.018 | -0.132 | -0.180 |
| 9 | 1 | 300 | MCAR | both | +0.000 | . | +0.002 | -0.019 | -0.043 | -0.077 | -0.074 | -0.115 |
| 10 | 0.5 | 300 | MCAR | both | -0.004 | . | -0.035 | -0.126 | -0.163 | -0.179 | -0.268 | -0.316 |
| 11 | 1 | 1000 | MCAR | both | +0.006 | . | +0.005 | -0.011 | -0.030 | -0.062 | -0.050 | -0.082 |
| 12 | 0.5 | 1000 | MCAR | both | +0.010 | . | -0.024 | -0.112 | -0.138 | -0.148 | -0.223 | -0.278 |
| 13 | 1 | 300 | MAR_phylo | both | -0.001 | . | +0.000 | -0.033 | -0.066 | -0.106 | -0.095 | -0.146 |
| 14 | 0.5 | 300 | MAR_phylo | both | -0.001 | . | -0.041 | -0.133 | -0.173 | -0.191 | -0.280 | -0.334 |
| 15 | 1 | 1000 | MAR_phylo | both | +0.002 | . | +0.004 | -0.019 | -0.044 | -0.080 | -0.058 | -0.099 |
| 16 | 0.5 | 1000 | MAR_phylo | both | +0.001 | . | -0.044 | -0.127 | -0.149 | -0.160 | -0.229 | -0.285 |

## Findings

1. **pigauto's current MI draws are not valid for PGLS.** In all 16 regimes the conformal
   draws, with either prediction route, bias the PGLS slope by -0.20 to -0.46 and the
   pooled 95% intervals cover the truth in 0 to 17% of replicates. MC-dropout draws are
   biased by -0.03 to -0.38. Single imputation is nearly unbiased when only x is missing
   but biased by about -0.23 when both traits are. The problem is not specific to the
   Mondrian scale or to clade-structured missingness.
2. **Proper conditional draws fix it when the Brownian model is right.** At lambda 1,
   `draw_cond_em` is unbiased in every regime (|bias| at most 0.009) and its coverage
   (0.84 to 0.93) is close to complete data's (0.88 to 0.98). The small shortfall is
   expected: the draws condition on one estimated covariance and do not propagate its
   uncertainty.
3. **The covariance estimate matters when several traits are missing.** With pigauto's
   plug-in covariance the same draws are biased by -0.08 to -0.11 in the both-missing
   regimes; with EM, by -0.009 to +0.053. The plug-in estimate shrinks cross-trait
   correlation under missingness (0.67 to 0.46 with 30% of both traits missing, six
   trees, n = 300; `max_iter = 50` gives the same value). That estimator is pigauto's
   joint-MVN baseline solver (`R/joint_mvn_solver.R`), so the same shrinkage reaches the
   baseline's own predictions.
4. **Misspecified phylogenetic signal remains.** At lambda 0.5, `draw_cond_em` assumes
   lambda 1 and is biased by +0.035 to +0.056 under PGLS with coverage 0.39 to 0.72
   (complete: 0.68 to 0.83). A lambda-aware conditional draw is the natural next step,
   and it connects directly to the lambda-estimating baseline being built on
   `feat/joint-lambda-default`.

## Recommendation

- Offer proper conditional draws as a `multi_impute()` option for continuous traits
  (`draws_method = "conditional"`), with the EM covariance and, once the lambda lane
  lands, the estimated lambda. Keep the default unchanged until Shinichi decides; the
  NEWS caveat on conformal MI under PGLS stays.
- Fix or document the plug-in covariance in the joint-MVN solver: its cross-trait
  shrinkage under missingness is a baseline issue, owned by the lambda lane.

## Scope

Two continuous traits, bivariate Brownian data with lambda 1 or 0.5, 30% missing, n 300
and 1000, one analysis model family (linear slope under PGLS and OLS). Not covered:
discrete traits, more than two traits, multi-observation data, non-Brownian processes,
analysis models other than a single slope, and trees with real rather than simulated
missingness patterns. `draw_conditional_bm()` is an internal prototype, not an exported
function.
