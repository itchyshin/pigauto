# Posterior multiple imputation: results

2026-09-24/25. Branch `arc/mi-posterior`, draft PR itchyshin/pigauto#189 (base `main`; #187 merged).

**Status:**
- Simulation complete.
- In-model G7 passes. In-model G6 fails on 3 rows, all of them the relative SE-ratio band, and the
  count matches what Monte Carlo noise predicts. Every other gated rule passes.
- Real data: 9 of 10 cells done; the FishBase cell is still running.
- Not ready to merge.

Every number below comes from a committed file:
- `sim_summary.csv` and `cell_coverage.csv` (from `script/mi_gls/03_summarise_v2.R`);
- `results_tables.md` (`06_results_tables.R`);
- `se_ratio_noise.txt` (`07_se_ratio_noise.R`);
- `evidence/` for G3, the smoke runs and the byte-identity checks;
- `diagnosis.md` for the round-1 failure analysis.

## Provenance

The simulation has 40 regimes x 200 reps = 8,000 cells, merged from two runs.
- **Commit 69670d4** (first campaign): 4,756 converged cells of regimes 1 to 24.
- **Commit 9597e18** (round 2), which adds automatic chain extension:
  - the 3,200 cells of the in-model twin regimes 25 to 40;
  - the 44 cells of regimes 1 to 24 that had not converged;
  - 14 converged cells re-run by accident (see "What did not go smoothly" in the after-task report).

The code change between the two commits affects a cell only when its chains would be extended. On
converged cells the new code is byte-identical to 69670d4:
- 16 of 16 converged cells re-run at 9597e18 match their 69670d4 files (results, cell detail, cell
  coverage, diagnostics; `evidence/smoke_round2/`);
- a pinned regression test in `test-mi-posterior.R` checks the same thing.

`03_summarise_v2.R` reports `MIXED_CODE_SHA` for the merged set and `SETTINGS 0 of 8000` non-campaign
files.

## What was built

`multi_impute(traits, tree, draws_method = "posterior")` for continuous traits. The missing cells are
drawn jointly from vec(Y) ~ N(1 mu', Sigma_P %x% R + Sigma_E %x% I_n), with full Sigma_P and
Sigma_E, R = cov2cor(vcv(tree)) and per-trait lambda_k = Sigma_P[k,k] / (Sigma_P[k,k] + Sigma_E[k,k]).
- **Proper imputation.** Sigma_P, Sigma_E and mu are redrawn by Gibbs data augmentation, with
  parameter expansion and collapsed Metropolis moves, so the m completed datasets carry parameter
  uncertainty.
- **Per-cell 95% intervals** come from 1,000 kept posterior-predictive draws.
- **Convergence.** Split R-hat and bulk ESS are computed internally. A fit that fails R-hat < 1.05
  and ESS > 400 is extended automatically, up to 4x its default length.
- **Downstream analysis.** The results pass through `with_imputations()` and `pool_mi()` under the
  provenance marker `pigauto_posterior_mi_v1`.
- The GNN is not used, and the default draws method is unchanged.

Design: `design.md`. Reviews: `review-design.md`, `review.md` (S5), and the round-2 review and the
claims audit (both recorded in the after-task report).

## Simulation design

- Two traits, n = 300 or 1000, 30% missing (x only, or both traits), MCAR or clade-biased
  (MAR_phylo) masks.
- Downstream analyses: `nlme::gls(corBrownian)` and `phylolm(model = "lambda")`, each against the
  complete data of the same replicate.
- **Gated, in the sampler's model family** (Shinichi, CP2 follow-up):
  - regimes 17 to 24: two lambdas, phylogenetic and residual correlations that differ;
  - regimes 25 to 40: in-model twins of regimes 1 to 16. They use the same trees, seeds and masks,
    with each tip's row rescaled so that the data come from Sig %x% (lambda R + (1 - lambda) I).
- **Stress test, reported only:** regimes 1 to 16, which simulate from the raw covariance of
  non-ultrametric trees. That is outside the model every pigauto imputation path fits
  (`diagnosis.md`).

## Headline 1: per-cell predictive coverage is about 95%

| Regimes | Mask | Posterior 95% interval coverage | Conformal (`impute(gnn = FALSE)`) coverage | Posterior width / conformal width |
|---|---|---|---|---|
| In-model 17 to 40 (20 regime x trait rows) | MCAR (gated) | 0.933 to 0.955 | 0.954 to 0.971 | 0.54 to 0.92 |
| In-model 17 to 40 (20 rows) | clade-biased (reported) | 0.941 to 0.952 | 0.944 to 0.966 | 0.57 to 0.93 |
| Stress 1 to 16 (12 rows) | MCAR | 0.943 to 0.957 | 0.954 to 0.970 | 0.57 to 0.76 |
| Stress 1 to 16 (12 rows) | clade-biased | 0.932 to 0.941 | 0.936 to 0.958 | 0.60 to 0.79 |

Every gated row lies inside the G7 band [0.92, 0.98]. The posterior intervals are 8% to 46% narrower
than conformal and still cover about 95%.

## Headline 2: the downstream slope

These are the in-model regimes 17 to 40, 48 regime x analysis rows (`results_tables.md`).

- **Paired bias** (MI pooled slope minus complete-data slope, same replicate):
  - twins: gls -0.014 to +0.009, phylolm -0.011 to +0.003;
  - two-lambda regimes: gls -0.008 to +0.006, phylolm -0.014 to +0.001;
  - all 48 rows are within the gate of max(0.02, 2.5 MCSE);
  - for comparison, conformal MI draws gave -0.20 to -0.46 on regimes 1 to 16 of the earlier sweep
    (same data-generating seeds; `arc/mi-gls-attenuation`, `docs/dev-log/mi-gls/results.md`).
- **Coverage of the pooled 95% interval:**
  - 0.845 to 0.975, against 0.76 to 0.97 for complete data under the same analysis model;
  - within 0.05 of complete data in 48 of 48 rows;
  - mean positive-part shortfall 0.0016.
- **SE ratio under phylolm, relative to complete data:**
  - 0.96 to 1.19, with mean 1.059 over the 24 gated rows;
  - so the pooled SEs are honest to slightly conservative (see the next section).
- **Proper vs plug-in** (reported, not gated). Holding Sigma_P and Sigma_E at their posterior means
  gives a lower SE ratio in 48 of 48 both-missing pairs (all 40 regimes). Coverage drops to
  0.77 to 0.945, against 0.88 to 0.975 for the proper draws. So propagating parameter uncertainty
  matters downstream, even though per-cell intervals barely change (widths differ by under 0.2%).
- **Convergence:**
  - 7,999 of 8,000 fits converged (4,799 of 4,800 in the gated regimes);
  - the 44 fits that failed at 69670d4 were all re-run at 9597e18, extended automatically, and all
    converged.

## Gate verdicts (`.unlazy/mi-posterior/GATES.md`)

| Gate | Result | Detail |
|---|---|---|
| G1 unit tests | met | two posterior test files: 118 + 196 expectations, 0 failures (includes the pinned 69670d4 regression test) |
| G2 exactness | met | `EXACTNESS_OK` |
| G3 calibration | met | `RECOVERY_OK`: 200 fits at 69670d4, 95% interval coverage of lambda and rho_P 0.92 to 0.98 |
| G4 smoke convergence | met | `CONVERGENCE_OK` |
| G5a suite / G5b solver / G5c check | met at 9059f87; re-verified at the final HEAD before merge | see the after-task report |
| **G6 simulation acceptance** | **not met (3 in-model rows, one rule)** | relative SE ratio above 1.15 under phylolm in twins 35 (1.186), 38 (1.167), 36 (1.153); every other gated rule passes |
| G7 per-cell coverage | **met** | `CELL_COVERAGE_PASS` |
| G8 real data | pending | FishBase cell still running |

### Why the three G6 failures are noise, and what to decide (`se_ratio_noise.txt`)

The relative rule divides the MI SE ratio by the complete-data SE ratio. Each is an empirical SD over
200 reps, so the quotient has a Monte Carlo SE of about 0.075 per row.

- Across the 24 gated phylolm rows the relative ratio averages 1.059.
- Its SD across regimes (0.058) is below its per-row MCSE (0.075), so the between-regime spread is no
  larger than noise.
- Noise alone predicts 3.1 of 24 rows outside [0.90, 1.15]. 3 were observed.
- The pooled mean is 1.059 with SE 0.015. So on average the MI SEs are about 6% conservative relative
  to complete data. That is a real but small effect, on the safe side.
- The plan's original absolute rule would have failed 5 rows in the other direction (MI ratio as low
  as 0.80). There, the complete-data analysis is itself miscalibrated (complete ratio 0.79 to 1.01).

**Proposal for Shinichi (not applied): one of**
- (a) gate the relative ratio with an MCSE-aware band (for example 1 +/- 2.5 MCSE per row);
- (b) gate the pooled mean relative ratio across the 24 rows (1.059, band [0.95, 1.10]);
- (c) accept G6 as failed on this rule and report it.

The gate is unchanged pending your decision.

## Stress test: regimes 1 to 16 (raw covariance of non-ultrametric trees; reported only)

G6 rule outcomes (`04_acceptance.R`, STRESS block): 5 violations in 5 of 16 regimes.
- Paired bias under phylolm in regime 1 (-0.026) and regime 9 (-0.020).
- Coverage shortfall in regime 3 (0.050, at the limit).
- Relative SE ratio in regime 8 (0.879) and regime 15 (0.897).

G7 per-cell coverage passes in every stress regime (0.943 to 0.957 MCAR). `diagnosis.md` explains the
bias:
- the dominant part is the model mismatch (tip-variance heterogeneity that `cov2cor` removes);
- the in-model twins remove most of it (twin phylolm bias -0.011 to +0.003, against -0.026 to
  +0.001 here).

For ultrametric (dated) trees `vcv(tree)` is proportional to `cov2cor(vcv(tree))`, so this mismatch
does not arise.

## How we got here (round 1)

The first campaign (24 regimes, commit 69670d4) failed G6 on 8 checks: non-convergence in three
n = 300 regimes, the lambda = 1 bias, one coverage shortfall and two SE ratios. `diagnosis.md`:
- **The lambda = 1 bias** was mostly the harness DGP lying outside the model family. Exact
  imputations under the true covariance were unbiased, and the in-model twin removed 86% of the bias.
- **The non-convergence** was an ESS shortfall only. Longer chains converged all 32 failures without
  moving the slopes.

Shinichi chose to add in-model twins and automatic chain extension (design.md 5e). The CP1 decisions
(design.md 5d) were also made before any result: the G3 calibration gate, the relative SE rule, and
reporting proper vs plug-in without gating it.

## Real data (preview: 9 of 10 cells; G8 pending the FishBase cell)

All real data ran at 69670d4, because G8 requires one code SHA.
- **Data:** PanTHERIA (4,027 species, 4 continuous traits; 3 MCAR and 3 clade-structured masks) and
  AVONET (1,500 species, 4 traits; 3 MCAR masks).
- **Convergence:** 7 of 9 cells converged. Two PanTHERIA cells reached min ESS 368 and 366 with R-hat
  at most 1.01; automatic extension would have extended them.
- **Per-trait coverage:** model-based 0.87 to 0.96 across 36 trait-cells, against split conformal
  0.90 to 0.98 and Mondrian 0.91 to 0.99. On real data the model-based intervals cover somewhat less
  than conformal.
- **Slopes:**
  - AVONET: pooled slopes within about 3% of the complete-row reference.
  - PanTHERIA: some cells differ by up to 3 reference SEs (for example longevity ~ body mass: 0.148
    against 0.180).
  - Likely reason: the PanTHERIA columns are already log values, and the default
    `log_transform = TRUE` logs them again. The imputation model is then linear on a log(log) scale
    while the pre-registered analysis is on the log scale, which is outside the stated congeniality
    scope.
  - A `log_transform = FALSE` sensitivity run is the natural follow-up.
- The full G8 report follows when FishBase finishes.

## What this does not cover

- Discrete, mixed and multi-observation traits.
- GNN blending.
- Covariates in the imputation model.
- Analyses with nonlinear terms, interactions or external covariates.
- Non-ultrametric trees: supported, but the imputation model treats tips as equal-variance, as all of
  pigauto does (stress test above).
- More than two traits in simulation (real data has 4 or 5).
- Trees much beyond 10,000 tips. FishBase (10,484 tips, 5 traits) takes about 1 s per sweep, so a
  default fit takes hours.

The default draws method is unchanged (Shinichi's decision).
