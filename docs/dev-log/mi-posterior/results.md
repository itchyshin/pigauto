# Posterior multiple imputation: results

2026-09-24/25. Branch `arc/mi-posterior`, draft PR itchyshin/pigauto#189 (base `main`; #187 merged).

**Status:**
- Simulation complete.
- In-model G7 passes. In-model G6 fails on 3 rows, all of them the relative SE-ratio band (twins 35,
  36 and 38 under phylolm). Monte Carlo noise does not explain them: noise alone predicts 0.6 rows
  outside the band (P(3 or more) = 0.02). In all three the MI SEs are well calibrated in absolute
  terms (MI SE ratio 0.98 to 1.06); they fail because the complete-data analysis is over-confident
  there (ratio 0.85 to 0.90). Every other gated rule passes. A decision is proposed below; no gate
  has been changed.
- Real data: 9 of 10 cells done; the FishBase cell is still running.
- Not ready to merge.

Every number below comes from a committed file, except the G1, G2 and G5 run outputs named in the
gate table:
- `sim_summary.csv` and `cell_coverage.csv` (from `script/mi_gls/03_summarise_v2.R`; the merge is
  logged in `evidence/summarise_round2.log`);
- `results_tables.md` (`06_results_tables.R`);
- `se_ratio_noise.txt` (`07_se_ratio_noise.R`);
- `04_acceptance.R` run on `sim_summary.csv` for the G6 rule outcomes;
- `evidence/` for G3, G4, the smoke runs and the byte-identity checks;
- `diagnosis.md` for the round-1 failure analysis;
- `real_preview/` for the real-data preview;
- `docs/dev-log/mi-gls/results.md` on branch `arc/mi-gls-attenuation` for the earlier conformal
  sweep.

## Provenance

The simulation has 40 regimes x 200 reps = 8,000 cells, merged from two runs: 4,756 + 44 + 3,200.
- **Commit 69670d4** (first campaign): the 4,756 cells of regimes 1 to 24 that converged.
- **Commit 9597e18** (round 2), which adds automatic chain extension:
  - the 44 cells of regimes 1 to 24 that had not converged at 69670d4;
  - the 3,200 cells of the in-model twin regimes 25 to 40.
- 14 of the 4,756 converged cells (regime 1 reps 1 to 10, regime 11 reps 1 to 4) were also re-run
  at 9597e18 by accident. The merge keeps the later copy, so 4,742 rep files come from 69670d4 and
  3,258 from 9597e18. The 14 re-runs are byte-identical to their 69670d4 files on every compared
  field. `evidence/smoke_round2/README.md` says how the accident happened.

`03_summarise_v2.R` logs the merge: 4,800 files from the 69670d4 folder; 3,258 from the 9597e18
folder, replacing 58 (the 44 plus the 14); `SETTINGS 0 of 8000` non-campaign files; and
`MIXED_CODE_SHA`. A re-run of the merge reproduces the committed `sim_summary.csv`, `sim_summary.md`
and `cell_coverage.csv` byte for byte (`evidence/summarise_round2.log`).

The code change between the two commits affects a cell only when its chains would be extended. On
converged cells the new code is byte-identical to 69670d4:
- 16 comparisons of 15 distinct converged cells re-run at 9597e18 match their 69670d4 files
  (results, cell detail, cell coverage, core diagnostics, missing-cell counts). Regime 1 rep 1 was
  compared twice, once in the smoke run and once among the 14 (`evidence/smoke_round2/`).
- A pinned regression test in `test-mi-posterior.R` checks the same thing.

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

Design: `design.md`. Reviews: `review-design.md` and `review.md` (S5). The round-2 review and the
claims audit are recorded only through their fixes (commits 9597e18, 79bf98e and 8d2f612); their
findings are not in a committed file. The M2 statistical review corrected the SE-ratio noise analysis
and several numbers in this document.

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

## Headline 1: per-cell predictive coverage is close to 95%

| Regimes | Mask | Posterior 95% interval coverage | Conformal (`impute(gnn = FALSE)`) coverage | Posterior width / conformal width |
|---|---|---|---|---|
| In-model 17 to 40 (20 regime x trait rows) | MCAR (gated) | 0.933 to 0.955 | 0.954 to 0.971 | 0.54 to 0.92 |
| In-model 17 to 40 (20 rows) | clade-biased (reported) | 0.941 to 0.952 | 0.944 to 0.966 | 0.57 to 0.93 |
| Stress 1 to 16 (12 rows) | MCAR | 0.943 to 0.957 | 0.954 to 0.970 | 0.57 to 0.76 |
| Stress 1 to 16 (12 rows) | clade-biased | 0.932 to 0.941 | 0.936 to 0.958 | 0.60 to 0.79 |

Every gated row lies inside the G7 band [0.92, 0.98]. The posterior intervals are 7% to 46% narrower
than conformal across every row of the table (8% to 46% in the gated in-model MCAR rows).

Coverage is close to 95% but not uniformly nominal (`cell_coverage.csv`):
- In the two-lambda regimes 17 to 24 all 16 regime x trait rows are below 0.95: 0.933 to 0.949
  (mean 0.944), 0.1 to 1.7 points short. Each row scores 17,811 to 60,235 cells, so the binomial SE
  is 0.001 to 0.002, and 15 of the 16 rows are more than 2 SE below 0.95. That SE treats cells as
  independent. Cells of one replicate are correlated, so it understates the MCSE somewhat, but the
  largest shortfall (row 17 y, 1.7 points, 9 binomial SE) is not noise.
- All 12 lambda = 0.5 twin rows are below 0.95 (0.942 to 0.949); in regimes 28, 34, 36, 38 and 40
  every row is more than 2 binomial SE below (26, 30 and 32 are within 2 SE).
- The lambda = 1 MCAR twins cover 0.950 to 0.955.
- Conformal over-covers (0.954 to 0.971 in the gated MCAR rows). Part of the width gap therefore
  comes from the two coverage offsets. Rescaled to equal coverage under a normal approximation, the
  posterior intervals are still 4% to 41% narrower in the gated MCAR rows.

## Headline 2: the downstream slope

These are the in-model regimes 17 to 40, 48 regime x analysis rows (`results_tables.md`), unless a
bullet says otherwise.

- **Paired bias** (MI pooled slope minus complete-data slope, same replicate):
  - twins: gls -0.014 to +0.009, phylolm -0.011 to +0.003;
  - two-lambda regimes: gls -0.008 to +0.006, phylolm -0.014 to +0.001;
  - all 48 rows are within the gate of max(0.02, 2.5 MCSE).
- **The in-model bias is small but not zero at lambda = 1.**
  - All 16 twin rows at lambda = 1 (regimes 25, 27, ..., 39; gls and phylolm) are negative: -0.004 to
    -0.014, that is 0.5% to 2% of the true slope of 0.7 and 2.7 to 8.7 MCSE from zero (gls -0.005
    to -0.014 at 3.8 to 8.7 MCSE; phylolm -0.004 to -0.011 at 2.7 to 6.3 MCSE).
  - None of these 16 rows would pass on 2.5 MCSE alone; they pass through the rule's 0.02 floor.
  - The lambda = 0.5 twins show no such negative bias: phylolm -0.0004 to +0.0025 (every |z| < 2).
    For gls the both-missing lambda = 0.5 twins 34 to 40 lean positive (+0.008 to +0.009, 2.1 to 3.6
    MCSE).
  - Attribution (partial, and an inference): on the stress regimes, `diagnosis.md` measured a
    sampler component of similar size (campaign minus the KL oracle, -0.002 to -0.009 in regimes 1,
    3, 5, 11 and 13). It attributed roughly half to all of it, in regimes 1 and 3, to the Sigma_E prior,
    which pulls the residual correlation towards 0 when Sigma_E is small (lambda near 1). The prior
    arm ran only on regimes 1 and 3, not on the twins. `diagnosis.md` also measured a finite-n
    component of about -0.003 at n = 300 with only x missing.
- **Comparison with the earlier conformal draws** (like for like). The stress regimes 1 to 16 use
  the earlier sweep's DGP and seed formula (`dgp_v2.R` header). There the posterior paired bias is
  -0.013 to +0.017 under gls, and -0.026 to +0.017 across gls and phylolm (`sim_summary.csv`). On the
  same regimes the earlier sweep's conformal MI draws gave a gls (PGLS) bias of -0.20 to -0.46
  against the true slope of 0.7 (120 reps; `mi_conf_exact` and `mi_conf_per_column`,
  `docs/dev-log/mi-gls/results.md` Table 1 on `arc/mi-gls-attenuation`). The estimands differ
  (paired bias against bias versus the truth), but the complete-data bias in that sweep was at most
  0.009 in absolute value, which is negligible next to this gap.
- **Coverage of the pooled 95% interval:**
  - 0.845 to 0.975, against 0.76 to 0.97 for complete data under the same analysis model;
  - within 0.05 of complete data in 48 of 48 rows;
  - mean positive-part shortfall 0.0016.
- **SE ratio under phylolm, relative to complete data:**
  - 0.96 to 1.19, with mean 1.059 over the 24 gated rows;
  - so the pooled SEs are honest to slightly conservative relative to the complete-data analysis.
    In absolute terms the MI SE ratio is 0.80 to 1.06, and below 0.90 in 5 rows where the
    complete-data ratio is just as low (see the G6 section).
- **Proper vs plug-in** (reported, not gated). Plug-in rows exist only in the 24 both-missing
  regimes (9 to 24 and 33 to 40), so this bullet also counts 16 stress-test pairs. Holding Sigma_P
  and Sigma_E at their posterior means gives a lower SE ratio in 48 of 48 pairs (32 of 32 in-model
  pairs, 16 of 16 stress-test pairs). Coverage drops to 0.77 to 0.945, against 0.88 to 0.975 for the
  proper draws; the ranges are the same over the in-model pairs alone. So propagating parameter
  uncertainty matters downstream, even though per-cell intervals barely change (plug-in widths
  0.002% to 0.18% narrower).
- **Convergence:**
  - 7,999 of 8,000 fits converged (4,799 of 4,800 in the gated regimes; the exception is in
    regime 33, 199 of 200);
  - the 44 fits that failed at 69670d4 were all re-run at 9597e18, extended automatically, and all
    converged.

## Gate verdicts

The gate ledger is `.unlazy/mi-posterior/GATES.md`, which is local and not committed. Each row names
its committed source where one exists.

| Gate | Result | Detail |
|---|---|---|
| G1 unit tests | met | two posterior test files: 118 + 196 expectations, 0 failures (includes the pinned 69670d4 regression test); `evidence/ledger/g1_counts.log` and `evidence/ledger/gate_reverify_2026-09-25.log` |
| G2 exactness | met | `EXACTNESS_OK` (`evidence/ledger/gate_reverify_2026-09-25.log`) |
| G3 calibration | met | `RECOVERY_OK` (`gate_calibration.R check` on `evidence/g3_calibration.rds`): 200 fits at 69670d4, 50 per setting. For the gated quantities (lambda in settings B to D; rho_P in A to C, where both lambda >= 0.3), 95% interval coverage is 0.92 to 0.98. Not gated: lambda = 1 in setting A sits on the boundary, so no interval can cover it (coverage 0.00; bias -0.004 and -0.003); rho_P in setting D (lambda_1 = 0.05) is weakly identified (bias -0.14, coverage 0.98). See `evidence/README.md` |
| G4 smoke convergence | met | `CONVERGENCE_OK` (`evidence/smoke/g4.log`) |
| G5a suite / G5b solver / G5c check | met at the final code | `gate-check --reverify` re-ran every runnable gate: `SUITE_GREEN`, `SOLVER_UNTOUCHED`, R CMD check --as-cran `ERR=0 WARN=0` (`evidence/ledger/gate_reverify_2026-09-25.log`; its `R/` and `tests/` trees, e99198b and 9f5ac75, equal the branch head's) |
| **G6 simulation acceptance** | **not met (3 in-model rows, one rule)** | relative SE ratio above 1.15 under phylolm in twins 35 (1.186), 38 (1.167) and 36 (1.153); more than Monte Carlo noise predicts; every other gated rule passes |
| G7 per-cell coverage | **met** | `CELL_COVERAGE_PASS` |
| G8 real data | pending | FishBase cell still running |
| M2 fresh review | pending final confirmation | the M2 statistical review's findings are fixed (f03837b, 1c8e4da) and a three-member completion panel (D-43) withheld no claim; the formal PROCEED is recorded after the real-data section is final |

### The three G6 failures are not Monte Carlo noise (`se_ratio_noise.txt`)

The relative rule divides the MI SE ratio by the complete-data SE ratio. Each ratio is a mean SE over
an empirical SD from the same 200 reps, and within a rep the MI and complete-data slopes are strongly
correlated (rho 0.64 to 0.93, recovered from `sim_summary.csv`). The Monte Carlo SE (MCSE) of the
relative ratio must allow for that pairing: rel x sqrt((1 - rho^2) / (R - 1)), which is 0.027 to
0.057 per row (mean 0.044). The earlier version of this section treated the two SDs as independent
(0.075 per row). That overstated the noise 1.8-fold on average and led to the wrong conclusion that
the failures were noise. On the committed first-campaign per-rep slopes (regimes 17 to 24), a paired
bootstrap agrees with the corrected formula (0.049 against 0.050; independence 0.071; header of
`07_se_ratio_noise.R`).

- **Mean shift.** The relative ratio averages 1.059 over the 24 gated phylolm rows. Its SE is 0.009
  from the per-row MCSEs, or 0.012 from the between-regime SD (95% CI 1.035 to 1.084). So, relative
  to their own empirical SDs, the MI SEs are on average about 6% larger than the complete-data SEs.
  That is the conservative direction.
- **Between-regime heterogeneity is real.** The SD across regimes (0.058) exceeds the mean MCSE
  (0.044). Cochran Q = 41.0 on 23 df (p = 0.012); the estimated between-regime SD is 0.035. The
  excess sits in the twins (25 to 32: p = 0.033; 33 to 40: p = 0.014), not in the two-lambda
  regimes 17 to 24 (p = 0.50).
- **The count.** If every row shared the mean of 1.059, noise would put 0.61 rows outside
  [0.90, 1.15]; 3 were observed (P(3 or more) = 0.021). Rows 35, 36 and 38 sit 1.8 to 2.4 MCSE above
  the mean and 3.0 to 3.5 MCSE above 1. They are the upper end of a real spread, all on the
  conservative side.
- **The failing rows are well calibrated in absolute terms.** Their absolute MI SE ratio is 0.98 to
  1.06. They fail the relative rule because the complete-data phylolm ratio there is only 0.85 to
  0.90: the complete-data analysis is over-confident in these rows, and the MI analysis is not.
- **The other direction.** The plan's original absolute rule would fail 5 other rows (21, 23, 24,
  27 and 30; MI ratio 0.80 to 0.90). There the complete-data analysis is miscalibrated to the same
  degree (complete ratio 0.79 to 0.87 in those rows), and the relative ratio is 0.98 to 1.09.

**Decision for Shinichi (no gate changed).** The consequence of each option on the current numbers
(`se_ratio_noise.txt`, OPTIONS block):
- **(a) Keep the relative rule and accept G6 as failed on it.** 3 rows fail (35, 36, 38), all on the
  conservative side; report them with the reading above.
- **(b) Gate the absolute MI SE ratio where the complete-data analysis is itself miscalibrated**
  (complete ratio outside [0.90, 1.15]: 14 of the 24 rows), and the relative ratio elsewhere. Rows
  35, 36 and 38 then pass, but 5 rows fail (21, 23, 24, 27, 30), where the MI and complete-data SEs
  are both too small by about the same amount. A variant that passes a row when either ratio is in
  the band fails no row.
- **(c) Other rules**, with their outcomes:
  - the relative ratio within 1 +/- 2.5 correlation-aware MCSE per row: 6 rows fail (25, 26, 29, 35,
    36, 38). The earlier version of this section proposed this rule; it looked like a pass only
    because of the inflated MCSE;
  - the pooled mean relative ratio in [0.95, 1.10]: 1.059 passes, and its 95% CI (1.035 to 1.084)
    lies inside the band. This rule cannot see a single bad regime.

The gate is unchanged pending your decision.

## Stress test: regimes 1 to 16 (raw covariance of non-ultrametric trees; reported only)

G6 rule outcomes (`04_acceptance.R`, STRESS block): 4 violations in 4 of 16 regimes.
- Paired bias under phylolm in regime 1 (-0.026) and regime 9 (-0.020).
- Relative SE ratio in regime 8 (0.879) and regime 15 (0.897).
- Regime 3 under phylolm has a coverage shortfall of exactly 10/200 = 0.050, so it meets the rule
  "coverage >= complete - 0.05". The earlier output counted it as a violation through a
  floating-point comparison; `04_acceptance.R` now allows 1e-9 there, as its other rules do.

G7 per-cell coverage passes in every stress regime (0.943 to 0.957 MCAR). `diagnosis.md` explains
the bias:
- for phylolm, the dominant part is the model mismatch (tip-variance heterogeneity that `cov2cor`
  removes);
- over the 16 regime pairs the in-model twins remove most of the phylolm part: mean paired bias
  -0.0106 here against -0.0030 in the twins (72%); twin range -0.011 to +0.003, against -0.026 to
  +0.001 here;
- at lambda = 1 the twins remove about half of it (54% on average; from 16% in regime 13 to 75% in
  regime 3). The remaining lambda = 1 bias is systematic (Headline 2);
- for gls the twins remove none of the lambda = 1 bias (mean -0.0067 here against -0.0090 in the
  twins); `diagnosis.md` found the model-family component small for gls.

For ultrametric (dated) trees `vcv(tree)` is proportional to `cov2cor(vcv(tree))`, so this mismatch
does not arise.

## How we got here (round 1)

The first campaign (24 regimes, commit 69670d4) failed G6 on 8 checks (`04_acceptance.R` on the
first-campaign `sim_summary.csv` at commit f12ed09): non-convergence in three n = 300 regimes, the
lambda = 1 bias, one coverage shortfall and two SE ratios. `diagnosis.md`:
- **The lambda = 1 bias under phylolm** was mostly the harness DGP lying outside the model family.
  Exact imputations under the true covariance removed nearly all of it (about -0.003 remained at
  n = 300 with only x missing). In a 20-rep check on regime 3 the in-model twin removed 86% of the
  phylolm bias. The full 200-rep twin campaign gives 75% for that pair (regime 3 against its twin
  27) and 54% on average over the lambda = 1 regimes. The gls bias was not reduced.
- **The non-convergence** was an ESS shortfall only. Chains 4x longer converged all 32 non-converged
  fits in regimes 5, 21 and 23 without moving the slopes. The other 12 non-converged fits, in
  regimes 1 to 13, were first re-run in round 2, where all 44 converged.

Shinichi chose to add in-model twins and automatic chain extension (design.md 5e). The CP1 decisions
(design.md 5d) were also made before any result: the G3 calibration gate, the relative SE rule, and
reporting proper vs plug-in without gating it.

## Real data (preview: 9 of 10 cells; G8 pending the FishBase cell)

Source: `real_preview/` (tables from `script/mi_realdata/02_summarise.R`). All nine finished cells ran
at 69670d4, because G8 requires one code SHA (`real_preview/receipts.csv`).
- **Data:** PanTHERIA (4,027 species, 4 continuous traits; 3 MCAR and 3 clade-structured masks) and
  AVONET (1,500 species, 4 traits; 3 MCAR masks).
- **Convergence:** 7 of 9 cells converged. Two PanTHERIA cells (seed 20260819, MCAR and structured)
  reached min ESS 368 and 366 with R-hat at most 1.01; automatic extension would have extended them
  (`convergence_table.csv`).
- **Per-trait coverage** (`coverage_table.csv`, 36 trait-cells):
  - PanTHERIA (24 trait-cells): model-based 0.87 to 0.96, split conformal 0.90 to 0.98, Mondrian
    0.91 to 0.99; means 0.931, 0.948 and 0.953. Here the model-based intervals cover less than
    conformal (below split in 17 of 24 trait-cells).
  - AVONET (12 trait-cells): model-based 0.91 to 0.99, split 0.91 to 0.99, Mondrian 0.92 to 1.00;
    means 0.957, 0.959 and 0.962, about equal.
- **Slopes** (`slope_table.csv`, `pair_summary.csv`):
  - AVONET: pooled slopes within 3.4% of the complete-row reference (at most 1.7 reference SEs).
  - PanTHERIA: some cells differ by up to 4.0 reference SEs (body mass ~ head-body length, MCAR seed
    20260820: 2.952 against 2.858). Longevity ~ body mass differs by 3.3 reference SEs (0.148
    against 0.180). The body mass ~ head-body length shift is positive under all three MCAR masks
    (+2.5 to +4.0 reference SEs) and negative under all three structured masks (-0.5 to -3.1).
  - One candidate reason (untested until the `log_transform = FALSE` run): the PanTHERIA columns are
    already log values, and the default `log_transform = TRUE` logs them again. The imputation model
    is then linear on a log(log) scale while the pre-registered analysis is on the log scale, which
    is outside the stated congeniality scope. This alone does not explain why the sign of the body
    mass ~ head-body length shift depends on the mask type.
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
- Trees much beyond 10,000 tips. FishBase (10,484 tips, 5 traits) took 688 s for 600 sweeps in the
  smoke run, about 1 s per sweep including set-up (`evidence/smoke/README.md`), so a default fit
  takes hours.

The default draws method is unchanged (Shinichi's decision).
