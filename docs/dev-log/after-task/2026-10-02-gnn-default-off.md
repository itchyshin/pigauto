# After-task: `gnn = FALSE` becomes the default

2026-10-02. Branch `feat/gnn-default-off` (worktree `pigauto-gnn-default-off`), based on origin/main b565cad. Platform: Claude Code. Follow-up to PR #199 (defaults audit), on Shinichi's instruction "flip gnn to FALSE in a follow-up PR".

## 1. Goal

Make the calibrated phylogenetic baseline, with no GNN, the default for `impute()`, `fit_pigauto()`, `multi_impute()` and `multi_impute_trees()`, with a message pointing users with covariates to `gnn = TRUE`, a NEWS entry, a version bump and a stated way to restore the old behaviour.

## 2. Implemented

- Signature defaults `gnn = FALSE` in `impute()`, `fit_pigauto()`, `multi_impute()`, `multi_impute_trees()` and its two internal helpers.
- Roxygen for `gnn` in all four functions states the new default, why, and that `gnn = TRUE` restores the old behaviour and is needed for covariates.
- The covariate warning under `gnn = FALSE` now ends "Set gnn = TRUE to use them"; the Mondrian error now offers `gnn = TRUE`.
- `compare_methods()` and `simulate_benchmark()` gain `gnn = TRUE` arguments, so they keep comparing the baseline with the GNN-corrected model.
- `plot(fit)` on a fit without training history shows the gates plot instead of erroring; `plot(fit, type = "history")` on such a fit explains that no GNN was trained. The `plot.pigauto_fit` example sets `gnn = TRUE`.
- DESCRIPTION: version 0.11.0.9001; the Description says the default is the calibrated baseline and the GNN is optional.
- NEWS: breaking-change entry with the rerun evidence and how to restore the old behaviour.
- Vignettes: getting-started (how it works, the pipeline sentence, Step 5 now sets `gnn = TRUE`, the covariates example sets `gnn = TRUE`), common-pitfalls (gate section applies to `gnn = TRUE`), mixed-types (fitting section now uses the default and says what it does), gnn-architecture (note that the GNN runs only with `gnn = TRUE`).
- Tests: new default tests in `test-gnn-off.R` (signature defaults, default equals explicit `FALSE`, covariate warning text, `plot()` fallback); `gnn = TRUE` added to eight existing tests that exercise GNN features.
- Evidence copied from PR #199 so the citations resolve whichever PR merges first: `docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`, `script/campaign_gnn_rerun_results/*.csv`, `script/campaign_gnn_rerun_merge.R`, the one-line runner change in `script/campaign_gnn_off_cell.R`. Identical to #199's copies.

## 3a. Decisions and Rejected Alternatives

- Shinichi: flip to `FALSE`, with a covariate message; `"auto"` rejected earlier in the session.
- Mine: the covariate notice stays a warning (it already was one), with the fix appended, because covariates being ignored is a silent change in what the user asked for. `compare_methods()` and `simulate_benchmark()` keep the GNN on through an explicit argument rather than following the new default, because their documented purpose is the baseline against GNN comparison; `cross_validate()` follows the default because it evaluates the default model. torch stays in Imports (moving it is a separate proposal in the audit). Version bumped to 0.11.0.9001 (dev); the minor bump to 0.12.0 is a release decision.

## 4. Files Touched

Modified: `DESCRIPTION`, `NEWS.md`, `R/evaluate.R`, `R/fit_pigauto.R`, `R/impute.R`, `R/multi_impute.R`, `R/multi_impute_trees.R`, `R/plot.R`, `R/simulate_benchmark.R`, `man/compare_methods.Rd`, `man/fit_pigauto.Rd`, `man/impute.Rd`, `man/multi_impute.Rd`, `man/multi_impute_trees.Rd`, `man/plot.pigauto_fit.Rd`, `man/simulate_benchmark.Rd`, `script/campaign_gnn_off_cell.R`, `tests/testthat/test-fit-predict.R`, `tests/testthat/test-gnn-off.R`, `tests/testthat/test-mondrian-conformal.R`, `tests/testthat/test-new-features.R`, `tests/testthat/test-phase9-integration.R`, `tests/testthat/test-shipping-coverage.R`, `vignettes/common-pitfalls.Rmd`, `vignettes/getting-started.Rmd`, `vignettes/gnn-architecture.Rmd`, `vignettes/mixed-types.Rmd`.
Created: `docs/dev-log/after-task/2026-10-02-gnn-default-off.md`, `docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`, `script/campaign_gnn_rerun_merge.R`, `script/campaign_gnn_rerun_results/{agg_gnn_per_trait,agg_gnn_summary_arm,agg_solver_summary_arm}.csv`.

## 5. Checks Run

- First full suite after the flip (`NOT_CRAN=true`): 8 failing tests, all exercising GNN features (training history, DAE context, covariate threading, multi-obs refinement, transformer convergence, Mondrian). That run hit testthat's 10-failure cap, so it was rerun with the cap lifted after the fixes.
- After fixes, the six touched test files: `test-gnn-off.R` PASS 84, `test-fit-predict.R` 87, `test-mondrian-conformal.R` 26, `test-new-features.R` 111, `test-phase9-integration.R` 28, `test-shipping-coverage.R` 49; FAIL 0, SKIP 0 in each.
- Full suite and `rcmdcheck --as-cran`: see the addendum.
- `slop_check.py` on all added prose lines: 0 findings, 0 em dashes.

## 6. Tests of the Tests

- The new default tests fail on origin/main by construction (`formals(impute)$gnn` is `TRUE` there).
- The `plot()` fallback was found by the suite (`plot(obj$fit, type = "history")` failing), not anticipated; the new test checks both the fallback and the explanatory error.
- The covariate-warning tests muffle unrelated warnings and match on the new text, so they fail if the hint is removed.

## 7a. Issue Ledger

- Found while doing this: PR #199's rerun note claimed `.rds` aggregates were in the repo, but `*.rds` is git-ignored and only the CSVs were committed. Fixed on #199 (commit 10816a3) and here.
- `multi_impute_trees()` keeps `draws_method = "mc_dropout"` as its default, so its default call now prints the "mc_dropout draws are BM-posterior draws" message every time. Left as is; aligning that default is part of the separate `draws_method` proposal.

## 8. Consistency Audit

- Every exported caller of `impute()` / `fit_pigauto()` checked: `cross_validate()` follows the default; `compare_methods()` and `simulate_benchmark()` pinned to `gnn = TRUE`.
- Roxygen examples scanned for covariates, Mondrian, MC-dropout and history plots: only the `plot.pigauto_fit` example needed `gnn = TRUE`.
- README makes no GNN claim and needs no change in this PR.

## 9. What Did Not Go Smoothly

- The first test rerun stopped at 10 failures despite `options(testthat.progress.max_fails)`; the environment variable `TESTTHAT_MAX_FAILS` was needed.
- A multi-edit script aborted on a pattern that matched twice; redone with line-anchored edits that check each target line.

## 10. Known Residuals

- Downstream users and other lanes: any script that calls `impute()` without `gnn` gets different results after this merges. The `arc/imputation-sim` campaign lane reinstalls pigauto from main; its "pigauto GNN on" arm must pass `gnn = TRUE` explicitly.
- The rerun evidence is single-observation MCAR on four data sets; multi-observation data with covariates is where the GNN's `obs_refine` path matters and it was not re-benchmarked.
- torch is still a hard Import, so installation still needs it even though the default path makes no torch call.

## 11. Team Learning

- When a default flips, run the whole suite with the failure cap lifted before reading results; the first ten failures are not the population.
- A default that changes which object fields exist (here, training history) breaks generic methods such as `plot()`; check every S3 method that reads the changed fields.

## 12. Cross-Product Coverage

Covers: the four imputation entry points, `compare_methods()`, `simulate_benchmark()`, `cross_validate()` (by inheritance), `plot.pigauto_fit()`, the covariate and Mondrian messages, vignettes and DESCRIPTION.
Does NOT cover: moving torch to Suggests; `draws_method` defaults; `multi_impute_trees()` message noise; re-benchmarking multi-observation data with covariates; `pigauto_report()` wording beyond what the 2026-09-18 gnn-off work already handled.

## Section 5 addendum: final check lines

- Full suite, `TESTTHAT_MAX_FAILS=Inf NOT_CRAN=true devtools::test()`: FAIL 0, PASS 3099, SKIP 8 (all environmental: BIEN cache absent; rgbif, terra and smcfcs installed so their missing-package branches cannot run).
- `rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"))`: 0 errors, 0 warnings, 1 note (development version number; a pre-existing GitHub link in `common-pitfalls` returned a transient HTTP 503).
