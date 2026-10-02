# After-task: defaults audit, documentation cleanup, multiple-imputation integration

2026-10-02. Branch `docs/defaults-audit` (worktree `pigauto-defaults-audit`), based on origin/main b565cad. Platform: Claude Code.

## 1. Goal

A user should be able to use pigauto's defaults well without knowing the alternatives (Shinichi, 2026-10-02). Three parts: audit every user-facing argument with alternatives and propose a verdict for each; lead the documentation with the defaults and caveat the alternatives; make `multi_impute()` -> `with_imputations()` -> `pool_mi()` work cleanly with worked examples for glmmTMB, lme4, brms, drmTMB and gllvmTMB. No default changes and no deprecations; proposals only.

## 2. Implemented

- `docs/dev-log/defaults-audit/2026-10-02-defaults-audit.md`: nine functions, one table each (default, alternatives, evidence with file:line, proposed verdict), five decisions and proposals, ten inconsistencies.
- `docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`: rerun of the 2026-09-19 with/without-GNN campaign and the joint-solver diagnostic at current main on Totoro (600 cell runs, 0 errors). Aggregates in `script/campaign_gnn_rerun_results/`.
- `pool_mi()` refuses `brmsfit` objects with a pointer to `brms::brm_multiple()` / `brms::combine_models()`; roxygen lists supported classes and explains posterior concatenation for Bayesian fits.
- `tests/testthat/test-pool-mi-backends.R`: end-to-end posterior MI -> `with_imputations()` -> `pool_mi()` for glmmTMB, lme4, drmTMB and gllvmTMB, plus the brms refusal.
- `vignettes/multiple-imputation.Rmd` (new pkgdown article): which draws can be pooled, worked examples for nlme, lme4, glmmTMB, drmTMB, gllvmTMB, brms (`brm_multiple`) and MCMCglmm (stacked `Sol`).
- README: "Defaults, and when to change them" section and article link. `_pkgdown.yml`: article entry. NEWS entry.
- `fit_pigauto()` roxygen: `gate_method` default is `"cv_folds"` (was documented as `"single_split"`); `min_val_cells` default is 20 (was documented as 10), with a corrected explanation of the conformal quantile.
- DESCRIPTION Suggests: brms, drmTMB, gllvmTMB.
- `script/campaign_gnn_off_cell.R`: the derived full-baseline arm now replays the per-trait `predict_method` route; `script/campaign_gnn_rerun_merge.R` (new).

## 3a. Decisions and Rejected Alternatives

Shinichi's decisions (2026-10-02, in chat):
- `draws_method` default unchanged for now; the audit proposes a follow-up `draws_method = "auto"`.
- `gnn`: the audit proposes `gnn = FALSE` by default with a message when covariates are supplied. Rejected: `gnn = "auto"` (on when covariates are given), because it hides a large runtime difference and a torch requirement and there is no evidence the GNN uses covariates well. Earlier pick "keep TRUE, report both" was revisited by Shinichi the same session.
- brms is refused in `pool_mi()`; no new exported concatenation helper (rejected as surface to maintain).
- brms, drmTMB, gllvmTMB into Suggests (gllvmTMB is not yet on CRAN; see residuals).
- The GNN-on rerun used a private torch install on Totoro rather than repairing the shared library (Shinichi's choice), so another lane's runs were not touched.

My calls: chunks in the vignette are `eval = FALSE` like the other vignettes, because the posterior sampler takes about ten minutes; BACE was dropped from the rerun (its numbers are not a defaults question and it was the costliest arm).

## 4. Files Touched

Modified: `DESCRIPTION`, `NEWS.md`, `R/fit_pigauto.R`, `R/pool_mi.R`, `README.md`, `_pkgdown.yml`, `man/fit_pigauto.Rd`, `man/pool_mi.Rd`, `script/campaign_gnn_off_cell.R`.
Created: `docs/dev-log/defaults-audit/2026-10-02-defaults-audit.md`, `docs/dev-log/arc/2026-10-02-campaign-gnn-rerun.md`, `docs/dev-log/after-task/2026-10-02-defaults-audit.md`, `script/campaign_gnn_rerun_merge.R`, `script/campaign_gnn_rerun_results/{agg_gnn_summary_arm.csv, agg_gnn_per_trait.csv, agg_solver_summary_arm.csv}` (the `.rds` aggregates stay on Totoro; `*.rds` is git-ignored), `tests/testthat/test-pool-mi-backends.R`, `vignettes/multiple-imputation.Rmd`.
Not tracked: `.unlazy/defaults-audit/` (gate ledger). Totoro: `~/defaults-audit/` (private library, raw cells, logs).

## 5. Checks Run

- `NOT_CRAN=true devtools::test()` (full suite): see section 5 addendum below for the final line.
- `rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"))`: see addendum.
- `testthat::test_file("tests/testthat/test-pool-mi-backends.R")` with `NOT_CRAN=true`: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 25 ]`. `test-multi-impute.R`: `[ FAIL 0 | WARN 32 | SKIP 0 | PASS 202 ]` (warnings from `stats::cor` in `evaluate_imputation()`, pre-existing).
- Vignette code purled and run end to end with shortened chains (`m = 3`, 2 chains, 200 burn-in, 400 iterations): every chunk ran, including `brm_multiple` and MCMCglmm; the only warning was the expected convergence warning from the short chains, plus a gllvmTMB 0.2.0 message about `latent()` including a unique variance by default.
- Gate ledger `.unlazy/defaults-audit/gates/leaf-main.md`: G1 audit covers nine functions, G23 backend tests 0 fail 0 skip, G411 `draws_method` and `gnn` defaults unchanged against origin/main: all met.
- `slop_check.py`: 0 findings on the audit, the rerun note, the vignette, README and this report.
- Totoro rerun: 400 pass-1 cell runs and 200 pass-2 cells, 0 errors; smoke cells read before each launch.

## 6. Tests of the Tests

- The brms refusal test failed before the stop was added (the old error came from the coefficient validator), then passed.
- Added `riv > 0` to the backend assertions after review: pooling a single fit, or identical datasets, gives `riv = 0`, so the tests now fail if between-imputation variance is lost.
- The first gate scripts were wrong twice (relative paths against the ledger folder; `$` expanded by the shell inside a double-quoted `Rscript -e`). Both were fixed in the gate, and G411 was checked to print the touched lines when a default line differs.
- The gllvmTMB backend test is a plumbing test only (jittered pseudo-replicates); it says so in a comment.

## 7a. Issue Ledger

From the claim-vs-evidence review (fresh context, Opus): 0 blocking, 15 required, 7 suggestions. All 15 required items fixed: `pool_mi()` does not refuse a hand-built list (docs now say `with_imputations()` refuses the draws); audit numbers (posterior in-model bias range, coverage regime, "two" not "three", runtime regime); line-number basis stated; stale solver rows updated from the rerun; ranking labelled as the auditor's proposal; two roxygen errors (`cv_folds` default date, the conformal-quantile explanation); vignette overgeneralisation of attenuation and of per-cell coverage; drmTMB `sigma` caveat (marked as not measured); runtime claim given its regime; `git add -f` for ignored `docs/`; README "sound answer" softened. Suggestions S1 to S5 and S7 done; S6 (add the Rubin-campaign limits to the vignette once `arc/rubin-freq-bace` lands) deferred.

## 8. Consistency Audit

- `with_imputations()` already refused conformal and MC-dropout output on main (PR #174); the audit and docs were reframed around that, including the usability trap it creates for the default route.
- The 2026-09-19 claim that the AVONET300 gap came from the liability step had already been withdrawn by PR #182; the audit cites the solver result instead.
- Other doc/default mismatches found and listed but not fixed here: `impute()` epochs 2000 against `fit_pigauto()` 3000; `phylo_signal_method` not validated in `impute()`; `multi_impute_trees()` defaults to `"mc_dropout"` and offers no `"posterior"`; old NEWS entries.
- A pre-existing roxygen link warning (`[0, 1]` in `R/bm_internal.R`) and the README's experimental banner em dash were left as they are.

## 9. What Did Not Go Smoothly

- The Totoro torch runtime was gone; the smoke cell caught it (GNN-on arm errored in 3 s). Installing torch and libtorch into a private library fixed it without touching the shared one.
- The audit agent briefly reported the brms refusal as already on main; it was reading the other agent's uncommitted edit in the same worktree.
- My own first draft of three origin/main line numbers was wrong and was corrected by reading `git show origin/main:R/pool_mi.R`.

## 10. Known Residuals

- **gllvmTMB is in Suggests but not on CRAN**, and there is no `Additional_repositories`. pigauto's next CRAN submission must wait until gllvmTMB is on CRAN, or drop it from Suggests.
- The vignette's chunks are not evaluated at build; they were run once with shortened chains, not with the default sampler.
- The rerun covers single-observation MCAR data on four DGPs; multi-observation, covariates, MNAR, `zi_count` and `multi_proportion` are not covered, and the GNN-on arms were not run with the Rphylopars solver.
- The `safety_floor` / `phylo_signal_gate` finding is simulation-only.
- Evidence for conformal-draw attenuation and the posterior n = 1000, lambda = 1 under-coverage lives on unmerged branches.

## 11. Team Learning

- Before writing user docs about a refusal, find which function refuses: here `with_imputations()` refuses and `pool_mi()` only warns on a bare list.
- A smoke cell per arm, read before launch, pays for itself on a shared server whose libraries change under you.
- When two agents share a worktree, tell the reviewer which uncommitted edits belong to whom.

## 12. Cross-Product Coverage

Covers: the nine named functions' user-facing arguments; `pool_mi()` with lm, gls, lme4, glmmTMB, drmTMB, gllvmTMB (tested end to end for four), brms and MCMCglmm (refused); `multi_impute(draws_method = "posterior")` as the inference route.
Does NOT cover: `multi_impute_analysis()`, `cross_validate()`, `simulate_benchmark()`, `compare_methods()` arguments; GNN hyper-parameters; the posterior route for discrete traits (not implemented); phylolm in `pool_mi()` (generic fallback only, untested); brms fits through `pool_mi()` other than the refusal; any default change.

## Section 5 addendum: final check lines

- Full suite, `NOT_CRAN=true devtools::test()` with `SummaryReporter`: no failures (no "Failed" section), 8 skips, all environmental (BIEN cache absent in 5; rgbif, terra and smcfcs installed so their missing-package branches cannot run). Run on the tree before the review fixes; after them `test-pool-mi-backends.R` was rerun alone (`PASS 25`) and the only other changes were roxygen text.
- `rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"))`: 0 errors, 0 warnings, 1 note. The note lists the development version number, `gllvmTMB` as a Suggests not in mainstream repositories, and the new article URL returning 404 until pkgdown is redeployed.
