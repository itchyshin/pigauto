# After-task report: `predict_method = "auto"` as the default (lane "go exact")

Date: 2026-09-26. Platform: Claude Code. PRs: #192 (merged bf7ca7a), #193 (merged 7608725).

## 1. Goal

Make cross-trait "exact" prediction pigauto's default baseline route, provided it is not worse than the previous
per-column default on the core simulation cells and the real-data cases, and document what it costs.

## 2. Implemented

- GLS-mean centring in the exact conditional solve, and a quiet once-per-session fallback to per-column above
  about 20,000 unknown cells.
- `predict_method = "auto"` (new default on `impute()`, `fit_pigauto()`, `fit_baseline()`): fits both routes and
  keeps, per trait, the one with lower validation loss. The choice is stored in
  `model_config$predict_method_by_trait` and replayed on the production refit via `predict_route`.
- Validation split between route choice and calibration, by held-out trait row (species in multi-obs data),
  only when a trait has at least 38 rows and the two fits differ.
- Under exact, discrete traits share `lambda_block`; fix for `em_iterations >= 1` ignoring `predict_method`.
- NEWS, roxygen, man pages and README updated; evidence in `docs/dev-log/exact-default/` and
  `docs/dev-log/ordinal-route/`.

## 3a. Decisions and Rejected Alternatives

- Exact alone as the default (round 1): rejected, z-RMSE +2.0 to +2.6% at lambda 1.
- Route choice on the same cells that calibrate (round 2): rejected after review, the choice biases the
  conformal residuals.
- Split every trait with at least 10 cells (round 3): rejected, coverage fell up to 0.032 at n <= 300.
- Split by latent cell (round 4): rejected after the second review, categorical rows leaked between halves.
- Ordinal route scoring by class error rate or by log-probability of the true class: tested, both lowered
  ordinal accuracy; squared error kept. Decisions on each step were Shinichi's (recorded in the plan file).

## 4. Files Touched

R: `fit_baseline.R`, `fit_pigauto.R`, `impute.R`, `joint_mvn_baseline.R`, `joint_mvn_solver.R`,
`joint_threshold_baseline.R`, `multi_impute.R`, `multi_impute_trees.R`, `ovr_categorical.R`.
man: `fit_baseline.Rd`, `fit_pigauto.Rd`, `impute.Rd`, `multi_impute.Rd`.
tests: `test-exact-centring.R`, `test-exact-default.R`, `test-route-choice.R` (new); `test-exact-conditional.R`,
`test-fit-predict.R`, `test-joint-lambda.R`, `test-joint-refine-iter.R`, `test-lambda-default.R`,
`test-lambda-dispatch.R`, `test-lambda-per-type.R`, `test-sigma-fisher-ml.R` (pins updated).
Docs: `NEWS.md`, `README.md`, `docs/dev-log/exact-default/*`, `docs/dev-log/ordinal-route/*`.
Unmerged provenance branches: `feat/ordinal-route-accuracy` (buggy first ordinal run),
`feat/ordinal-route-logloss` (benchmark switch).

## 5. Checks Run

- Full suite on the merged code: FAIL 0 / PASS 2735 (`S5e-suite.log`).
- `R CMD check --as-cran --no-manual`: 0 errors, 0 warnings, 1 note (development version) (`S5e-check.log`).
- CI on #192: macOS release, Ubuntu release and devel all passed; main post-merge R-CMD-check and pkgdown
  passed at bf7ca7a.
- Benchmarks on Totoro, GNN off, 200 seeds, 18 cells: rounds 1 to 5 plus ordinal runs (3,600 jobs each);
  13 real-data cases at 5 seeds.

## 6. Tests of the Tests

- The row-split test uses a categorical trait with at least 38 held-out rows and asserts no row appears in both
  halves; before the fix, the cell split failed this on AVONET (Trophic.Level calibrated on 10 of 29 rows).
- The no-choice test asserts a single-trait fit keeps every validation cell for calibration; before the fix
  its conformal score was 1.929 instead of 1.180.
- A mixed-route rebuild test replays stored lambdas over 12 seeds and requires mu and se to match within 1e-10.

## 7a. Issue Ledger

| issue | found by | status |
|---|---|---|
| validation cells used twice | first review | fixed (d2b0c57) |
| mixed-route rebuild did not reproduce | first review | fixed (d2b0c57) |
| split starved conformal at n <= 300 | round 3 | fixed (13bd3e6) |
| categorical rows leaked across halves | second review | fixed (2be557f) |
| needless split when routes agree | second review | fixed (2be557f) |
| refit warned every time under default | second review | fixed (2be557f) |
| NEWS overclaims | both reviews | fixed |
| ordinal scorer decoded 1..K | own check | re-run with 0..K-1 |
| LP ordinal candidate decodes 1..K | own check | open, not in scope |

## 8. Consistency Audit

Checked every caller of the split and route functions (`fit_pigauto`, `impute`, `multi_impute`,
`multi_impute_trees`, `cross_validate`); roxygen for `lambda_mode`, `predict_method`, `predict_route` and
`seed`; README examples; every number in NEWS, the round-4/5 notes and the PR bodies against the CSVs.

## 9. What Did Not Go Smoothly

- My own first NEWS draft claimed coverage was "equal or higher"; it was not.
- The first ordinal test clamped classes to 1..K and produced a false negative that reached a PR (#193) before
  I caught it; the PR was put back to draft and corrected before merging.
- Remote `pkill -f` patterns matched their own SSH command twice.
- #193 merged as soon as it was marked ready because docs-only changes run no R check.

## 10. Known Residuals

- Coverage at lambda 1 with 100 species falls by 0.006 to 0.018 per trait (main already 0.86 to 0.92).
- z-RMSE +0.4% at lambda 1 with 1,000 species; GlobTherm full +0.8%.
- The label-propagation ordinal candidate drops lowest-class observations (commit a3b89e6); unmeasured.
- `multi_impute_trees()` offsets the seed per tree, so tree 1 does not reproduce the reference fit's split.

## 11. Team Learning

- Check an encoding against the preprocessing docs before scoring on decoded values; the 0..K-1 ordinal coding
  is documented in `preprocess_traits()`.
- A validation split must use the unit the holdout uses (the trait row), not the storage unit (the latent cell).

## 12. Cross-Product Coverage

Covers: single-observation mixed-type data (continuous, count, proportion, binary, three-class categorical,
ordinal) under Brownian motion with Pagel's lambda 0.3, 0.7 and 1, n 100 to 1,000, 30% MCAR, GNN off; and 13
continuous-trait real-data cases, GNN off.

This lane does NOT cover:
- the GNN-on default of `impute()` (not re-benchmarked);
- multi-observation data (unit-tested only, not benchmarked);
- non-Brownian processes (OU, early burst) and missing-not-at-random masks;
- zero-inflated count and multi-proportion traits, which are not in the core cells;
- Monte Carlo error on the real-data runs (5 seeds only);
- the Rphylopars joint solver arm, which was not re-run under `"auto"`.
