# Rose review: feat/exact-default (predict_method = "auto" default)

Reviewer: Rose (fresh adversarial completion reviewer, D-43). I did not build this branch.
Scope: `git diff origin/main...HEAD -- R/ tests/ NEWS.md` (commits b29ac47, e4e2b6b on top of eccb298),
the evidence in this folder, and two cheap probes run with `devtools::load_all()` (described below).
Read-only apart from this file.

Overall verdict: NEEDS-CHANGES. The code is mostly sound and I found no leakage from test cells. The
problems are in the public claims (NEWS.md and the round-2 report) and in two documented guarantees the
code does not keep.

## 1. Validation leakage and double use of the validation split: WEAK

What holds (CONFIRMED):

- The route choice reads only `splits$val_idx` (R/fit_baseline.R:1328 to 1357). Test cells are never read.
- Both candidate fits are made with validation and test cells masked (R/fit_baseline.R:441 to 442;
  R/joint_threshold_baseline.R:22 to 23 and 41 to 42), so the scored cells were not used for fitting.
- Truth scales are right: squared error on the z-scored latent value for continuous, count, proportion,
  ordinal and zi magnitude; log-loss of `plogis(mu)` against 0/1 for binary and zi gate; negative log-prob
  of the true class for categorical (mu holds log-probs on both LP and OVR paths). Multi-obs maps
  observation rows to species rows correctly.
- The choice runs inside `fit_baseline()`, before gate calibration and conformal scoring in
  `fit_pigauto()`.

What does not hold:

- The same validation cells are used three times: to choose the route, to calibrate the gate, and to
  compute conformal scores. Choosing the route that minimises validation error makes the validation
  residuals optimistic, which shrinks the conformal quantile. The effect is small with two options and many
  cells, but largest where validation cells are few (the minimum is 5 cells per trait). Gate calibration
  already double-used the split before this branch, so this is not new in kind, but it is new in size.
- The round-2 table does not show "no problem". From `core_lambda_auto_off_agg_summary.csv` against
  `core_lambda_main_off_agg_summary.csv`, 27 of 72 trait-by-cell coverage values fall. At true lambda 1 and
  n = 100, all 8 fall (c1 0.889 to 0.883, c2 0.915 to 0.907, cnt 0.876 to 0.858, prp 0.892 to 0.878), and
  7 of 8 fall at n = 300. At n = 100 coverage is already below 0.95, so lower is worse. Each drop is within
  about 1.2 unpaired standard errors, but the sign is consistent across traits, which fits the
  selection-optimism mechanism (inference, not proven).
- The test at tests/testthat/test-route-choice.R:64 ("without hurting its validation MSE") scores the
  choice on the same validation cells that made it. It cannot fail except by a coding error, so it is not
  evidence that the choice helps.

## 2. Production refit, multi_impute, multi_impute_trees, predict(): CONFIRMED, with two edge cases

- `fit_pigauto()` (R/fit_pigauto.R:555 to 561) and `impute()` with `gnn = FALSE` (R/impute.R:541 to 550)
  replay `baseline$predict_method_by_trait` via `predict_route`. `multi_impute()` goes through `impute()`.
  `multi_impute_trees()` replays the reference route for `gnn = FALSE` and re-decides per tree on
  `splits_ref` for `gnn = TRUE` (R/multi_impute_trees.R:559 to 583). This is consistent.
- `predict()` uses `baseline_full` in production under `gnn = FALSE` and `baseline` otherwise
  (R/predict_pigauto.R:353 to 358). This is the same object that was fitted. CONFIRMED.
- Edge case A: under a concrete `predict_method = "exact"`, `predict_method_by_trait` is one aggregated
  value for every trait (R/fit_baseline.R:1237). If one joint fit falls back (for example one OVR class
  fit), every trait is labelled "per_column", and `baseline_full` is then forced to per-column for ALL
  traits while `baseline`, used for calibration, was mostly exact. The two baselines differ.
- Edge case B: under "auto", when exact falls back everywhere (oversized problem), both fits are
  identical, the losses tie, and every trait is reported as "exact" although per-column ran. The baseline
  is also fitted twice for no gain.

## 3. GLS-mean centring in the exact path: CONFIRMED

`.mvn_gls_mean_at_lambda()` (R/joint_mvn_solver.R:228 to 246) uses only non-NA cells of each column, and
`L` arrives with validation and test cells already NA, so the mean comes only from observed training
cells. It is on the latent (z or liability) scale at `lambda_block`. The variance returned is `ec$var`
unchanged (R/joint_mvn_solver.R:823), and the mean is added back to every cell. Small point: the dense
Cholesky per column (R/joint_mvn_solver.R:805) runs before `exact_conditional_mvn()` applies its size gate
(line 810), so oversized problems pay for centring and then fall back.

## 4. Tests pinned to predict_method = "per_column": mostly legitimate, one pin hides a real behaviour

Pins added in this branch (all in tests/testthat/): test-joint-lambda.R (3 tests), test-joint-refine-iter.R
(1), test-lambda-default.R (1), test-lambda-per-type.R (4), test-lambda-dispatch.R (5, plus a per_column
comparison inside "lambda_per_trait populated"), test-sigma-fisher-ml.R (1), and comparison arms in
test-exact-centring.R, test-exact-default.R and test-route-choice.R.

Legitimate per-column contracts: the ab02e31 reference fixtures, the equivalence to direct
`bm_impute_col()` calls, the discrete-stays-at-lambda-1 contract, and the partial `lambda_fixed` "other
columns match fixed_1" property. Several of these gained exact-route counterparts, which is good.

Problems:

- test-joint-refine-iter.R:38 to 53. The pin's own comment says that under exact, `joint_refine_iter` has no
  effect, because the exact branch returns before the EM loop. That was true on main too, but exact is now
  chosen for most traits by default, so `joint_refine_iter > 0` is silently ignored for those traits. The
  roxygen for `joint_refine_iter` (R/fit_baseline.R:178 to 184, and the fit_pigauto/impute copies) does not
  say this. This general property was pinned, not fixed or documented.
- test-lambda-dispatch.R "lambda_per_trait populated" asserts that b1 reports `lambda_block` under the
  default. That holds only if "auto" happens to choose exact for b1 on that fixture. It is fragile.
- test-lambda-dispatch.R "lambda_fixed rebuild" passes only because that fixture's routes do not mix. See
  finding 6.

## 5. Default-change reach: CONFIRMED

`impute()`, `fit_pigauto()`, `fit_baseline()`, `multi_impute()` (through `impute()`),
`multi_impute_trees()` (`baseline_arg("predict_method", "auto")`), `simulate_benchmark()`
(R/simulate_benchmark.R:106) and the cross-validation loop (R/evaluate.R:564) all reach "auto". The internal
helpers now default to "exact", but every call site I found passes `predict_method` explicitly
(R/joint_mvn_baseline.R:152, R/joint_threshold_baseline.R:391 and 642, R/ovr_categorical.R:131, 281, 319).
The EM paths now forward it, which matches the NEWS "Fixed" line.

## 6. Claims against evidence: WRONG

- "interval coverage is equal or higher": false. See finding 1: 27 of 72 coverage values fall, including
  every trait at lambda 1 and n = 100. The round-2 report's narrower sentence ("Coverage (c1, c2) is at or
  above main in every cell") is also false: c1 at lambda 1, n = 100 goes from 0.889 to 0.883 (rho 0) and
  0.882 (rho 0.5), and c2 from 0.915 to 0.907.
- "A baseline rebuilt from stored lambdas ... now reproduces the original fit under exact": only when every
  trait uses one route. Probe: 4 traits (2 continuous at lambda 0.3, binary, ordinal), n = 150, 12 seeds,
  default "auto", then `fit_baseline(..., lambda_fixed = bl$lambda_per_trait)` with the same splits. In
  every seed where a discrete trait chose per_column and another trait chose exact (7 of 12), the rebuild
  reproduced the same routes but not the same mu: max |d mu| from 0.019 to 0.16 on the latent scale. Cause:
  `.fit_baseline_auto()` combines `lambda_per_trait` from both fits, so a per_column discrete column carries
  lambda 1 into the rebuild's exact fit, where the original exact fit used `lambda_block`.
- Known limitation: NEWS says the ordinal drop happens "because the per-trait choice scores ordinal traits by
  squared error"; the round-2 report says "Likely cause". State it as likely.
- "Discrete accuracy rises by up to 0.044" leaves out that it falls at lambda 1 (-0.001 and -0.002 on the
  cell means, -0.0067 for ordinal at n = 1000).
- Headline z-RMSE numbers are traceable (3.1 to 9.8% at lambda 0.3 and 0.7, 0.1 to 1.2% at lambda 1; 16.8 to
  50.5% on AVONET and PanTHERIA; worst real case +0.1%). They are means over five traits. Per trait, 11 of 72
  cells get worse, up to +1.1% (cnt at n = 300, lambda 1: 0.727 to 0.734; prp at n = 100: 0.552 to 0.558).
  The mean hides this. It should be one clause.
- "18 core simulation cells": the round-2 table shows 9 rows, with rho 0 and 0.5 pooled. Say so.
- "1.4 to 1.5 times as long": traceable to S5b-route-choice-report.md:225 to 231, which records single
  unreplicated runs on two fixtures. Acceptable if worded as "about".

## 7. Package and CRAN reviewer concerns: WEAK

- The one-time fallback message still says `predict_method = "exact" (the default)`
  (R/joint_mvn_solver.R:71). The default is now "auto". It will mislead users.
- `predict_method_explicit` is a documented argument of the exported `fit_baseline()` marked "Internal use
  only; do not set directly". It works, but a package reviewer will object to an internal flag in the public
  signature. A non-exported wrapper, or an option read internally, would be cleaner.
- `predict_route` is documented in man/fit_baseline.Rd. Unknown names are silently ignored and missing names
  default to "exact" (R/fit_baseline.R:1440). A misspelled trait name gives no warning.
- Runtime: every default fit now makes two baseline fits (measured 1.42x and 1.50x). Examples and tests
  inherit this. I did not measure suite wall time. It is not a CRAN blocker at these sizes.
- The message is printed with `message()` regardless of `verbose = FALSE`. It can be suppressed, so this is
  acceptable, but consider respecting `verbose`.
- Evidence logs S5b-check.log and S5b-suite.log are one-line summaries, not full logs. They are newer than
  the R/ and tests/ files, so they cover this code, but the "--as-cran" status cannot be read from them.

## Required changes

1. NEWS.md and S5-benchmark-round2.md: remove "interval coverage is equal or higher" and the c1/c2 coverage
   sentence. Report the lambda 1 coverage drops (up to -0.018 at n = 100, below nominal) plainly.
2. Fix the `lambda_fixed` rebuild under mixed "auto" routes, or narrow the NEWS claim to single-route fits.
   A fix: store the exact fit's own `lambda_per_trait` (for example as a second attribute) so a rebuild
   replays each route's own lambdas. Add a mixed-route rebuild test; the probe above gives a fixture.
3. Update the fallback message text at R/joint_mvn_solver.R:71 to name "auto" as the default.
4. Document in the `joint_refine_iter` roxygen (fit_baseline, fit_pigauto, impute) that it has no effect on
   traits predicted by the exact route, so under the default it applies only to traits where "auto" picks
   per_column.
5. NEWS.md: state the ordinal cause as "likely"; add the per-trait z-RMSE worsening (up to 1.1% for count and
   proportion at lambda 1) and the small discrete-accuracy falls at lambda 1; say "9 settings by 2 rho values"
   rather than "18 cells" unless the table shows 18.
6. Either use a separate half of the validation split for the route choice (as `conformal_split_val` does
   for conformal), or state in the docs that route choice, gate calibration and conformal scores share the
   validation cells. Replace the tautological assertion at test-route-choice.R:64 with one scored on test
   cells.
7. Recommended, not blocking: make `predict_method_by_trait` report per-trait outcomes under a concrete
   "exact" (edge case A), report "per_column" when exact fell back under "auto" (edge case B), skip the
   duplicate fit when exact is known to fall back, and warn on unknown `predict_route` names.
