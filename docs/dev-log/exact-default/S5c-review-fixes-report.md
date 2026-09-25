# S5c: fixes to Rose's review of predict_method = "auto"

Branch: feat/exact-default (worktree /Users/z3437171/Dropbox/Github
Local/pigauto-exact-default). Fixes docs/dev-log/exact-default/rose-review.md
against the "auto" route-choice machinery from S5b
(docs/dev-log/exact-default/S5b-route-choice-report.md).

## Fix 1: split validation cells (route vs calibration/conformal)

Files: `R/fit_baseline.R` (`.pigauto_split_route_score()`, new;
`.fit_baseline_auto()`), `R/fit_pigauto.R`
(`.pigauto_calibration_val_idx()`, new; both `val_mask_mat`/`val_mask_score`
construction sites, gnn = TRUE around line 1274 and gnn = FALSE around line
601).

Under `predict_method = "auto"` with real validation cells,
`.pigauto_split_route_score()` splits each trait's validation cells into a
ROUTE half (route choice only) and a SCORE half, seeded via a new
`fit_baseline(..., seed = NULL)` argument (`seed + 29L + <trait index>`,
distinct from the existing `seed + 17L` / `seed + 23L` offsets used by
`calibrate_gates()` / `split_val_cal_conf()`). Below 10 total cells for a
trait, no split happens: the route falls back to `"exact"` (the existing
`n_total < 5` tie-break already does this once the trait sees zero route
cells) and ALL of that trait's cells stay available for scoring. The route
half is surfaced only internally (fed into
`.pigauto_choose_predict_route()`); the score half is returned as
`$score_val_idx`, and per-trait counts as `$route_val_n` / `$score_val_n`.

`fit_pigauto()` (both the `gnn = TRUE` and `gnn = FALSE` paths) and
`impute()` (via `fit_pigauto()`) now build a SEPARATE `val_mask_score`
matrix from `.pigauto_calibration_val_idx(splits, baseline, predict_method)`
and feed THAT (not the full `val_mask_mat`) into `split_val_cal_conf()`,
which then does its own cal/conf halving on top. The safety-floor mean and
`val_rmse` reporting still use the FULL `val_mask_mat` (every held-out cell
stays excluded from the training mean regardless of which half scored it)
-- only the cells `calibrate_gates()`/`compute_conformal_scores()` see are
restricted. Under an explicit `"exact"`/`"per_column"` request (no routing
decision made), `.pigauto_calibration_val_idx()` returns the full
validation set unchanged.

`model_config$route_val_n` / `$score_val_n` (both `fit_pigauto()` branches)
record the counts.

Example (n = 300, 5 mixed traits, `impute(..., gnn = FALSE, seed = 123)`):

```
predict_method_by_trait: c1=per_column c2=exact cnt=exact b1=per_column k1=exact
route_val_n:  c1=16 c2=14 cnt=16 b1=15 k1=47
score_val_n:  c1=15 c2=13 cnt=16 b1=14 k1=46
```

(k1 is categorical, K = 3, so its cell count spans 3 latent columns.)

Tests: `tests/testthat/test-route-choice.R` new blocks (g) (disjointness +
`model_config` fields) and (c) rewritten to score the chosen route on
`bl_auto$score_val_idx` instead of the cells that chose it.

## Fix 2: mixed-route lambda_fixed rebuild

Files: `R/fit_baseline.R` (`.pigauto_lambda_fixed_for_route()`, new;
`.fit_baseline_auto()`, `.fit_baseline_route()`).

`.fit_baseline_auto()` now attaches `attr(lambda_per_trait,
"lambda_by_route") <- list(exact = fit_exact$lambda_per_trait, per_column =
fit_pc$lambda_per_trait)` (the two UNMIXED, full-length lambda vectors that
actually produced `bl`), alongside the pre-existing `lambda_block`
attribute. `.pigauto_lambda_fixed_for_route(lambda_fixed, route)` reads
this attribute (when present) to pick the matching per-route vector instead
of the merged one; both `.fit_baseline_auto()` (for a plain `lambda_fixed =
bl$lambda_per_trait` rebuild under the default "auto") and
`.fit_baseline_route()` (for an explicit `predict_route` replay) now
resolve `lambda_fixed` through this helper before calling
`.fit_baseline_core()`. Root cause (per Shinichi's decision in the brief):
a mixed-route `bl$lambda_per_trait` is a MERGE of two different fits'
lambda values (e.g. a per_column-routed discrete trait's lambda pinned at 1
from `fit_pc`, vs `lambda_block` for the same-type column inside
`fit_exact`); feeding that merged vector uniformly into a fresh candidate
fit changes that candidate's Sigma estimate (via the per-column BM init
step) even for OTHER, unrelated traits sharing the same joint fit. Feeding
each candidate its OWN original vector reproduces it exactly.

Measured: `tests/testthat/test-route-choice.R` test (h), 4-trait fixture
(2 continuous at true lambda 0.3, binary, ordinal), n = 150, 12 seeds. 9 of
12 seeds produced a genuinely mixed route; all 9 reproduced
`predict_method_by_trait` identically and `mu`/`se` within 1e-8 (the
required tolerance) when rebuilt via `fit_baseline(pd, tree, spl,
lambda_fixed = bl$lambda_per_trait, seed = <same seed>)`. The same `seed`
must be passed to both calls -- the route/score split (fix 1) is the only
source of randomness in `fit_baseline()`, and matching it exactly is what
makes the ROUTE decision itself reproduce, not just the per-route lambda
values.

## Fix 3: predict_method_by_trait is now genuinely per-trait

Files: `R/fit_baseline.R` (`trait_route` tracker in `.fit_baseline_core()`,
set at every `col_path` joint-dispatch site: threshold_joint's
continuous/binary/ordinal populated-column loop, the ordinal
path-selection bm_mvn/lp override, the continuous-only joint_mvn branch,
and the OVR categorical loop; final `predict_method_by_trait` construction
now reads `trait_route` instead of broadcasting the single aggregate
`predict_method_used`). `.pigauto_choose_predict_route()`'s per-trait
default now comes from `fit_exact$predict_method_by_trait` (itself now
accurate) instead of a hardcoded `"exact"`, fixing Edge case B (a tie/
insufficient-evidence trait now reports what `fit_exact` actually ran for
it, not a blanket "exact").

Test: `tests/testthat/test-route-choice.R` (i) mocks
`exact_conditional_mvn()` to fail on exactly the SECOND call (the first of
two OVR class fits for a 3-level categorical trait `k1`, called after one
successful continuous-only joint fit for `c1`/`c2`) under an explicit
`predict_method = "exact"`. Before this fix, `c1`/`c2`/`k1` were ALL
mislabelled `"per_column"` (the whole-fit aggregate). After: `c1` and `c2`
report `"exact"`, only `k1` reports `"per_column"`.

## Fix 4: predict_method_explicit removed from the exported signature; predict_route warns

Files: `R/fit_baseline.R` (`fit_baseline()` now a thin wrapper with no
`predict_method_explicit` formal; new unexported `.fit_baseline_dispatch()`
carries the old body plus the `seed` argument), `R/fit_pigauto.R`,
`R/impute.R` (both now call `.fit_baseline_dispatch()` directly instead of
`fit_baseline()`, passing their own resolved `predict_method_explicit`).
`man/fit_baseline.Rd` no longer documents `predict_method_explicit`
(verified via `devtools::document()`).

`.fit_baseline_route()` now warns once per call
(`'predict_route' has name(s) that match no trait ... (ignored): ...`) on
any `predict_route` name absent from `data$trait_map`; values outside
`c("exact", "per_column")` already errored (unchanged). Roxygen for
`joint_refine_iter` (fit_baseline, fit_pigauto, impute) now states it has
no effect on any trait predicted by the exact route.

Test: `tests/testthat/test-route-choice.R` (j).

## Fix 5 (messages): joint_mvn_solver.R fallback text

`R/joint_mvn_solver.R:.pigauto_exact_fallback_notify()` no longer says
`predict_method = "exact" (the default)`; it now says the exact route "was
not usable ... either while predict_method = \"auto\" (the default) is
comparing its \"exact\" candidate, or when \"exact\" is itself the resolved
route and fell back." No test pinned the old string
(`grep -rn` confirmed).

## Deviations from the brief

- The review's `test-lambda-dispatch.R` "ordinal mu is finite and reports
  lambda_block under exact estimate" test (not "lambda_per_trait
  populated", a SEPARATE pre-existing test that was also fragile) called
  `fit_baseline(..., lambda_mode = "estimate")` with no `predict_method`,
  i.e. it ran under the new default "auto" while its own title and
  assertions are about the exact route's specific ordinal-lambda contract.
  Before fix 3, `fit_exact`'s aggregate mislabelling made "auto" default to
  "exact" for `o1` on a tie even when `fit_exact` itself had internally
  chosen bm_mvn/lp for that column (lambda 1, not lambda_block); fix 3
  correctly exposes that as "per_column", which broke this test on all 3
  of its seeds. Pinned with `predict_method = "exact"` (matching the
  test's own stated intent) rather than weakened.
- `route_val_n` / `score_val_n` model_config counts are recorded per S1's
  precedent (`lambda_per_trait`, `predict_method_by_trait`), sourced from
  `baseline` (never `baseline_full`, which has no validation cells of its
  own).
- `tests/testthat/test-fit-predict.R` "conformal_split_val=TRUE changes
  conformal scores (no double-dipping)" (a pre-existing C.3 regression
  test, not on the review's list) broke under the full-suite run: fix 1's
  restriction of calibration/conformal cells to the SCORE half roughly
  halves the ~40 val cells/column this test relies on, pushing both fits
  below `2 * min_val_cells = 40`, so `conformal_split_val`'s own
  per-column split stopped firing for EITHER fit and they became
  identical for a reason unrelated to what the test checks. Pinned
  `predict_method = "per_column"` on both `fit_pigauto()` calls (no
  routing decision is made, so the full validation set is used,
  restoring the original cell-count assumption) -- this test is about
  the C.3 split, not about "auto".

## VERIFY

SUITE FAIL 0 | PASS 2709
CHECK errors 0 warnings 0 notes 1 (CRAN-incoming-feasibility dev-version-number note only, pre-existing)
