# S5b: per-trait predict_method route choice ("auto")

Branch: feat/exact-default (worktree /Users/z3437171/Dropbox/Github
Local/pigauto-exact-default). Not committed, not pushed, per instructions.

STATUS AT HANDBACK: implementation complete, own tests written and green,
the 28 pre-existing files that mention `predict_method`/`lambda` re-run
green individually. **The full-suite run and `rcmdcheck` were forced to
hand back before completion; this is NOT a green-suite claim.** See
"Outstanding verification" below for the exact resume commands and log
paths.

## What changed

### R/fit_baseline.R
- `fit_baseline()` (exported) is now a thin dispatcher. Its own body
  (unchanged code) was renamed `.fit_baseline_core()` (not exported,
  `predict_method` restricted to `c("exact", "per_column")`); it gained one
  new field on its return value, `predict_method_by_trait` (every trait
  name mapped to that single fit's own aggregate `predict_method_used`).
- New exported `fit_baseline()` signature: `predict_method = c("auto",
  "exact", "per_column")` (default flipped from `"exact"` to `"auto"`) and
  a new `predict_route = NULL` argument (named character vector, values
  `"exact"`/`"per_column"`, one entry per trait; forces the route,
  bypassing `"auto"` and overriding `predict_method`).
- New internal helpers:
  - `.pigauto_route_val_loss(tm, val_row, species_row, X_truth, mu)`:
    per-trait validation loss for one route's `mu`. Continuous / count /
    proportion / ordinal -> MSE on the z-scored latent scale. Binary ->
    mean log-loss of `plogis(mu)` against 0/1 truth. Categorical -> mean
    multinomial log-loss (true class from the one-hot truth,
    `-log_prob` at that class). zi_count -> gate log-loss + magnitude MSE,
    summed, with `n` = gate cells + magnitude cells. Anything else
    (multi_proportion, which never reaches a joint fit -- both routes are
    identical for it -- or any other type) -> `n = 0`, no comparison
    possible.
  - `.pigauto_choose_predict_route(data, splits, fit_exact, fit_pc)`:
    decodes `splits$val_idx` into (row, col) the same way the pre-existing
    per-trait ordinal path-selection block does (`n_rows_sp <-
    nrow(data$X_scaled)`, `val_col`/`val_row` via integer division), maps
    obs row to species row via `data$obs_to_species` for multi-obs, then
    for each `trait_map` entry compares `.pigauto_route_val_loss()` on
    `fit_exact$mu` vs `fit_pc$mu`. Ties, non-finite losses, or fewer than 5
    validation cells (`max(n_exact, n_pc) < 5`) -> `"exact"` (already the
    vector's default); `per_column` strictly lower -> `"per_column"`.
  - `.fit_baseline_auto(...)`: no `splits`/no val cells -> delegates to
    `.fit_baseline_core(predict_method = "exact", predict_method_explicit
    = FALSE)` once (byte-identical fallback behaviour to the pre-S5b
    default), stamps `predict_method_used <- "auto"` on the result but
    keeps that core call's own `predict_method_by_trait` (so a trait whose
    exact route is itself unusable, e.g. K < 2, still shows its real
    fallback outcome, not a hardcoded `"exact"`). With val cells:
    fits `.fit_baseline_core()` twice (`predict_method = "exact"` and
    `"per_column"`, both `predict_method_explicit = FALSE` so an internal
    fallback message-once, not a warning), calls
    `.pigauto_choose_predict_route()`, then builds the output by copying
    each trait's `mu`/`se`/`path`/`lambda_per_trait` columns from whichever
    of the two full fits was chosen (reuses the two existing code paths
    verbatim; no solver duplication).
  - `.fit_baseline_route(...)`: resolves `predict_route` against
    `trait_map` (unnamed trait defaults to `"exact"`), fits `.fit_baseline_core()`
    once per DISTINCT route actually needed (once if `predict_route` is
    uniform, twice if mixed), then mixes trait columns the same way as
    `.fit_baseline_auto()`.
- Roxygen: `predict_method`'s `@param` documents `"auto"` (default) fully
  (per-trait loss types, tie-break, no-split fallback); `predict_route` is
  a new `@param`; `@return` documents the new `predict_method_by_trait`
  field and that `predict_method_used` can now be `"auto"`.

### R/fit_pigauto.R
- `predict_method = c("auto", "exact", "per_column")` (default flipped).
- The `baseline_full <- fit_baseline(..., splits = NULL, ...)` call (inside
  the `gnn = FALSE` branch) now passes `predict_route =
  baseline$predict_method_by_trait`: this refit has no validation cells of
  its own, so it REUSES the choice `baseline` (the splits-based fit, just
  above) made on real evidence, instead of `"auto"` silently defaulting
  every trait to `"exact"` from zero evidence. Because a forced
  `predict_route` always reports `predict_method_used = "auto"` inside
  `fit_baseline()` (it cannot tell "genuine auto request" apart from "a
  uniform forced route" from where it's called), the call site then
  restores `baseline_full$predict_method_used <- baseline$predict_method_used`
  so an explicit `predict_method = "exact"`/`"per_column"` request is not
  misreported as `"auto"` once replayed onto `baseline_full` -- this was
  caught by the smoke test in "Verification" below (initially failed:
  explicit `"exact"` request reported `model_config$predict_method_used ==
  "auto"`; fixed by this relabel).
- `model_config` (both the `gnn = FALSE` and `gnn = TRUE` branches) gained
  `predict_method_by_trait`, sourced with the same `baseline_full`-over-
  `baseline` precedence as `lambda_per_trait`/`lambda_block`/
  `predict_method_used` for `gnn = FALSE`, and from `baseline` alone for
  `gnn = TRUE` (predict() always uses `object$baseline` for GNN-on fits,
  never `baseline_full` -- confirmed by reading `R/predict_pigauto.R`
  lines ~330-360 before making this change).
- Roxygen for `predict_method` rewritten to describe `"auto"`.

### R/impute.R
- Same enum flip. `impute()` builds `baseline` and `baseline_full` ITSELF
  (not inside `fit_pigauto()`, since it passes both in pre-computed) -- so
  the identical `predict_route` reuse + `predict_method_used` relabel fix
  had to be duplicated here, at `impute()`'s own `baseline_full <-
  fit_baseline(pd, tree, splits = NULL, ...)` call site. Roxygen updated
  identically to `fit_pigauto()`'s.

### R/multi_impute_trees.R (`run_shared_gnn()`)
- `baseline_arg("predict_method", "exact")` fallback default flipped to
  `"auto"` (this is the "multi_impute_trees fallback" named in the brief
  as one of the dispatcher default-flip points).
- New `predict_route_arg`: for `gnn = FALSE` (this function's per-tree
  baseline replay uses `splits = NULL` there, matching
  `baseline_full`/production semantics), reuses
  `fit_ref$model_config$predict_method_by_trait` (via the existing
  `baseline_arg()` config-fallback helper) exactly as the `fit_pigauto()`/
  `impute()` fix above. For `gnn = TRUE` (this function uses `splits_ref`,
  i.e. real validation cells, for the per-tree replay), `predict_route` is
  left `NULL` so `"auto"` re-decides per tree on its own real evidence, as
  the existing per-tree-refit design already does for every other option
  (`lambda_mode`, `joint_solver`, etc.).

### R/joint_mvn_solver.R, R/joint_mvn_baseline.R, R/joint_threshold_baseline.R, R/ovr_categorical.R
- **Not touched.** Per the brief: "internal solver functions can keep
  'exact' as their own default, only the dispatcher needs 'auto'." These
  are only ever called with an explicit, already-resolved `predict_method`
  from `fit_baseline()`/`.fit_baseline_core()`, so their own defaults are
  dead code paths in practice; changing them was out of scope and would
  have been a no-op.

### tests/testthat/test-route-choice.R (new)
Seven `test_that()` blocks (`[route] ...`), covering (a)-(f) from the
brief (test (a) folded into one block with both formals checks and the
default-behaviour smoke, matching the S3 test file's own precedent):
1. Defaults are `"auto"` at `fit_baseline`/`fit_pigauto`/`impute`
   (`formals()`, no fit).
2. DGP (b): lambda = 0.7, 3 correlated (rho = 0.6) continuous traits,
   n = 200, simulated directly from the exact route's own generative
   model (`vec(L) ~ MVN(0, Sigma %x% R(lambda))`, via Cholesky:
   `L = t(chol(R_lambda)) %*% Z %*% chol(Sigma)`) -- auto chooses `"exact"`
   for all 3.
3. DGP (c): 2 independent continuous traits at true lambda = 1, plus a
   near-white-noise (lambda = 0.05) count trait -- confirms
   `lambda_block < 1` (the count pulls it down), then asserts `"auto"`
   picks `"per_column"` for at least one continuous trait AND that trait's
   auto MSE is never worse than the route it didn't pick. See "Route
   choices and MSEs" below for the numbers.
4. `predict_route` forces the route and reproduces both a uniform-exact,
   uniform-per_column, AND a mixed-route forced fit within `1e-10` of the
   matching `.fit_baseline_core()` calls.
5. The `fit_pigauto()` production refit (`impute(..., gnn = FALSE)`)
   reuses the validation-split route: `baseline_full$predict_method_by_trait`
   is identical to `baseline$predict_method_by_trait`, and
   `model_config$predict_method_by_trait` matches both.
6. `splits = NULL`, and separately a `splits` object with `val_idx =
   integer(0)`, both give every trait `"exact"`.

Also fixed 5 pre-existing assertions in `tests/testthat/test-exact-default.R`
that assumed the old `"exact"` default (per the brief's pin-or-generalize
rule):
- Test 1 (`[exact-default] ... default to predict_method = 'exact'`):
  this test IS about the default value itself -- a genuine behaviour
  change, not a per-column-specific quirk -- so its assertion was updated
  to expect `"auto"` (not pinned).
- Test 2 (`model_config$predict_method_used is 'exact' on a small fit`):
  pinned, added `predict_method = "exact"` to the `impute()` call -- this
  test is about the exact route's own bookkeeping, not about what "auto"
  happens to pick on that particular random draw.
- Test 3 (`oversized/mocked exact under the default: message, ...`): left
  UNPINNED (no `predict_method` argument) because its whole point is the
  DEFAULT's fallback message-vs-warning behaviour, and with `splits = NULL`
  (the fixture here has none), `"auto"`'s no-val branch runs the exact
  route internally with `predict_method_explicit = FALSE` -- byte-for-byte
  the same message behaviour as the old default. Updated its two
  assertions: `bl$predict_method_used` is now `"auto"` (the constant
  top-level marker) rather than `"per_column"`; the real per-trait outcome
  of the mocked fallback is checked via the new
  `bl$predict_method_by_trait` field instead (`all(... == "per_column")`).
- Test 5 (`discrete accuracy under exact is not worse than per_column`):
  pinned, added `predict_method = "exact"` to the "default" call -- this
  test specifically compares the concrete exact route to the concrete
  per_column route (docs/dev-log/exact-default/S5-benchmark-round1.md's
  measured gap), not "auto" to per_column.

Also caught and fixed one over-eager design decision of my own during this
slice (see "Deviations" below): `.fit_baseline_route()`'s uniform-route
case initially always stamped `predict_method_used <- "auto"`, which made
`fit_pigauto(predict_method = "exact")` report `model_config$
predict_method_used == "auto"` after the `baseline_full` reuse fix --
wrong, since the user asked for a concrete route. Fixed via the
`baseline$predict_method_used`-restore relabel described above.

### man/, NEWS.md
- `devtools::document()` regenerated `fit_baseline.Rd`, `fit_pigauto.Rd`,
  `impute.Rd` (the three exported functions whose roxygen changed). The
  one pre-existing roxygen warning (`joint_mvn_solver.R:1002: @param Could
  not resolve link to topic "0, 1"`) reproduced again, confirming S3/S4's
  finding that it predates this lane.
- `NEWS.md`: added a new `## Amendment (S5b): the default is now
  predict_method = "auto", not "exact"` section immediately above the
  existing `## Default flip: predict_method = "exact"` entry (retitled
  that entry `(superseded by "auto" above)` rather than rewriting it, per
  "amend" in the brief) -- states the measured regressions that motivated
  the change (lambda = 1 n >= 300, three small real datasets, AmphiBIO
  habitat), the per-trait loss/tie-break rules, `predict_route` reuse, and
  the measured ~1.4-1.5x cost. No em dashes (checked).

## Route choices and MSEs (DGP c, test 3)

From the test's own `message()` output, one run (seed-controlled,
deterministic): 2 independent continuous traits (true lambda = 1) plus a
near-white-noise count trait, n = 300, `lambda_block` (exact route)
estimated at 0.915 (pulled below 1 by the low-signal count, confirming the
mechanism reported in S5-benchmark-round1.md):

```
c1: chosen=exact       exact=0.1817  per_column=0.2215
c2: chosen=per_column  exact=0.3039  per_column=0.2961
```

`c2` (per_column) satisfies "at least one continuous trait chooses
per_column and its MSE is not worse than exact's" (0.2961 < 0.3039). `c1`
legitimately chose `exact` here (its own val MSE happened to favour it),
which is consistent with the brief's per-trait design -- the two
continuous traits need not agree.

## Timing (item 4)

`fit_baseline()` wall time, `system.time()`, single run each (not
replicated; order of magnitude only):

| fixture | exact | auto | ratio |
|---|---|---|---|
| n = 300, 5 mixed traits (2 continuous, 1 count, 1 binary, 1 categorical) | 0.428 s | 0.607 s | 1.42x |
| n = 2000, 3 continuous traits | 1.936 s | 2.904 s | 1.50x |

Both ratios are well below the naive 2x (fitting the joint baseline twice)
because per-trait dispatch, label propagation, and the shared
phylogenetic-correlation/graph setup are NOT duplicated -- only the joint
solver call itself runs twice.

## Verification

Each of the 28 pre-existing files that `grep -rl "predict_method\|lambda"
tests/testthat/*.R` matches, run individually with `NOT_CRAN=true`: all
green, zero failures (test-active-impute.R, test-community-surface.R,
test-compute-corner-loss.R, test-clamp-outliers.R, test-exact-centring.R,
test-exact-conditional.R, test-exact-default.R [after the 5 pinned-test
fixes above], test-joint-baseline.R, test-joint-lambda.R,
test-honesty-warnings.R, test-joint-threshold-baseline.R,
test-lambda-covariates.R, test-joint-refine-iter.R,
test-lambda-per-type.R, test-lambda-dispatch.R, test-lambda-default.R,
test-monomorphic-discrete.R, test-multi-impute.R, test-pagel-bayes.R,
test-pagel-cv.R, test-pagel-eigendecomp.R, test-pagel-lambda.R,
test-phylo-signal-multiobs.R, test-phylo-signal-gate.R, test-pmm.R,
test-property-invariants.R, test-sigma-fisher-ml.R,
test-zi-count-conformal-mi.R).

`tests/testthat/test-route-choice.R` (new, 7 tests): all green,
individually run.

Manual smoke checks (not part of the automated suite, but exercised while
debugging):
- `fit_baseline(pd, tr, splits)` (auto, real val cells) -> per-trait mix of
  `"exact"`/`"per_column"`, `predict_method_used == "auto"`.
- Forced `predict_route` (all `"exact"`) reproduces an explicit
  `predict_method = "exact"` fit with `max(abs(mu_forced - mu_exact)) == 0`.
- `splits = NULL` -> every trait `"exact"`.
- `impute(..., gnn = FALSE, predict_method = "exact")` (explicit, not
  auto) -> `model_config$predict_method_used == "exact"` (this caught and
  drove the fix to the `.fit_baseline_route()` mislabel described above).

## Outstanding verification (NOT done -- do not claim green from this report)

**Both required commands from the brief were launched but NOT completed by
handback; this report was forced to hand back mid-run.**

1. **Full suite, `NOT_CRAN=true`, `testthat::test_dir()`**: launched in the
   background BEFORE `tests/testthat/test-route-choice.R` existed on disk,
   so that run does NOT include the new file even once it finishes --
   it only re-validates that this slice caused no NEW regressions
   elsewhere. Log: `/tmp/full_suite.log` on the authoring machine (PID
   55399 at handback, ~4 minutes elapsed, still running -- S4's own full
   run took long enough to also be forced to hand back, so this one likely
   needs comparable or longer). Live output at handback showed only
   `Auto-detected trait type` / `Tracing function` informational lines, no
   visible failure markers, but the reporter is `"silent"` so nothing
   would print until the very end regardless -- **this is not evidence of
   FAIL 0**.
2. **`rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"))`**: not
   started yet (was waiting on (1) per the brief's own sequencing).

**Resume commands** (from the worktree):
```r
devtools::load_all(".")
r <- testthat::test_dir("tests/testthat", reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(r)
cat(sprintf("SUITE FAIL %d | PASS %d\n", sum(df$failed), sum(df$passed)))
print(unique(df[df$failed > 0, c("file", "test")]))
```
```r
r <- rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"),
                           error_on = "never", quiet = TRUE)
cat(sprintf("CHECK errors %d warnings %d notes %d\n",
            length(r$errors), length(r$warnings), length(r$notes)))
```
Update this report's "Verification" section with both final counts before
treating S5b as done. Given every individually-run file (28 pre-existing +
the new test-route-choice.R, 29 files total covering every
`predict_method`/`lambda` touch point) is green, a clean full-suite FAIL 0
is the expected outcome, but per this repo's own D-43 discipline that
expectation is not itself a claim of one.

## Deviations from the brief

1. `.fit_baseline_route()`'s uniform-route case initially mislabelled
   `predict_method_used` as `"auto"` even for an explicit
   `predict_method = "exact"`/`"per_column"` request replayed through
   `baseline_full`; caught by a manual smoke test (not the automated
   suite) during this slice and fixed at the `fit_pigauto.R`/`impute.R`
   call sites (restore `baseline$predict_method_used` after the
   `predict_route` replay) rather than inside `fit_baseline()` itself,
   since `fit_baseline()` has no way to distinguish "genuine auto request"
   from "a uniform forced route" from where it's called.
2. `R/multi_impute_trees.R`'s `gnn = TRUE` per-tree replay does NOT get a
   `predict_route` override (left `NULL`, "auto" re-decides per tree) --
   this matches the existing per-tree design (every other baseline option
   is also re-resolved per tree against real `splits_ref` validation
   cells there), and the brief's "no validation cells of its own" argument
   for forcing `predict_route` applies only to the `gnn = FALSE`
   (`splits = NULL`) branch of that same function.
3. Both required verification commands (full suite, `rcmdcheck`) were
   launched but not completed by handback -- see "Outstanding
   verification" above. This is the primary deviation.
4. Lane-check hooks flagged 4-11 other branches per file carrying
   unmerged work on every file touched in this slice (same residual risk
   S3/S4 already flagged for this worktree); not exhaustively diffed
   against every branch.
