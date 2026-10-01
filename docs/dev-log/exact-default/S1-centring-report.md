# S1 centring report -- predict_method = "exact" mean-model fix

Branch: feat/exact-default (off origin/main eccb298). Files touched:
R/joint_mvn_solver.R, tests/testthat/test-exact-centring.R (new).
R/exact_conditional.R was read but not changed -- the centring lives
entirely in the caller, not in exact_conditional_mvn() itself.

## What changed

`R/joint_mvn_solver.R`, `fit_mvn_bm_inhouse()`:

- New argument `exact_centre = TRUE` (signature line ~536, after
  `lambda_cols = NULL`).
- Roxygen-style comment block above the function (~509-527) documents
  the argument and its rationale.
- The exact-conditional branch (~684-720, was ~664-684 before this
  edit): for each column k, when `exact_centre` is TRUE, computes
  `mu_exact[k] <- .mvn_gls_mean_at_lambda(L[, k], R, lambda_block,
  nugget = eps)`, subtracts it from `L[, k]` to form `L_exact`, calls
  `exact_conditional_mvn(L_exact, Sigma, henderson_bar, ...)` on the
  centred matrix, then adds `mu_exact` back to every cell of the
  returned mean via `sweep(ec$mu, 2L, mu_exact, \`+\`)`. Variance
  (`ec$var`) is returned unchanged and untouched by centring.
- `fit_joint_solver()` (~904-953): new `exact_centre = TRUE` argument,
  forwarded to both `fit_mvn_bm_inhouse()` call sites (the
  `joint_solver = "inhouse"` path and the rphylopars-failure fallback
  path). Not consulted on the `"rphylopars"` success path. New
  `@param exact_centre` roxygen entry; function stays `@noRd`, so no
  man/ regeneration was needed or done.

No change to `R/exact_conditional.R`: the centring is entirely a
caller-side transform (subtract before the call, add back after), so
the solver itself, its tests, and its contract are untouched.

## Item 1: does centring at lambda_block = 1 change existing exact output?

Yes, measurably. Reference fixture: n = 20, K = 3, `ape::rcoal(20)`
(seed 42), `Sigma = [[1,.5,-.2],[.5,1,.3],[-.2,.3,1]]`, 30% MCAR
missingness (seed 43), default `lambda = "fixed_1"` so
`lambda_block = 1`. Measured GLS means at lambda = 1 on the observed
subset of each column: column a = -0.630, column b = -1.312,
column c = +0.654 -- all far from the model's assumed zero. Running
`fit_mvn_bm_inhouse(..., predict_method = "exact")` before and after
the edit (`exact_centre = FALSE` vs `TRUE`) on this fixture:

- `exact_centre = FALSE` vs the pre-edit (unmodified) function: max
  abs diff = 0 (exact reproduction, as designed).
- `exact_centre = TRUE` vs the pre-edit function: max abs diff =
  1.3346.
- `anc_var` is identical between `exact_centre = TRUE` and `FALSE`
  (`all.equal` TRUE), confirming the variance-unchanged design.

Decision per the task brief: apply centring for all lambda_block
values, including 1, because it is the correct model (the exact
solver's zero-mean assumption is wrong whenever the data's own mean
isn't zero, independent of lambda) -- `exact_centre` defaults to TRUE.

### Which existing tests move

None of the existing numeric-value tests move, because none of them
exercise `fit_mvn_bm_inhouse(predict_method = "exact")`'s literal
output values with a pinned expectation:

- `test-exact-conditional.R`'s numeric checks all call
  `exact_conditional_mvn()` directly (bypassing `fit_mvn_bm_inhouse()`
  and hence this centring step entirely).
- `test-exact-conditional.R`'s two tests that DO go through
  `fit_mvn_bm_inhouse()` / `fit_baseline()` / `impute()` only assert
  finiteness, dimension match, and "differs from per_column" /
  "differs from per_column" -- none pin exact values, so they still
  pass with the new default.

Confirmed by running the full suite below: FAIL 0 everywhere.

## Item 3b: RMSE, centred vs uncentred, mean-shifted DGP

DGP (see `.centring_dgp()` in the new test file): n = 200, K = 3
correlated traits (Sigma as above), `lambda_true = 0.3` (chosen
because Rose's review previously measured up to 0.25 mean-model bias
at this same lambda on the per-column path), every column shifted by
+1.5 after drawing `L ~ MVN(0, Sigma %x% R(lambda))`, MCAR missingness
with a floor of >= 5 observed cells per column, 20 seeds (1:20).

Missingness fraction is 0.85, not a lower value -- measured while
designing the test: at miss = 0.3, mean RMSE was 0.7921 (centred) vs
0.7938 (uncentred), a 10/20-seed coin flip; at miss = 0.5, 14/20 wins;
at miss = 0.85, 16/20 wins with a clearer mean gap. The mechanism: at
low missingness the UNCENTRED exact conditional partially self-corrects
the mean by averaging over many weakly-correlated observed cells (a
law-of-large-numbers effect); high missingness removes that crutch,
which is exactly the regime this fix targets.

Reported numbers (miss = 0.85, 20 seeds):

- Mean RMSE, centred: 0.8930
- Mean RMSE, uncentred: 0.9290
- Seeds where centred RMSE < uncentred RMSE: 16/20

Test asserts `mean(rmse_centred) < mean(rmse_uncentred)`; passes.

## Item 3c: zero-mean DGP, centred vs uncentred agreement tolerance

Same DGP with `shift = 0`, `lambda_true = 0.3`, n = 200, K = 3, miss =
0.3 (lower missingness here, so the GLS mean estimate itself
concentrates closer to zero under the true-zero DGP), 20 seeds.
Measured `max(abs(centred - uncentred))` per seed: max over seeds =
0.283, mean over seeds = 0.087. This nonzero gap is sampling noise in
the realised GLS mean estimate (finite n_obs even under a genuinely
zero-mean DGP), not a design flaw. Test asserts
`max(max_diff across seeds) < 0.35` (measured max 0.283 plus margin);
passes.

## Item 3d: fallback path

`fit_mvn_bm_inhouse()` does not expose `exact_conditional_mvn()`'s own
`max_cells` argument -- that function's own "refuses oversized
problems" behaviour is already covered directly in
`test-exact-conditional.R`. What needed protecting here is the
SURROUNDING fallback logic in `fit_mvn_bm_inhouse()` (`ec <-
exact_conditional_mvn(...); if NULL, warn and fall through to
per_column`), now that a centring step runs before that call. Tested
via `testthat::local_mocked_bindings(exact_conditional_mvn = function(...)
NULL)`, which reproduces exactly the oversized-problem return value.
Confirms: a warning fires ("falling back to the per-column path"), and
the fallback result (`anc_recon`, `anc_var`) is identical to calling
`fit_mvn_bm_inhouse(..., predict_method = "per_column")` directly, and
finite. Passes. Deviation from the literal instruction ("set max_cells
tiny") noted above and in the reply.

## Test summary lines (verbatim, NOT_CRAN=true)

```
test-exact-centring.R               FAIL 0 | WARN 0 | SKIP 0 | PASS 12
test-exact-conditional.R            FAIL 0 | WARN 0 | SKIP 0 | PASS 22
test-joint-lambda.R                 FAIL 0 | WARN 0 | SKIP 0 | PASS 27
test-joint-solver.R                 FAIL 0 | WARN 1 | SKIP 0 | PASS 14
test-lambda-dispatch.R              FAIL 0 | WARN 0 | SKIP 0 | PASS 43
test-lambda-per-type.R              FAIL 0 | WARN 0 | SKIP 0 | PASS 20
```

The WARN 1 in test-joint-solver.R is pre-existing: confirmed by
stashing this slice's changes and re-running the same file against
unmodified origin/main HEAD (eccb298) -- same FAIL 0 / WARN 1 / PASS
14. Not caused by this slice.

As predicted, `test-lambda-per-type.R` is unaffected: it exercises the
default `predict_method = "per_column"` path throughout, which this
slice never touches.

Also ran (not on the required list, but they call
`fit_mvn_bm_inhouse()`/`fit_joint_solver()` and were checked as a
walk-around per the Rose principle): `test-joint-baseline.R` (FAIL 0 /
WARN 0 / PASS 26), `test-joint-refine-iter.R` (FAIL 0 / WARN 1 / PASS
12, same pre-existing warning pattern), `test-sigma-fisher-ml.R`
(FAIL 0 / WARN 0 / PASS 19). None of these three files pass
`predict_method` at all (grep confirmed), so they exercise
`predict_method = "per_column"` only and are structurally unreachable
by this change.

## Deviations from the brief

1. Item 3d ("set max_cells tiny"): `fit_mvn_bm_inhouse()` has no
   `max_cells` parameter to set -- it isn't threaded through from
   `exact_conditional_mvn()`. Used `local_mocked_bindings()` to force
   the same NULL-return / fallback code path instead of driving a
   genuinely oversized problem. Scope of the fix (add `exact_centre`
   only) did not include adding a new `max_cells` passthrough
   parameter, which was not requested.
2. Item 3b/3c missingness fractions (0.85 and 0.3 respectively) were
   chosen by measurement, not handed a fixed value in the brief beyond
   n/K/lambda/shift/seeds -- reported and justified above.
3. Pre-existing lane-check hooks flagged `origin/codex/active-recovery-evidence`
   and `origin/ci-run/25974486015` as carrying unmerged work on
   R/joint_mvn_solver.R. Checked both: `active-recovery-evidence`'s
   diff on this file is centered on removing `.mvn_gls_mean_at_lambda`
   entirely relative to a much older merge-base (7d775c0) predating the
   lambda work in this checkout -- a stale/archived branch, not active
   competing work on `exact_centre`. `ci-run/25974486015` is an
   archived CI run branch. Neither touches the exact-conditional
   centring area; no fork of in-flight work occurred.
