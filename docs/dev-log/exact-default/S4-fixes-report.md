# S4 fixes report -- clearing the failures left after the S3 default flip

Branch: feat/exact-default (worktree /Users/z3437171/Dropbox/Github
Local/pigauto-exact-default). Not committed, not pushed, per instructions.

STATUS: DONE. All six named test files individually re-run FAIL 0, a
`NOT_CRAN=true` full-suite run produced **SUITE FAIL 0 | PASS 2637**, and
`rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"))` completed with
**CHECK errors 0 warnings 0 notes 1** (the single NOTE is
"checking CRAN incoming feasibility" flagging the maintainer email and
that the version string "contains large components" (0.11.0.9000) --
routine for a dev version, not caused by this slice).

`devtools::document()` was run: no `.Rd` changes resulted (no roxygen tags
were touched in this slice, only code comments and test files). The one
pre-existing roxygen warning
(`joint_mvn_solver.R:1002: @param Could not resolve link to topic "0, 1"`)
reproduced again, confirming S3's finding that it is not caused by this
lane.

## Root cause A (code fix): lambda_fixed rebuild did not reproduce the fit

**Problem.** Under `predict_method = "exact"`, `fit_baseline(...,
lambda_fixed = bl$lambda_per_trait)` did not reproduce `bl$mu` / `bl$se`.
Two independent bugs combined to cause this:

1. `.mvn_resolve_lambda()`'s numeric-vector branch (a full per-column
   `lambda` vector supplied by the caller) always set
   `lambda_block <- mean(lambda_vec)`. But the ORIGINAL "estimate" fit's
   `lambda_block` came from the argmin of the SUMMED profile-REML NLL over
   `lambda_cols`, which is not generally the mean of the resulting
   per-column `lambda_hat` values -- and, under `predict_method =
   "exact"`, the non-`lambda_cols` (discrete) columns are already fixed AT
   `lambda_block` by the S3 exception, which biases `mean(lambda_vec)`
   toward `lambda_block` without reproducing it exactly. Since
   `lambda_block` drives the shared `henderson_bar`, the Sigma M-step, and
   the exact conditional's per-column GLS-mean centring, any drift in it
   changes `mu`/`se` at every imputed cell.

2. `fit_joint_threshold_baseline()`'s own `lambda_fixed` handling
   (R/joint_threshold_baseline.R:352-363, before this fix) only read
   continuous-family (`lambda_family_idx`) columns from `lambda_fixed`;
   every other column (the discrete liability columns fit in the SAME
   threshold-joint call, e.g. binary) was hardcoded to 1 regardless of
   what `lambda_fixed` actually carried for that column. Under the
   original "estimate" + "exact" fit, those discrete columns are set to
   `lambda_block` (not 1) by `.mvn_resolve_lambda()`'s S3 exception, so a
   rebuild that forces them back to 1 changes the per-column-BM init
   (`.mvn_init_per_column()`) that feeds Sigma, again moving `mu`/`se`
   away from the original fit.

**Fix (file:line).**
- `R/joint_mvn_solver.R`, `.mvn_resolve_lambda()`'s numeric-vector branch
  (around what was lines 382-397, now expanded): read
  `attr(lambda, "lambda_block")` and use it as `lambda_block` when finite;
  fall back to `mean(lambda_vec)` only when the attribute is absent.
- `R/fit_baseline.R`, just before building the `out <- list(...)` return
  value (~line 1127): `attr(lambda_per_trait, "lambda_block") <-
  lambda_block_out` when `lambda_block_out` is finite, so a caller that
  later replays `bl$lambda_per_trait` as a fresh `lambda_fixed` argument
  carries the exact block value forward.
- `R/joint_mvn_baseline.R` (~line 132-144) and
  `R/joint_threshold_baseline.R` (~line 352-373): capture
  `attr(lambda_fixed, "lambda_block")` BEFORE subsetting by column name
  (subsetting a numeric vector drops custom attributes, keeping only
  `names`), then reattach it to the resolved `lam_arg`/`la` vector after
  subsetting, so the attribute survives the hop into
  `fit_joint_solver()` / `fit_mvn_bm_inhouse()`.
- `R/joint_threshold_baseline.R` (~line 352-373): additionally changed the
  `lambda_fixed` lookup itself to read EVERY `X_fit` column by name
  (`unname(lambda_fixed[colnames(X_fit)])`, defaulting missing names to
  1) instead of restricting the lookup to `lambda_family_idx`. This is
  the general fix for bug 2 above: a discrete column simply absent from
  `lambda_fixed` still defaults to 1 (the pre-existing, still-correct
  per-column default), but a discrete column that IS present in
  `lambda_fixed` (as happens when replaying a previous "estimate" +
  "exact" fit's own `$lambda_per_trait`, which now legitimately carries
  `lambda_block` for those columns) is honoured rather than silently
  overridden.

After this fix, `fit_baseline(..., lambda_fixed = bl$lambda_per_trait)`
reproduces `bl$mu` / `bl$se` to `tolerance = 1e-8` under
`predict_method = "exact"` -- confirmed by the now-passing
"[lambda-dispatch] lambda_fixed rebuild" and "[lambda-default] predict
rebuild" tests below, both left UNPINNED (general property).

## Per-test disposition

### test-exact-conditional.R
- **"[exact] fit_baseline default is unchanged and 'exact' differs"**:
  already fixed by the previous builder (S3's own update, per its report)
  before this slice started. Re-ran clean, no change needed.

### test-joint-lambda.R
- All three named tests ("fixed_1 is identical...", "mean model:...",
  "block lambda: a column excluded...") were already pinned to
  `predict_method = "per_column"` by the previous builder before this
  slice started (per S3's UPDATE section). Re-ran clean, no change
  needed.

### test-joint-refine-iter.R
- "joint_refine_iter = 3L runs ... differs from 0L": already pinned by
  the previous builder before this slice started. Re-ran clean (1
  pre-existing warning, confirmed in S1, unrelated to this slice), no
  change needed.

### test-lambda-default.R
- **"predict rebuild: predict() reproduces the fit-time baseline at
  estimated lambda"**: CODE FIX (root cause A above). Left UNPINNED
  (general property: a baseline rebuilt from its own reported lambdas
  must reproduce the fit). Now passes under the exact default.
- **"estimated lambda is applied"**: PINNED
  (`predict_method = "per_column"` added to both `impute()` calls, plus a
  comment). Part (b) of this test compares the estimate-mode prediction
  directly against an INDEPENDENT `bm_impute_col()` call -- a per-column
  equivalence that does not hold under `"exact"`, where the joint
  conditional predicts every column together from one shared Sigma /
  R(lambda_block). Part (a) (estimate vs fixed_1 differ) still holds
  either way; pinning both calls keeps the whole test internally
  consistent.

### test-lambda-dispatch.R
- **"lambda_per_trait populated"**: CODE FIX (per the coordinator's
  instruction A) -- rewrote the assertion rather than the underlying
  behaviour. Split `disc_names` into `bin_names` (the binary column `b1`,
  fit IN THE SAME threshold-joint call as the continuous columns, and
  therefore subject to the S3 exact-route exception: its
  `lambda_per_trait` entry now must equal `bl$lambda_block` exactly, not
  merely "not 1") and `cat_names` (the categorical column `k1`'s OVR
  latent columns, which R/ovr_categorical.R always fits at
  `lambda_mode = "fixed_1"` regardless of the caller's `predict_method` --
  confirmed by direct measurement: these stayed at 1 even under exact).
  Also added a `predict_method = "per_column"` companion assertion that
  every discrete column (binary AND categorical) stays at 1 on that
  route, unchanged from pre-S3 behaviour.
- **"lambda_fixed rebuild"**: CODE FIX (root cause A above). Left
  UNPINNED. Now passes under the exact default.
- **"partial lambda_fixed defaults missing columns to 1 (threshold-joint)"**
  and **"... (joint MVN)"**: PINNED (`predict_method = "per_column"`
  added to every `fit_baseline()` call in both tests), with the comment
  "per-column contract; under the exact default this does not hold
  because the exact conditional couples every joint column through one
  shared Sigma / R(lambda_block), so fixing ONE column's lambda still
  changes the OTHER columns' predictions." Added two new exact-route
  counterpart tests ("... runs under exact and uses the named lambda
  (threshold-joint)" / "(joint MVN)") asserting the partial call runs
  without error under the default and that the NAMED column
  (`lambda_per_trait["c1"]`) is fit at exactly the given lambda, per the
  coordinator's instruction.
- **"fixed_1 dispatcher reference"**: PINNED
  (`predict_method = "per_column"` added), with a comment: the reference
  fixture (`lambda_fixed1_reference_ab02e31.rds`) was captured before
  `predict_method = "exact"` existed as a default, and the exact
  conditional legitimately produces different imputed-cell numbers even
  at `lambda_mode = "fixed_1"` (lambda = 1 is still the LAMBDA used, but
  "exact" vs "per_column" is a different PREDICTION step entirely).
- **"ordinal stays at lambda = 1 under estimate"**: PINNED
  (`predict_method = "per_column"` added to both `fit_baseline()` calls),
  with a comment explaining the S3 exception (ordinal is fit inside the
  same threshold-joint call as continuous columns, via B3's
  `estep_liability_ordinal`, so it is subject to the same exact-route
  lambda_block exception as binary). Added a new exact-route counterpart
  test ("ordinal mu is finite and reports lambda_block under exact
  estimate") asserting the ordinal column's mu is finite and
  `lambda_per_trait` equals `lambda_block` (when finite) under the
  default, across 3 seeds.

### test-sigma-fisher-ml.R
- **"sigma_method = 'single_pass' numeric regression on a fixed-seed
  L/tree"**: PINNED (`predict_method = "per_column"` added to the single
  `fit_mvn_bm_inhouse()` call), with a comment: the pinned
  `sum(fit$anc_recon)` / `sum(fit$anc_var)` values were computed via the
  per-column prediction; `predict_method = "exact"` (now the default)
  returns the exact matrix-normal conditional mean/variance instead.
  `fit$pars$phylocov` (Sigma) is unaffected by `predict_method` and was
  left unpinned/unchanged -- confirmed still passing at its original
  tolerance.

## Files changed in this slice

- `R/joint_mvn_solver.R` -- `.mvn_resolve_lambda()` numeric-vector branch
  (lambda_block attribute read).
- `R/fit_baseline.R` -- attach `lambda_block` attribute to the returned
  `$lambda_per_trait`.
- `R/joint_mvn_baseline.R` -- carry the attribute through the
  `lambda_fixed` name-subset.
- `R/joint_threshold_baseline.R` -- carry the attribute through the
  `lambda_fixed` name-subset; removed the `lambda_family_idx` restriction
  so every `X_fit` column (not just continuous-family ones) is read from
  `lambda_fixed` by name, defaulting missing names to 1.
- `tests/testthat/test-lambda-dispatch.R` -- pinned 4 tests, rewrote 1
  assertion, added 3 new exact-route counterpart tests.
- `tests/testthat/test-lambda-default.R` -- pinned 1 test
  (`predict_method = "per_column"` on both `impute()` calls).
- `tests/testthat/test-sigma-fisher-ml.R` -- pinned 1 test.

No deviations from the coordinator's brief. No failure matched category C
(all fell under A or B as described in the brief).

## Verification

Each of the six named files individually, `NOT_CRAN=true`:

```
test-exact-conditional.R  FAIL 0
test-joint-lambda.R       FAIL 0
test-joint-refine-iter.R  FAIL 0
test-lambda-default.R     FAIL 0
test-lambda-dispatch.R    FAIL 0
test-sigma-fisher-ml.R    FAIL 0
```

Full suite (`NOT_CRAN=true`):

```
SUITE FAIL 0 | PASS 2637
```

`rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"))`:

```
CHECK errors 0 warnings 0 notes 1
[1] "checking CRAN incoming feasibility ... [4s/24s] NOTE
Maintainer: 'Shinichi Nakagawa <itchyshin@gmail.com>'

Version contains large components (0.11.0.9000)"
```

This NOTE is routine for a dev-version maintainer/incoming-feasibility
check and is not caused by this slice.
