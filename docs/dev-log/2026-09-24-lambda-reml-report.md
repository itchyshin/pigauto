# Full-REML fix for the Pagel's lambda profile NLL (2026-09-24)

Branch: `fix/lambda-reml`, worktree `pigauto-lambda-reml`, off `origin/main` at `d86ebc2`.
Owner scope: `R/pagel_lambda.R` and `tests/testthat/{test-pagel-lambda.R,
test-joint-lambda.R, test-lambda-covariates.R, test-pagel-eigendecomp.R,
test-pagel-bayes.R, test-pagel-cv.R}`.

## Problem

`build_pagel_nll_cache()` in `R/pagel_lambda.R` minimised a profile NLL
that omitted the REML determinant correction for profiling out the GLS
mean, `0.5 * log(1' R_oo(lambda)^-1 1)`. This biased the per-column
lambda estimate low at weak-to-moderate signal.

## Change

`R/pagel_lambda.R`:

- Header comment above `build_pagel_nll_cache()` (line ~142-160): rewritten
  to describe the full-REML NLL, including the third determinant term, and
  citing the measured bias reduction.
- Intercept-only branch `nll_at` (line ~211-242): added
  `+ 0.5 * log(sum_a)` (with a `sum_a <= 0` guard mirroring the existing
  `sigma2` guard) to the returned NLL. `sum_a` was already computed for
  `mu_hat`, so no extra work. `mu_at()` (line ~243-250) is untouched --
  it already returns the GLS mean, unaffected by the NLL's REML term.
- X-aware branch `nll_at_x` (line ~258-281): added the matrix analogue,
  `+ 0.5 * log(det(XtDinvX))` via
  `determinant(XtDinvX, logarithm = TRUE)$modulus`, with a `tryCatch`
  guard returning `.Machine$double.xmax` on a singular/error determinant
  (mirrors the existing `solve()` guard on `beta`).
- `bayes_lambda_for_col()` header comment (line ~35-42): added a note that
  the posterior is now a REML posterior (weights `exp(-nll)` unchanged in
  form; `nll` itself changed).
- `lambda = "fixed_1"` never calls `build_pagel_nll_cache()` (confirmed by
  reading `bm_impute_col()` / `.mvn_resolve_lambda()`), so that path is
  bit-identical, as required.

`R/joint_mvn_solver.R::.mvn_resolve_lambda()` (NOT edited, per ownership
boundary): confirmed at line ~377-406 that it builds one
`build_pagel_nll_cache(L[, j], R, ...)` per lambda_cols column and sums
`cache$nll(lam)` over them for the `lambda_block` search, and calls
`.pagel_lambda_from_cache(caches[[idx]])` for the per-trait estimate. Both
consume the shared closure generically, so the REML fix propagates there
with no code change needed.

## Before / after bias (own measurement, real package code, not the
reference scratch script)

Regime matches the PROBLEM statement: 200 seeds, n = 300, 30% MCAR, BM DGP
at Pagel lambda, observed mean shifted +1.5. Measured directly via
`build_pagel_nll_cache()` + `.pagel_lambda_from_cache()` (not a
reimplementation), using `git stash` to get the pre-fix code for the
"before" row and popping it back immediately after:

| lambda_true | before (bias) | after (bias) |
|---|---|---|
| 0.3 | -0.0614 | -0.0219 |
| 0.7 | -0.0383 | -0.0271 |
| 1.0 | -0.0050 | -0.0050 |

These match the PROBLEM statement's numbers (-0.061/-0.038 before,
-0.022/-0.027 after, -0.005 unchanged at lambda = 1) to within Monte
Carlo noise from the seed used for `git stash`/`pop` re-runs (same RNG
stream both times, so numbers here are essentially exact reproductions).

## Consumer check (task step 2)

- `.pagel_lambda_from_cache()`: point-estimate search over `cache$nll` --
  behaviour unchanged in structure, only the surface it searches shifted
  (now a REML surface). No code change needed.
- `bayes_lambda_for_col()`: `exp(-nll(lambda_i))` weights are now REML
  weights, not profile-likelihood weights. Comment updated; no code
  change.
- `.mvn_resolve_lambda()` in `R/joint_mvn_solver.R`: confirmed it sums
  `cache$nll` across `lambda_cols_idx` columns for the block search and
  reuses `.pagel_lambda_from_cache()` per column -- both automatically
  pick up the fix. File not edited (out of ownership scope).

## Tests changed

### `tests/testthat/test-lambda-covariates.R`
- `[lambda-cov] build_pagel_nll_cache(y, R) with X = NULL is unchanged`:
  the captured reference `cache$nll(0.4)` changed from
  `-16.2152757282217` (pre-fix) to `-15.2177731244722` (post-fix, same
  seed/data). Updated the literal and added a comment recording the old
  value.
- `[lambda-cov] lambda = 1 is bit-identical ...` (numeric `lambda = 1.0`,
  never calls `build_pagel_nll_cache`): unaffected, verified by reading
  `bm_impute_col_with_cov()` (`R/bm_internal.R` line ~309-336) -- the
  `"estimate"` branch is the only caller of the cache, and this test uses
  the fixed numeric path.
- `[lambda-cov] recovers lambda with covariates: mean bias at n = 300`
  (600 seeds, bound 0.05, untouched): re-measured bias = -0.0303,
  MAE = 0.1226. Passes with margin (~0.02 headroom). No change needed.

### `tests/testthat/test-pagel-eigendecomp.R`
- `.nll_dense_reference()` (used by two "cached matches dense" tests):
  added the same REML term, `+ 0.5 * log(sum_a)` where
  `sum_a = sum(chol_solve(ones))`, matching the eigenbasis cache's
  `sum_a = sum(c_1^2 * inv_d)` term exactly (same quantity, different
  basis). Tolerances (`1e-8`, `1e-7`) left unchanged; both tests pass.

### `tests/testthat/test-joint-lambda.R`
- `[joint-lambda] recovery: per-trait lambda_hat is not badly biased`
  (40 seeds, was bound 0.08): re-measured with the real
  `fit_mvn_bm_inhouse(lambda = "estimate")` path:
  - seeds 1:40 (the actual production window): bias =
    -0.0272 (lambda=0.3), **-0.0553** (lambda=0.7), -0.0050 (lambda=1)
  - seeds 1:150: bias = -0.0145, -0.0338, -0.0050
  - **Deviation from the task's literal ask**: instructions said "tighten
    to 0.05 if the new estimator meets it with margin." At n = 150 it
    does (max |bias| 0.0338, margin ~0.016), but the actual 40-seed
    window this test runs does NOT -- its lambda=0.7 bias is -0.0553,
    which breaches 0.05 outright, not just "no margin." I additionally
    ran seven further non-overlapping 40-seed windows (seeds 41-320) as a
    sanity check on Monte Carlo noise at n=40; the same lambda=0.7 column
    ranged from -0.0146 to -0.0553 across all eight windows, i.e. the
    production window happens to be the worst one sampled. Tightening
    straight to 0.05 would make this test fail deterministically as
    written. I tightened to **0.07** instead (down from 0.08) -- still a
    real tightening reflecting the improved estimator, but clearing every
    measured 40-seed window (including the production one) with ~0.015
    margin. Comment block above the test rewritten to record all of this.
  - Negative-control test (`lambda = 1` vs `lambda = 0.3` distinguishability,
    threshold 0.05, unrelated to the bias bound) re-run: unaffected, still
    passes (not itself testing bias magnitude).

### `tests/testthat/test-pagel-lambda.R`
- Added `[pagel] full REML term reduces the low bias at weak signal`:
  100 seeds, n = 300, lambda_true = 0.3, 30% MCAR, mean shift +1.5,
  estimator = `pigauto:::ml_lambda_for_col(y, R)`. Assertion:
  `abs(mean(lambda_hat) - 0.3) < 0.04`.
  - With the fix: mean(lambda_hat) = 0.2900, bias = -0.0100. Passes.
  - Negative control: temporarily reverted the one-line REML addition in
    the intercept-only `nll_at` (restored immediately after, verified
    `diff` against a backup showed zero difference from the pre-revert
    state), re-ran this exact test: mean(lambda_hat) = 0.2563,
    bias = -0.0437, which FAILS the 0.04 gate
    (`Expected abs(bias) < 0.04. Actual: 0.0437 >= 0.0400`). Confirms the
    test is discriminating, as required.

### `tests/testthat/test-pagel-bayes.R`, `test-pagel-cv.R`
No hardcoded NLL-derived literals in either file (checked via grep for
`expect_equal` against long decimals and for direct `$nll(` calls); all
assertions are threshold/shape-based (`lambda_post_mean` bounds, weight
sums, `expect_gt`/`expect_lt`). Both files ran unmodified with FAIL 0.

## Verification (task step 5)

Owned six files + the three additional dispatch/default/per-type files,
run individually with `NOT_CRAN=true`:

```
test-pagel-lambda.R:      FAIL 0 | PASS 37
test-joint-lambda.R:      FAIL 0 | PASS 27
test-lambda-covariates.R: FAIL 0 | PASS 18
test-pagel-eigendecomp.R: FAIL 0 | PASS 132
test-pagel-bayes.R:       FAIL 0 | PASS 19
test-pagel-cv.R:          FAIL 0 | PASS 8
test-lambda-dispatch.R:   FAIL 0 | PASS 43
test-lambda-default.R:    FAIL 0 | PASS 13 (1 pre-existing, unrelated WARN:
                           small-validation-set conformal-coverage message)
test-lambda-per-type.R:   FAIL 0 | PASS 20
```

Full suite:

```
NOT_CRAN=true Rscript -e 'devtools::load_all(); r <- testthat::test_dir("tests/testthat", reporter="silent", stop_on_failure=FALSE); df <- as.data.frame(r); cat(sprintf("SUITE FAIL %d | PASS %d\n", sum(df$failed), sum(df$passed)))'
SUITE FAIL 0 | PASS 2599
```

## Deviations from the literal task

1. `test-joint-lambda.R` bias bound tightened to **0.07**, not 0.05 (see
   above -- 0.05 does not clear the actual 40-seed production window with
   margin; it fails outright at -0.0553). Reported both the 40-seed and
   150-seed measurements as instructed.
2. No other deviations. `NEWS.md`, `man/`, `R/multi_impute.R`,
   `R/henderson_s_inv.R` untouched. No commits made, nothing pushed.

## Files touched

- `R/pagel_lambda.R`
- `tests/testthat/test-pagel-lambda.R`
- `tests/testthat/test-joint-lambda.R`
- `tests/testthat/test-lambda-covariates.R`
- `tests/testthat/test-pagel-eigendecomp.R`

(`test-pagel-bayes.R`, `test-pagel-cv.R` read and verified, not modified.)

## Note on a lane-check hook warning

The pre/post-edit hooks flagged that branch `chore/dedupe-agent-files`
also carries commits touching `R/pagel_lambda.R` and
`tests/testthat/test-pagel-lambda.R`. Inspected via
`git diff HEAD..chore/dedupe-agent-files -- R/pagel_lambda.R`: that
branch's diff *removes* the X-aware branch and `mu_at` entirely (appears
to be based on an older pre-X-aware commit, not a REML fix) -- it is not
a competing implementation of this fix. Not merged with or built on top
of; flagging here per the task's instruction to report deviations/context
rather than silently acting on another lane's branch.
