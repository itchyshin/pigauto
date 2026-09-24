# S5: lambda_mode default flip and plumbing

Lane: feat/joint-lambda-default. Files owned/touched: `R/fit_pigauto.R`, `R/impute.R`,
`R/multi_impute.R`, `R/multi_impute_trees.R`, `NEWS.md`, `DESCRIPTION`, `man/` (regenerated),
`tests/testthat/test-lambda-default.R` (new). Did NOT touch `R/fit_baseline.R`,
`R/joint_mvn_solver.R`, `R/joint_*_baseline.R`, `R/ovr_categorical.R`, or `R/predict_pigauto.R`
(explained below).

## Important finding: S4 landed in this shared worktree while this task ran

This is a shared worktree with other concurrent agents. Partway through this task, S4's edits to
`R/fit_baseline.R`, `R/joint_mvn_solver.R`, `R/joint_mvn_baseline.R`, `R/joint_threshold_baseline.R`
appeared as uncommitted modifications (confirmed via `git status`, not something I touched). This
matters for two reasons:

1. My initially-written tests used `skip_if(is.null(bl$lambda_per_trait), "S4 not landed")` guards on
   the assumption those fields would be `NULL`. Once S4 landed, `fit_baseline()` actually returns
   `lambda_per_trait`, `lambda_block`, `lambda_mode`, and accepts `lambda_fixed` -- so those guards
   stopped triggering and the tests now exercise the REAL contract (all pass; see Test summary).
2. **The default-flip's actual mechanism changed underneath me.** I originally measured (pre-S4)
   that `lambda_mode = "estimate"` on 2+ continuous-family columns forced the per-column BM path
   (`force_per_column <- lambda_mode %in% c("estimate", "cv", "bayes")`). After S4 landed,
   `force_per_column <- lambda_mode %in% c("cv", "bayes")` -- "estimate" no longer forces per-column;
   the joint MVN baseline now estimates lambda internally and stays active
   (`bl$path` reports `"joint_mvn"` on both `"fixed_1"` and `"estimate"` runs). I re-measured after
   this landed and corrected `NEWS.md` accordingly (see below) -- the numeric difference is still
   ~0.38 max abs diff on the smoke fixture, but the reason is different, and I do not want to
   misattribute S4's mechanism in `NEWS.md` (S4's own report, not yet written as of this report, is
   the authoritative source for its internals).

## What changed (task by task)

### 1. Default flip, enum kept, per the contract

- `R/impute.R`: `lambda_mode = c("fixed_1", "estimate", "cv", "bayes")` ->
  `c("estimate", "fixed_1", "cv", "bayes")` (line, now ~341). Roxygen updated.
- `R/fit_pigauto.R`: same reorder (line ~339). Roxygen updated with the `predict_method = "exact"` /
  `joint_refine_iter > 0` -> `lambda_block` note.
- `R/multi_impute.R`: had NO `lambda_mode` argument before this change (only reachable via `...`
  forwarded to `impute()`). Added an explicit
  `lambda_mode = c("estimate", "fixed_1", "cv", "bayes")` argument, `match.arg()`'d, and threaded
  into BOTH `impute()` call sites (`mc_dropout` and `conformal` branches).
- `R/multi_impute_trees.R:562` (line moved from 559 after edits):
  `baseline_arg("lambda_mode", "fixed_1")` -> `baseline_arg("lambda_mode", "estimate")`. This is
  only the FALLBACK used when neither `dots$lambda_mode` nor `res_ref$fit$model_config$lambda_mode`
  is set; in practice `model_config$lambda_mode` is always set by `fit_pigauto()` now, so this
  fallback is a defensive no-op in the common case, but it is what the "per tree" test exercises when
  the config path is followed (it now defaults through to "estimate" either way, since
  `impute()`'s own default flipped).

Gate check (verbatim):
```
> Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); cat(eval(formals(impute)$lambda_mode)[1], eval(formals(fit_pigauto)$lambda_mode)[1], eval(formals(multi_impute)$lambda_mode)[1], "\n")'
estimate estimate estimate
```
`grep -rn 'c("fixed_1", "estimate"' R/` finds only `R/fit_baseline.R` (explicitly out of scope --
"fit_baseline's default flips in S4" per the contract; confirmed untouched by me).

### 2. `model_config$lambda_per_trait` / `$lambda_block`

Added to BOTH `model_config` construction sites in `R/fit_pigauto.R`:
- Main (gnn = TRUE) path (~line 1452): `lambda_per_trait = baseline$lambda_per_trait %||% NULL`,
  `lambda_block = baseline$lambda_block %||% NULL` -- sourced from `baseline` (not `baseline_full`)
  because `predict.pigauto_fit()` documents that GNN-on fits always use `object$baseline`.
- `gnn = FALSE` early-return path (~line 694): sourced from `baseline_full` with fallback to
  `baseline`, matching `predict.pigauto_fit()`'s documented `baseline_full`-for-production /
  `baseline`-for-evaluation split.

Verified populated post-S4-landing:
```
res$fit$model_config$lambda_per_trait  ->  x1: 0.005, x2: 0.0113
res$fit$model_config$lambda_block      ->  0.0101
```

### 3. `predict.pigauto_fit`: carried, not recomputed

`grep -n "model_config\$lambda\|lambda_mode\|lambda_per_trait\|lambda_block\|lambda_fixed"
R/predict_pigauto.R` returns nothing. Confirmed by reading the function: it selects between
`object$baseline` and `object$baseline_full` (or a `baseline_override`) and uses their `mu`/`se`
matrices directly -- there is no `fit_baseline()` call inside `predict.pigauto_fit()` at all. The
baseline (and hence its lambda) is entirely CARRIED on the fit object, never recomputed at predict
time. No code change was needed in `R/predict_pigauto.R` for this task; per the contract's
conditional instruction I document this rather than add a rebuild call that doesn't exist yet.

### 4. Roxygen + NEWS + DESCRIPTION

- `impute()`, `fit_pigauto()`, `multi_impute()` roxygen `@param lambda_mode` updated: default
  `"estimate"`, discrete traits stay at lambda = 1, and the `predict_method = "exact"` /
  `joint_refine_iter > 0` -> `lambda_block` interaction. `devtools::document()` regenerated
  `man/fit_pigauto.Rd`, `man/impute.Rd`, `man/multi_impute.Rd` only (verified via
  `git status --short man/`).
- `NEWS.md`: new `# pigauto 0.11.0.9000 (dev)` header above the existing `# pigauto 0.11.0` entry,
  documenting the default flip, the real (not just cosmetic) numeric effect, and the new
  `model_config` / `fit_baseline()` fields. Revised once during this task after the S4-landing
  finding above (see "Important finding").
- `DESCRIPTION`: `Version: 0.11.0` -> `0.11.0.9000`.

### 5. `tests/testthat/test-lambda-default.R` (new)

Five `test_that()` blocks, all titled `"[lambda-default] ..."`:
- (i) `formals()` defaults for `impute`/`fit_pigauto`/`multi_impute` all evaluate to `"estimate"`.
- (ii) `impute(df, tree, gnn = FALSE)` stores `model_config$lambda_mode == "estimate"` and a numeric
  `lambda_per_trait` of length `ncol(X_scaled)` (guarded with `skip_if(is.null(...), "S4 not
  landed")`; now passing for real).
- (iii) "predict rebuild" test: `fit_baseline(..., lambda_fixed = bl$lambda_per_trait)` reproduces
  `bl$mu`/`bl$se` to 1e-8 (guarded the same way; now passing for real).
- (iv) "per tree" test: `multi_impute_trees()` on `tree1` (`ape::rcoal`) and
  `transform_tree_pagel(tree1, 0.2)` stores DIFFERENT `lambda_per_trait` per tree. Required a
  non-boundary DGP: a pure-iid or pure-BM(lambda=1) trait pushes the ML estimate to 0.005 or 0.995
  on BOTH trees regardless of tree identity (checked explicitly), which defeats the test's purpose.
  Used `lambda_true = 0.5` mixed phylo/iid signal (BM simulated under `tree1`'s own correlation
  matrix), which lands in the interior for `tree1` (x1: 0.235, x2: 0.403) and is pulled toward the
  boundary for `tree2` (x1: 0.78, x2: 0.995) -- clearly different. Required a small addition to
  `R/multi_impute_trees.R`'s `run_shared_gnn()`: a new `lambda_per_trait_by_tree` list accumulator
  (one entry per tree, captured right after each tree's `fit_baseline()` call) added to the returned
  `pigauto_mi_trees` object -- this field did not exist before and nothing else currently reads or
  writes it, so this is additive only. Only added to `run_shared_gnn()` (the `share_gnn = TRUE`
  default path), which fits ONE shared GNN with per-tree baseline-only refits, so per-tree lambda
  was otherwise unobservable from the returned object. `run_per_tree()` (`share_gnn = FALSE`) needed
  no change: it already returns `fits` (one full `pigauto_fit` per tree), and each fit's own
  `model_config$lambda_per_trait` (task 2's addition) already carries this information per tree.
- (v) `NEWS.md` mentions "lambda" in its first 60 lines (`readLines(..., n = 60)` +
  `grepl(..., ignore.case = TRUE)`).

## Test summary (verbatim)

```
> Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); r <- testthat::test_file("tests/testthat/test-lambda-default.R", reporter="summary"); df <- as.data.frame(r); cat(sprintf("FAIL %d | PASS %d\n", sum(df$failed), sum(df$passed)))'
lambda-default: ...........
FAIL 0 | PASS 11
```
All 11 expectations pass; 0 skipped (S4 landed during this task, so the "S4 not landed" guards do
not trigger -- see "Important finding" above).

Required smoke files (NOT_CRAN=true):
```
test-lambda-default.R      FAIL 0 | PASS 11
test-multi-impute.R        FAIL 0 | PASS 202 | WARN 32
test-fit-predict.R         FAIL 0 | PASS 87  | WARN 18
test-multi-impute-trees.R  FAIL 0 | PASS 38  | WARN 29
```
`test-multi-impute-trees.R`'s 29 warnings are all pre-existing "Small validation set" /
few-tree-sensitivity-draws messages from its own tiny test fixtures (2 trees, epochs = 5),
unrelated to this change. The run also printed two unrelated stderr lines from an unloaded
`gllvmTMB` finalizer (`"FreeADFunObject" not available`) during garbage collection between test
files -- a different package's teardown noise in this shared R session, not a pigauto failure
(exit code 0, FAIL 0 throughout).

Also checked (not explicitly required, but touched by this task's `R/multi_impute_trees.R` edit and
`R/fit_pigauto.R`'s `model_config` edits): `test-gnn-off.R` (NOT_CRAN=true) FAIL 0 | PASS 73 | WARN 14
(warnings are pre-existing small-validation-set messages, unrelated to this change);
`test-share-gnn.R` (NOT_CRAN=true) FAIL 0 | PASS 35 | WARN 1.

## Deviations / notes

1. **NEWS.md was rewritten mid-task** after discovering S4's actual dispatch mechanism differs from
   what I could observe before S4 landed (see "Important finding"). The numeric claim (~0.38 max abs
   diff) is unchanged and re-verified after S4 landed; the mechanism description was corrected from
   "forces the per-column path" to "the joint solver estimates lambda internally, still dispatching
   to joint_mvn" to avoid misattributing S4's design.
2. **`R/multi_impute_trees.R` gained one new field** (`lambda_per_trait_by_tree` on the
   `pigauto_mi_trees` object) beyond the literal line-559 fallback-string change the contract named.
   This was necessary to make the "per tree" test meaningful at all -- without it there was no way
   for a caller (or a test) to observe the per-tree lambda estimates that `run_shared_gnn()` already
   computes internally and discards. It is purely additive (a new list field, `NULL` per entry until
   `fit_baseline()` populates `lambda_per_trait`, which it now does).
3. No changes to `R/predict_pigauto.R`: the baseline is carried on the fit, never recomputed at
   predict time, so there was nothing to thread `lambda_mode` through there. Documented in task 3's
   report above rather than adding unrequested rebuild logic.
