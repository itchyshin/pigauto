# Fixes for Rose's final review (PR #187, feat/joint-lambda-default)

Reviewer input: `docs/dev-log/lambda-default/rose-final-review.md` (2026-09-23). This report
covers the 6 fix items assigned, mapped to Rose's "Required changes" 1, 2, 3, 4 (partial), 5, 6, 7, 8.

## 1. Ordinal leak

**What changed.** The per-trait ordinal path-selection block in `fit_baseline()` computed its
"BM-via-MVN alternative" candidate with `lambda = bm_lambda` (the CONTINUOUS-FAMILY setting, which
tracks `lambda_mode`), instead of the documented `lambda = 1` contract for ordinal columns. Changed
the literal argument to `lambda = 1.0`.

- File:line: `R/fit_baseline.R:626` (was using `bm_lambda`, now hardcoded `1.0`, with a comment
  explaining why `bm_lambda` must not be substituted here).
- Test: `tests/testthat/test-lambda-dispatch.R`, `"[lambda-dispatch] ordinal stays at lambda = 1
  under estimate"` -- 2 continuous BM(lambda = 0.3) traits + 1 three-level ordinal trait, 6 seeds
  (900:905), `expect_identical()` (tolerance 0) on `mu`/`se` for the ordinal column between
  `fixed_1` and `estimate`, plus `lambda_per_trait["o1"] == 1`.
- Mutant check: reverted the fix (`lambda = bm_lambda`) and re-ran the test -- it failed on 1 of the
  6 seeds (differences up to ~0.9 on `se`, confirming the test is load-bearing, not vacuous). Restored
  the fix and re-ran clean.

## 2. Partial `lambda_fixed`

**What changed.** `unname(lambda_fixed[colnames(...)])` produced `NA` for any column absent from
`lambda_fixed`, which the downstream solver's numeric-vector validation then rejected
(`"numeric 'lambda' vector must have all entries in [0, 1]"`), contradicting the documented "columns
not present in `lambda_fixed` keep their lambda = 1 default." Both joint call sites now default those
`NA` entries to `1.0` before calling the solver. The per-column path (`R/fit_baseline.R:857-860`)
already implemented the documented fallback and needed no change.

- File:line: `R/joint_mvn_baseline.R:128-134` (continuous-only joint MVN path).
- File:line: `R/joint_threshold_baseline.R:348-356` (threshold-joint path, only the
  `lambda_family_idx` continuous-family subset is affected).
- Tests: `tests/testthat/test-lambda-dispatch.R`,
  `"[lambda-dispatch] partial lambda_fixed defaults missing columns to 1 (threshold-joint)"` (3
  continuous + binary + categorical fixture, `lambda_fixed` names only `c1`) and
  `"[lambda-dispatch] partial lambda_fixed defaults missing columns to 1 (joint MVN)"` (2
  continuous-only fixture, dispatches to `fit_joint_mvn_baseline()`). Both assert the call no longer
  errors and that the unnamed column(s) match a `lambda_mode = "fixed_1"` fit within `1e-8`.

## 3. Speed: redundant eigendecompositions

**What changed.** Under `lambda_mode = "estimate"` on the joint MVN / threshold-joint path, each
continuous-family column's `R_oo` was factorised THREE times: once in `.mvn_resolve_lambda()`'s
`lambda_block` cache build, once again inside `ml_lambda_for_col()` (a second, separate
`build_pagel_nll_cache()` call for the SAME column), and a third time (dense Cholesky) in
`.mvn_init_per_column()`'s Henderson-centering step via `.mvn_gls_mean_at_lambda()`.

Fix: `build_pagel_nll_cache()` (`R/pagel_lambda.R`) now also returns `$mu_at(lambda)`, a closure
reusing the same eigenbasis to compute the GLS mean at any lambda (the intercept-only branch already
computed this value internally per NLL evaluation; it just wasn't returned). The grid-then-refine
search itself was factored out into `.pagel_lambda_from_cache()` (`R/pagel_lambda.R:285`) so
`.mvn_resolve_lambda()` can call it directly on the SAME cache it already built for the
`lambda_block` search, instead of calling `ml_lambda_for_col()` (which used to rebuild the cache).
`ml_lambda_for_col()` and the new `ml_lambda_and_mu_for_col()` (`R/bm_internal.R:240`) now both
delegate to `.pagel_lambda_from_cache()`, so there is exactly one code path for the search. The
resulting `mu_hat_vec` (one GLS mean per lambda_cols column, `R/joint_mvn_solver.R:410-424`) is
threaded into `.mvn_init_per_column()` (`R/joint_mvn_solver.R:238-273`), which uses it instead of
calling `.mvn_gls_mean_at_lambda()` whenever it is available (finite).

Net effect: one eigendecomposition per continuous-family column instead of three, under the default
`lambda_mode = "estimate"`.

- File:line: `R/pagel_lambda.R:172-268` (`$mu_at` added to `build_pagel_nll_cache()`),
  `R/pagel_lambda.R:285` (`.pagel_lambda_from_cache()`), `R/bm_internal.R:216-256`
  (`ml_lambda_for_col()` / `ml_lambda_and_mu_for_col()`), `R/joint_mvn_solver.R:368-424`
  (`.mvn_resolve_lambda()` reuses `caches` for the per-column estimate), `R/joint_mvn_solver.R:541`,
  `589-593` (`mu_hat_vec` threaded to `.mvn_init_per_column()`), `R/joint_mvn_solver.R:238-273`
  (`.mvn_init_per_column()` consumes `mu_hat_vec`).
- Timing (`fit_baseline(pd, tree, lambda_mode = "estimate")`, K = 3 continuous BM(lambda = 0.3)
  traits, `ape::rcoal(n = 3000)`, 30% missing per column, single run, same tree/data both times):
  - **Before** (speed fix reverted via `git stash` on the 3 touched files, review fixes 1/2/4/5/6
    unaffected): **10.387 s**.
  - **After** (fix applied): **5.993 s**.
  - This is a real but partial speedup (~1.7x), not the full theoretical 3x, because the timing
    also includes the Henderson `S^{-1}` builds (one per distinct rounded lambda value), the Sigma
    initial-value computation, and other per-fit overhead untouched by this change.
- Correctness: `tests/testthat/test-joint-lambda.R`, `"[joint-lambda] fixed_1 is identical within
  1e-12 to the origin/main reference"` still passes (fixed_1 never touches the new code paths, since
  `lambda_vec` is all-1 and the Henderson-centering branch is guarded by `lam_j != 1`).
  `"[joint-lambda] mean model: solver mu matches bm_impute_col at the estimated lambda"` (the
  mutant-catching test from Rose's own review) still passes after the refactor.

## 4. New/rewritten dispatcher-level tests

Two changes to `tests/testthat/test-lambda-default.R`:

- **(a) Added** `"[lambda-default] estimated lambda is applied"`: a lambda_true = 0.3 BM DGP
  (n = 200, 2 continuous traits, ~30% NA), `phylo_signal_gate = FALSE` and `safety_floor = FALSE` (to
  isolate `lambda_mode`'s effect from that separate low-signal safety mechanism -- an
  earlier draft without these flags failed nondeterministically because the phylo-signal gate routed
  the trait to the grand-mean corner regardless of `lambda_mode`, which is a documented, unrelated
  behaviour, not a bug). Asserts `impute(gnn = FALSE)` default ("estimate") and `lambda_mode =
  "fixed_1"` give different `$completed` values, AND that the estimate-mode prediction equals an
  INDEPENDENT `bm_impute_col(y, R, lambda = lambda_hat)` call within `1e-3`.
- **(b) Rewrote** `"[lambda-default] predict rebuild"` (renamed
  `"predict() reproduces the fit-time baseline at estimated lambda"`): now actually calls
  `predict(fit)` on a `fit_pigauto(gnn = FALSE)` fit and compares `pred$imputed_latent` at
  originally-missing cells against the gnn = FALSE blend (`r_cal_bm * baseline_full$mu + r_cal_mean *
  mean_baseline_per_col`, the same formula `test-gnn-off.R` uses) computed from the fit's own stored
  fields. The old version only replayed `lambda_fixed` through `fit_baseline()` against itself and
  never called `predict()`, so it could not catch a mutant that estimates lambda but never applies it
  to a live prediction. The original `lambda_fixed` replay check is retained as a second assertion in
  the same test.

- File:line: `tests/testthat/test-lambda-default.R:74-99` (rewritten "(iii) predict rebuild"),
  `tests/testthat/test-lambda-default.R:155-199` (new "(v) estimated lambda is applied").
- Mutant check (item 8 of Rose's list): confirmed via the ordinal-leak mutant test above (section 1)
  and via Rose's own pre-existing `test-joint-lambda.R` mean-model test, both of which now fail on a
  broken lambda-application path.

## 5. Roxygen corrections

Rewrote the four flagged blocks to match current (S4/S6) dispatch behaviour: ordinal is explicitly
excluded from the continuous-family list wherever it appears; the joint MVN / threshold-joint paths
are documented as USING (not discarding) their own continuous-column output under `"estimate"`; the
covariate path is documented as accepting `"estimate"` / a numeric lambda (only `"cv"`/`"bayes"` fall
back to lambda = 1); and `joint_solver = "inhouse"`'s "byte-identical to prior releases" claim is
scoped to `lambda_mode = "fixed_1"`.

- File:line: `R/fit_baseline.R:27-42` ("Per-type lambda dispatch" Details paragraph).
- File:line: `R/fit_baseline.R:116-135` (`@param joint_solver`, was ~124-128 pre-edit).
- File:line: `R/fit_pigauto.R:226-247` (`@param lambda_mode`, was ~240 pre-edit).
- File:line: `R/impute.R:61-83` (`@param lambda_mode`, was ~75 pre-edit).
- Ran `devtools::document()`: regenerated `man/fit_baseline.Rd`, `man/fit_pigauto.Rd`, `man/impute.Rd`
  (no NAMESPACE change; no new roxygen warnings introduced -- the one pre-existing `[0, 1]` markdown
  link warning at `R/joint_mvn_solver.R`'s `@param lambda` was confirmed present on `HEAD` before any
  of this session's edits via `git stash` + re-run, so it is out of scope here).
- Verified rendering with `tools::Rd2txt()` on both edited sections; no malformed Rd.

## 6. NEWS.md (0.11.0.9000 entry)

Extended the existing "Default flip" section (did not create a new heading) with, and verified each
claim by `grep` before writing it (commands and file:line evidence in parentheses):

- Fixed the existing paragraph's list of discrete types to include ordinal explicitly, and replaced
  the vague "unless the joint delegate inherits the continuous block's `lambda_block`" clause (which
  read as ordinal sometimes NOT being at lambda = 1) with the accurate statement that the opt-in
  `exact` / `joint_refine_iter > 0` paths use `lambda_block` only for cross-trait computation, never
  overriding a discrete/ordinal column's own lambda = 1.
- **Covariate fits now estimate lambda**: `bm_impute_col_with_cov()` accepts `"estimate"`
  (`R/fit_baseline.R:826-838`, and the S4 comment there dates this to feat/joint-lambda-default).
- **`joint_solver = "rphylopars"` runs `model = "lambda"` under `"estimate"`**: verified
  `R/joint_mvn_solver.R:845-846` ("it selects `model = "BM"` vs `model = "lambda"`") and the
  plausibility-guard fallback at `R/joint_mvn_solver.R:853` (per Rose's review citation, unchanged by
  me). Noted the guard also fires under `fixed_1` (`model = "BM"`), matching Rose's finding 1.
- **`cross_validate()` / `compare_methods()` / `simulate_benchmark()` / `multi_impute_trees()` inherit
  the default**: verified by grep --
  - `cross_validate()` (`R/cross_validate.R:100-108`) calls `fit_pigauto(...)` with no `lambda_mode`
    argument.
  - `compare_methods()` (`R/evaluate.R:529` def, `R/evaluate.R:564` call) calls
    `fit_baseline(data, tree, splits = sp)` with no `lambda_mode`.
  - `simulate_benchmark()` (`R/simulate_benchmark.R:106`) calls `fit_baseline(pd, tree, splits = spl)`
    with no `lambda_mode`, then (`R/simulate_benchmark.R:107-109`) passes that already-fitted baseline
    into `fit_pigauto(..., baseline = bl, ...)`; documented the resulting gap (a `lambda_mode` passed
    through `simulate_benchmark(...)`'s `...` reaches `fit_pigauto()`'s recorded
    `model_config$lambda_mode` but never reaches the already-fitted `bl`, so the two can disagree).
  - `multi_impute_trees()` was already covered by the pre-existing top paragraph
    (`lambda_mode = baseline_arg("lambda_mode", "estimate")`, `R/multi_impute_trees.R:562`).
- **Known limitation**: `suggest_next_observation()` (`R/active_impute.R:473`, read-only per the
  "do not touch `R/active_impute.R`" instruction) reads `fit$graph$R_phy` directly (lambda = 1) and
  never `model_config$lambda_per_trait`.
- No em dashes used (checked with `grep -n "—"`, no hits in the touched NEWS section).

## Verification (verbatim)

```
SUITE FAIL 0 | PASS 2594
CHECK errors 0 warnings 0 notes 1
```

`docs/dev-log/lambda-default/S6-suite-final.log` and `docs/dev-log/lambda-default/S6-check-final.log`
hold these two lines. The suite command was run with `NOT_CRAN=true` via `devtools::load_all()` +
`testthat::test_dir()`; the check command was `rcmdcheck::rcmdcheck(args = c("--as-cran",
"--no-manual"), error_on = "never")`, which internally re-runs the full test suite a second time
(without `NOT_CRAN=true`, so slow/integration tests are skipped there instead: that inner run showed
`FAIL 0 | WARN 308 | SKIP 70 | PASS 2295`, consistent with the outer summary). The one CHECK note is
the standard dev-version note (`Version contains large components (0.11.0.9000)`), not new evidence of
a problem.

## Deviations from the brief

- Item 3's timing measurement used `ape::rcoal(n = 3000)` with a single seed/run each for before/after
  (not multiple reps), matching the brief's ask for "the timing numbers" rather than a full benchmark;
  wall-clock on a shared, busy machine, so treat the ~1.7x figure as directional, not a precise ratio.
- Item 4(a)'s test needed `phylo_signal_gate = FALSE, safety_floor = FALSE` beyond what the brief's
  one-line description specified, to avoid a real (unrelated, pre-existing) confound: pigauto's
  phylo-signal gate can route a trait to the grand-mean corner independent of `lambda_mode` when
  estimated signal is weak, which would make `completed` values coincidentally match between
  `"estimate"` and `"fixed_1"` for reasons having nothing to do with the bug this test targets.
- Item 6 fixed one additional pre-existing inaccuracy in the SAME paragraph the brief pointed at (the
  "unless the joint delegate inherits `lambda_block`" clause) while adding the requested content,
  since leaving it would have contradicted the new text directly above it.
- Did not touch the "profile REML" vs "per-column ML estimator" naming inconsistency Rose flagged
  (her required change 4, naming) or the internal "NULL until S4 adds these" comments (her finding 7)
  -- neither was in the four specific line ranges or four bullet points assigned to me for items 5/6.
- Did not modify `R/active_impute.R`, `docs/`, `BACE/`, or `script/campaign_*`, per the task's explicit
  fence. `docs/dev-log/lambda-default/benchmark.md`, `prerun.md`, and `real-data.md` show diffs in
  `git status` from a concurrent lane in this shared worktree (correcting the claims Rose's review
  item 4 flagged); I did not make or review those edits.
- Did not commit or push; per instructions, that is the orchestrator's job.
