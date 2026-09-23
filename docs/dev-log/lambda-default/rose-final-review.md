# Rose final review: PR #187 (feat/joint-lambda-default)

Reviewer: Rose (fresh adversarial pass, did not build this). Date: 2026-09-23.
Reviewed commit: HEAD 6bd25b6 (merge-base ab02e31; origin/main 2470ac2 differs only in docs).
All probes ran on clean `git archive` copies of ab02e31 and HEAD in a scratch directory; the worktree was
not modified except for this file. Scratch scripts: `fixed1_cases.R`, `est_cases.R`, `ord.R`, `mo.R`,
`se_check.R`, `attrib.R`, `timing.R` (session scratchpad, not committed).

## Overall verdict: NEEDS-CHANGES

The core is sound. The fixed_1 identity holds, and I checked it more widely than the PR's fixture does.
The estimator and the predictor share one mean model, and the solver-level lambda path is correct. The
changes needed are about truthfulness and edges: ordinal traits leak through the new default while the
docs and `lambda_per_trait` say they do not, `lambda_fixed` breaks on a partial vector, several roxygen
paragraphs now describe the old behaviour, and the suite and check evidence predate the last R commit.

## 1. fixed_1 identity: CONFIRMED

- The fixture is a real oracle, not a self-consistency check. I re-ran `fit_baseline(pd, ref$tree,
  splits = ref$splits)` and `fit_mvn_bm_inhouse(ref$L, ref$tree)` on a `git archive ab02e31` copy and got
  max |diff| = 0 against `tests/testthat/fixtures/lambda_fixed1_reference_ab02e31.rds` for mu, se,
  anc_recon, anc_var and phylocov. It is used in `tests/testthat/test-joint-lambda.R:38-66` (solver) and
  `tests/testthat/test-lambda-dispatch.R:154-164` (dispatcher). Both pass at HEAD.
- Fixture coverage is narrow. Its `df` routes to `threshold_joint` + `ovr_categorical` only. So I compared
  ab02e31 against HEAD at `lambda_mode = "fixed_1"` on seven more configurations: continuous-only
  `joint_mvn`, single-trait `per_column_bm`, covariates, zi_count, multi_proportion, multi-obs
  (`species_col`), and end-to-end `impute(gnn = FALSE)` (completed, se, conformal bounds). Max |diff| = 0
  in every one.
- New code that runs under fixed_1: the `.mvn_resolve_lambda()` shortcut (R/joint_mvn_solver.R:316-322),
  the Henderson cache keyed "1.000" that reuses the untransformed tree (R/joint_mvn_solver.R:538), and the
  extra return fields. None of it changes numbers. One real behaviour change under
  fixed_1: the Rphylopars plausibility guard (R/joint_mvn_solver.R:853, `fit_joint_solver()`) also fires for
  `model = "BM"`. So `joint_solver = "rphylopars", lambda_mode = "fixed_1"` now falls back to in-house
  when phylopars returns tip predictions more than 10x the observed range. That is defensible, but NEWS
  does not mention it.

## 2. Correctness of the lambda path: WEAK (one WRONG sub-finding)

Confirmed:
- Mean model. At the estimated lambda, sparse (Henderson on `transform_tree_pagel`, GLS-centred) and
  dense (`bm_impute_col`) mu agree to 1e-7 to 1e-9 for lambda-hat in {0.19, 0.35, 0.58, 0.68, 0.995}, on
  data shifted by +1.5 (`se_check.R`).
- `lambda_cols` mapping (R/joint_threshold_baseline.R:346-363). Only
  continuous/count/proportion columns are eligible. On a mixed fixture, b1, o1 and k1 report 1, and
  binary/categorical mu and se are identical (max |diff| = 0) between "estimate" and "fixed_1".
- `lambda_fixed` full-vector rebuild reproduces mu/se exactly (max |diff| = 0) on the mixed
  threshold-joint fixture and on the covariate per-column path.
- Multi-obs `impute(species_col =, gnn = FALSE)` runs under the default: no NA in completed, lambda 0.46
  on the signal trait. zi_count runs, with the magnitude column estimated via the per-column path and the
  gate at 1. multi_proportion runs, rows still sum to 1, and per-component lambda is estimated. The
  non-finite conformal bounds on the binary column are by design and identical under fixed_1.

WRONG: ordinal traits are not held at lambda = 1 under the new default.
- The per-trait ordinal path selection calls `bm_impute_col(X[, col], R_phy_local, lambda = bm_lambda)`
  at R/fit_baseline.R:606, and `bm_lambda` is now "estimate" by default. On 8 seeds (n = 100, one
  continuous trait plus one K = 3 ordinal, `ord.R`), the ordinal baseline mu changed in 5 of 8 (max
  |delta mu| 1.6 to 2.05 z-units). In 3 of those 5 the chosen path flipped (bm_mvn against lp).
  `lambda_per_trait["o1"]` still reports 1 in every case.
- This contradicts the `fit_baseline()` roxygen for `lambda_mode` ("NOT ordinal") and the
  `fit_joint_threshold_baseline()` roxygen ("Binary and ordinal liability columns always stay at lambda
  = 1"). It also contradicts the "Per-type lambda dispatch" Details. The block only runs when `splits` is
  non-NULL and data are single-obs, which is the held-out baseline that GNN-on fits train and calibrate
  on. `lambda_fixed` does not reach this call either.
- A second route: when multi_proportion is present, ordinal falls to the per-column path and gets an
  estimated lambda (0.995 in `est_cases.R`). That is pre-existing, but it is now the default.

WRONG (small): a partial `lambda_fixed` errors on every joint path.
- The roxygen says "Columns not present in `lambda_fixed` keep their lambda = 1 default". But
  `fit_baseline(pd, tree, lambda_fixed = c(c1 = 0.5))` on a two-continuous-trait fixture (and on one with
  an added binary) stops with "fit_mvn_bm_inhouse: numeric 'lambda' vector must have all entries in [0,
  1]". The cause is `unname(lambda_fixed[colnames(L_in)])` giving NA at R/joint_mvn_baseline.R:129
  and at R/joint_threshold_baseline.R:352. Only the per-column branch (R/fit_baseline.R:849-850)
  implements the documented fallback.
- A fit made with `joint_solver = "rphylopars"` stores `lambda_per_trait = NA`, so it cannot be replayed
  through `lambda_fixed` either, because `fit_baseline()` rejects non-finite entries.

## 3. Uncertainty: CONFIRMED, with one gap

- Per-column se at lambda < 1 comes from the Henderson precision on the Pagel-transformed tree
  (`get_henderson_at(lambda_vec[j])`), and the dense fallback transforms R. Both are on R(lambda).
  Conformal scores come from validation residuals and do not depend on lambda.
- Sparse and dense se differ by a constant factor per column: 0.86 to 1.52 at lambda-hat of about 1,
  0.97 to 1.28 at 0.6 to 0.7, and 1.00 to 1.10 at 0.2 to 0.35. The cause is `henderson_bm_predict()`
  using `var(y_obs)` for sigma^2 (R/henderson_s_inv.R) while `bm_impute_col()` uses the GLS REML sigma^2.
  This predates the PR (it is present at lambda = 1) and the gap shrinks as lambda falls. It is not a new
  bug, but it should be recorded as a known limitation.
- Opt-in paths (`predict_method = "exact"`, `joint_refine_iter > 0`) use `lambda_block` for every
  column, including binary columns. This is documented.
- Gap: `suggest_next_observation()` (R/active_impute.R:473-515) computes variance reduction from
  `fit$graph$R_phy` at lambda = 1 and never reads `model_config$lambda_per_trait`. Under the new default,
  the imputation baseline uses R(lambda-hat) but the sampling-design helper (selling point (e)) scores
  candidates under R(1). Before this PR both used lambda = 1 by default. The mismatch is now the default.

## 4. Tests that can fail: WEAK

I built a mutant that estimates and reports lambda but predicts at lambda = 1 (Henderson at 1, no
centring, in R/joint_mvn_solver.R) and a covariate mutant that computes lambda-hat but never applies the
transform (R/bm_internal.R:332). Then I ran the four new test files.
- Caught: `[joint-lambda] mean model: solver mu matches bm_impute_col at the estimated lambda`
  (test-joint-lambda.R) and `[lambda-cov] lambda = 0 matches OLS-with-covariates` (test-lambda-covariates.R).
- Not caught: all 6 tests in test-lambda-dispatch.R and all 5 in test-lambda-default.R pass on the mutant.
  In particular `[lambda-default] predict rebuild` (test-lambda-default.R) only checks that
  `lambda_fixed` replays itself, which is trivially true if lambda is silently 1. Despite its name, it
  never calls `predict()`. `[lambda-dispatch] lambda_per_trait populated` checks the reported value, not
  that it reaches mu. `[lambda-dispatch] covariates keep lambda` would also pass if the column were
  absent (`isTRUE(NA == 1)` is FALSE).
- So one solver-level test guards "lambda reaches the prediction". No dispatcher-level test checks that
  `fit_baseline(lambda_mode = "estimate")` mu differs from fixed_1 on lambda < 1 data by more than a
  bound. No test covers the ordinal leak in finding 2.

## 5. Scope fence: CONFIRMED

`git diff --stat origin/main...HEAD` touches nothing under `script/campaign_*`, `BACE/` or
`docs/dev-log/arc/2026-09-2*`. `script/bench_lambda_datasets.R` is new and lives in the build-ignored
`script/`.

## 6. Claims against evidence: WEAK (one WRONG)

Traced and correct:
- The real-data z-RMSE table matches `lambda_datasets_summary.csv`. "Never worse by more than 1%" holds
  (max +0.9%).
- The simulation no-GNN table re-derives exactly from `core_lambda_fast_agg_summary.csv` via
  `compare_benchmark.R`, using the committed yardstick in the imputation-sim worktree. That yardstick is
  not in this PR, so the "committed" column cannot be reproduced from the PR alone.
- The NEWS figure "~0.38" traces to S5-default-report.md.

Problems:
- WRONG: "Discrete accuracy is identical under every lambda setting" (PR body and real-data.md). The
  report's own table (`lambda_datasets_report.txt`, discrete accuracy block) shows AVONET 2000 migration
  at 0.820 under fixed_1 and 0.823 under the other three modes. prerun.md also says "Ordinal accuracy
  moves by at most 0.004". This agrees with finding 2. Say that binary and categorical are identical and
  that ordinal can move.
- The suite (FAIL 0 / PASS 2565) and R CMD check (0/0/1) logs were committed with b8eb384. Commit 379a95d
  then changed `R/joint_mvn_solver.R` (the Rphylopars guard) and added a test, so neither log is evidence
  for HEAD. My re-run at HEAD: the eight lambda and joint test files pass (FAIL 0, 250 expectations,
  about 27 s). Full-suite re-run: see the addendum at the end of this file. R CMD check was not re-run.
- "12 of 13 helps or ties": 11 improve, 1 is +0.1% (called a tie), 1 is +0.9%. Acceptable, but "11
  improve, 1 tie, 1 loss" is plainer.
- "Coverage stays at 0.957 to 0.960 for n >= 300": for gnn_off at n >= 300 the range is 0.957 to 0.969
  (the lambda = 1 cells are 0.968 and 0.969). The n = 100 cells are 0.86 to 0.90, unchanged from the
  committed run, but worth one clause.
- The pre-registered gate failed on two counts, not one. The PR names the lambda 0.3 half-gap failure,
  but the lambda = 1 |delta| < 0.01 guard also failed (in the improving direction), and
  `compare_benchmark.R` prints BENCH_FAIL. benchmark.md says so; the PR body should too.
- Naming of the estimator. The PR says "profile REML" and NEWS says "per-column ML estimator". The code
  (R/pagel_lambda.R, `build_pagel_nll_cache()`) is a profile likelihood with an (n - p) variance
  denominator that drops the log|X' R^-1 X| term. Its own comment says "not full REML". Use one accurate
  name everywhere.
- The flip changes `impute()`'s default (`gnn = TRUE`), but the GNN-on evidence is one 20-seed pre-run
  cell. The PR says the GNN wave is still running, so this is not hidden. Merge should wait for that wave.
- An attribution check I ran (`attrib.R`, pure continuous BM, n = 300, 40 seeds) found no mean-model
  confound at lambda = 1. There, "estimate" is about 2% worse than fixed_1 (0.170 against 0.166), so the
  campaign's lambda = 1 gains come from elsewhere in the types_mixed DGP. That is not a defect, but "the
  lambda = 1 cells also improve" should not be read as "estimate is free at lambda = 1".

## 7. What a CRAN or package reviewer would reject: WEAK

- Stale documentation that now states the opposite of the code:
  - R/fit_pigauto.R:239-242 and R/impute.R:74-77 (and man/*.Rd) say the covariate path "has no lambda
    argument and always fits at lambda = 1". It now estimates lambda.
  - R/fit_baseline.R:27-42 (Per-type lambda dispatch) says the threshold-joint continuous output "is
    discarded in favour of the lambda-aware per-column BM fit". That is no longer true for "estimate". The
    same paragraph lists ordinal both as lambda-governed and as held at 1.
  - R/fit_baseline.R:124-128 (`joint_solver`) says `lambda_mode != "fixed_1"` disables the joint MVN path
    "(it has no lambda argument)". It also says the in-house solver "is byte-identical to prior releases",
    which is not true under the new default.
  - Internal comments still say "NULL until S4 adds these" (R/fit_pigauto.R:697 and the matching
    comment in the GNN-on model_config block; R/multi_impute_trees.R:527, 571, 676).
- NEWS omits user-visible changes. The covariate baseline now estimates lambda. `joint_solver =
  "rphylopars"` now runs `model = "lambda"` by default, which benchmark.md measures as 26 times slower
  with occasional explosive seeds. `cross_validate()`, `compare_methods()` and `simulate_benchmark()` call
  `fit_baseline()` without `lambda_mode` (R/evaluate.R:564, R/simulate_benchmark.R:106), so their
  baselines change silently. In `simulate_benchmark()`, passing `lambda_mode = "fixed_1"` through `...`
  now gives an "estimate" baseline paired with a fit labelled fixed_1. NEWS also has no known-limitations
  line (downward lambda bias, the weak-signal gate shortfall, LepTraits flight duration).
- Scaling and usability (D-139). "estimate" replaces the O(n) sparse default with dense O(n^3)
  eigendecompositions of R_oo, done twice per continuous column: once for the block caches in
  `.mvn_resolve_lambda()` (R/joint_mvn_solver.R:356) and again inside `ml_lambda_for_col()`. Measured for `fit_baseline()` with K =
  3 and 30% missing: n = 1000 went from 0.1 s to 1.1 s, and n = 3000 from 0.3 s to 8.3 s. Cubic growth
  puts n = 10,000 at minutes per baseline, with multi-GB memory, repeated per tree in
  `multi_impute_trees()`. The timing evidence stops at 2,000 species. Under the default `per_column`
  predict route, the block lambda is used only for Sigma, which that route does not consume, so the first
  pass is mostly wasted.
- Minor: `bm_impute_col_with_cov(lambda = "estimate")` uses a bare `optimize()` on [0.01, 0.99]
  (R/bm_internal.R:323) and skips the grid-then-refine guard that `ml_lambda_for_col()`
  added against boundary plateaus. The covariate and non-covariate estimators therefore differ slightly.
- Not a problem: the new tests run in about 27 s, the heavy recovery tests are `skip_on_cran()`, the
  fixture is 7 KB, and `docs/`, `script/` and `dev/` are build-ignored.

## Required changes

1. Stop the ordinal leak or document it. Either pass `lambda = 1` (or the column's `lambda_fixed` value)
   to the ordinal bm_mvn alternative at R/fit_baseline.R:606 so ordinal really stays at lambda = 1, or
   change the roxygen and report the lambda actually used in `lambda_per_trait`. Add a test that ordinal
   mu is identical under "estimate" and "fixed_1" (or not, matching the chosen contract) on a fixture
   where bm_mvn wins.
2. Make a partial `lambda_fixed` default missing continuous-family columns to 1 on the joint paths
   (R/joint_mvn_baseline.R, R/joint_threshold_baseline.R), as documented. Add a test.
3. Correct the stale roxygen and Rd listed in finding 7, and the internal "NULL until S4" comments.
4. Correct the claims. Say "binary and categorical accuracy identical; ordinal can move (AVONET 2000
   migration 0.820 to 0.823)". Give the coverage range as 0.957 to 0.969. State both gate failures. Use
   one accurate name for the estimator in the PR, NEWS and roxygen.
5. Re-run the full suite and `R CMD check --as-cran` at the final HEAD and commit those logs. The
   committed logs predate 379a95d.
6. Extend NEWS: the covariate path now estimates lambda, `rphylopars` now uses `model = "lambda"` by
   default (cost and the fallback guard, which also applies under fixed_1),
   `cross_validate()`/`compare_methods()`/`simulate_benchmark()` inherit the new default, and a
   known-limitations line.
7. Either make `suggest_next_observation()` use the fit's per-trait lambda, or document that it scores
   candidates under lambda = 1 even when the baseline estimated lambda.
8. Add one dispatcher-level test that fails if lambda is estimated but not applied: `fit_baseline(...,
   lambda_mode = "estimate")` mu on a lambda = 0.3 fixture must differ from fixed_1 by more than a stated
   bound, or match `bm_impute_col(lambda = lambda-hat)`.
9. Before merge (not before review): wait for the GNN-on wave to land in benchmark.md. Record the
   large-n cost (my timing above, or a measured n = 5,000 point) and consider reusing the per-column
   caches so each column's eigendecomposition happens once.

## Addendum: full-suite re-run at HEAD

Full `testthat::test_dir()` on a `git archive HEAD` copy (R, tests, DESCRIPTION, NAMESPACE, inst, data,
NEWS.md): FAIL 4, PASS 2809, 66 files. All 4 failures are in `test-community-surface.R`, which reads
README.md, `_pkgdown.yml`, the vignettes and man/*.Rd. Those files were not in my archive copy. Re-run
read-only in the full worktree at HEAD, that file passes 8 of 8. So the suite is green at HEAD, with
that caveat. R CMD check was not re-run at HEAD (required change 5 stands). This review file sits under
a gitignored path (`git status --ignored` shows it as `!!`), so it needs `git add -f` if it is to be
committed.
