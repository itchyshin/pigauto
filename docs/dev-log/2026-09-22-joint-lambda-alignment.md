# Pagel's lambda in the joint baseline: symbolic alignment and the shared-vs-per-trait decision

Lane: feat/joint-lambda-default (brain D-278). Spec: `specs/2026-05-18-pagel-lambda-baseline-design.md`,
decisions 1 and 3. Written before any solver code, as the lane GOAL requires.

## 1. What exists and what is missing

| piece | file | lambda today |
|---|---|---|
| per-column BM kernel | `R/bm_internal.R::bm_impute_col(y, R, nugget, lambda)` | numeric, `"estimate"`, `"cv"`, `"bayes"` |
| profile REML NLL cache | `R/pagel_lambda.R::build_pagel_nll_cache(y, R)` | one eigendecomposition of `R_oo`, then O(n_o) per lambda |
| tree transform | `R/pagel_lambda.R::transform_tree_pagel(tree, lambda)` | internal edges x lambda, terminal edges + (1 - lambda) h_parent |
| sparse precision | `R/henderson_s_inv.R::build_henderson_S_inv(tree)` | Hadfield and Nakagawa (2010) eq. 29, O(n) |
| per-column sparse predictor | `R/henderson_s_inv.R::henderson_bm_predict(y, henderson, cor_scale = TRUE)` | none (lambda = 1) |
| joint solver | `R/joint_mvn_solver.R::fit_mvn_bm_inhouse(L, tree, ...)` | none, and no likelihood to profile one |
| joint dispatcher | `R/fit_baseline.R:357 force_per_column <- lambda_mode %in% c("estimate", "cv", "bayes")` | any estimated lambda drops continuous columns off the joint path |
| covariate BM | `R/bm_internal.R::bm_impute_col_with_cov(y, X, R, nugget, ridge, lrt_threshold)` | none; warning at `fit_baseline.R:731` says so |

The gap is therefore not the lambda estimator (it exists, per column) but the joint solver and the two
paths (joint, covariate) that never receive one.

## 2. Symbolic model and its code counterpart, term by term

One BM-eligible latent column k, observed at tips o, missing at tips m, on the preprocessed (z-scored) scale.

| symbol | meaning | code |
|---|---|---|
| A = vcv(tree) | shared root-to-MRCA path lengths | `ape::vcv(tree)` |
| R = cov2cor(A) | phylogenetic correlation | `phylo_cor_matrix(tree)` |
| R(lambda) = lambda R + (1 - lambda) I | Pagel's transform on the correlation scale | `lambda * R; diag <- lambda * diag + (1 - lambda)` (`bm_impute_col` lines 107-111) |
| Var(y_k) = sigma_k^2 [lambda_k R + (1 - lambda_k) I] | phylogenetic variance sigma_k^2 lambda_k beside residual variance sigma_k^2 (1 - lambda_k) | new `lambda` argument of `fit_mvn_bm_inhouse` |
| mu_k(lambda) = (1' R(lambda)^-1 1)^-1 1' R(lambda)^-1 y_o | GLS phylogenetic mean | `build_pagel_nll_cache`: `mu_hat <- sum_b / sum_a` |
| sigma_k^2(lambda) = e' R_oo(lambda)^-1 e / (n_o - 1) | REML variance at fixed lambda | same closure, `sigma2` |
| NLL_k(lambda) = 0.5 [(n_o - 1) log sigma_k^2(lambda) + log det R_oo(lambda)] | profile REML | `cache$nll(lambda)` |
| lambda_k = argmin NLL_k over [0.01, 0.99] | per-trait estimate | `ml_lambda_for_col` (exists); reused |
| lambda_bar = argmin sum_k NLL_k(lambda) | block estimate over continuous-family columns | new `.block_lambda_ml()` summing the K caches; one `optimize()` |
| E[y_m given y_o] = mu + R_mo(lambda) R_oo(lambda)^-1 (y_o - mu) | conditional mean | `henderson_bm_predict` on the transformed tree, or `bm_impute_col(lambda=)` dense |
| Var[y_m given y_o] = sigma^2 (1 - h_i), h_i = diag(R_mo R_oo^-1 R_om) at lambda | conditional variance | same functions, `se` |
| vec(L) ~ MVN(0, Sigma (x) R(lambda_bar)) | joint block with one shared lambda | `build_henderson_S_inv(transform_tree_pagel(tree, lambda_bar))` feeds `.mvn_sigma_kron_M`, `exact_conditional_mvn`, `.mvn_estep_refine` |
| with covariates: y = X beta + u, beta(lambda) = (X' R(lambda)^-1 X)^-1 X' R(lambda)^-1 y | GLS in the profile | extend the cache with `c_X = U' X`; `bm_impute_col_with_cov(lambda=)` |

## 3. The tree transform is exact on the correlation scale, for any tree

`transform_tree_pagel` scales every internal edge by lambda and lengthens each terminal edge by
(1 - lambda) times the root-to-parent depth. For tips i and j with most recent common ancestor at depth
d_ij, every edge on the root-to-MRCA path is internal, so `vcv(T(lambda))[i, j] = lambda d_ij = lambda A[i, j]`.
For the diagonal, root-to-tip depth is `lambda h_parent + t + (1 - lambda) h_parent = h_parent + t = A[i, i]`.
Hence

    cov2cor(vcv(T(lambda)))[i, j] = lambda A[i, j] / sqrt(A[i, i] A[j, j]) = lambda R[i, j]   (i != j),

and 1 on the diagonal. So `cov2cor(vcv(T(lambda))) = lambda R + (1 - lambda) I` exactly, ultrametric or
not. The comment in `R/pagel_lambda.R` calling this "a close approximation for non-ultrametric trees" is
right on the covariance scale (`vcv(T(lambda)) = lambda A + (1 - lambda) diag(A)`, not `lambda A + (1 - lambda) I`)
but over-cautious on the correlation scale the solver actually uses; this lane corrects the comment.

Consequence: Hadfield and Nakagawa's sparse Q built on T(lambda) gives R(lambda)^-1 b at tips in O(n) via
`henderson_R_inv_apply(..., cor_scale = TRUE)`, with the stored `tip_sqrt_d` unchanged because tip depths are
preserved. No dense n x n matrix is needed at any lambda.

## 4. Shared lambda per block versus per-trait lambda: the decision

Three candidates were weighed.

1. **Shared lambda_bar for the whole block.** Keeps `Sigma (x) R(lambda_bar)` exactly; one Henderson build.
   Cost: a weak-signal trait shrinks a strong-signal trait toward the grand mean and vice versa. The spec's
   motivating case (BIEN: sla weak, height strong) is exactly where this hurts, and spec 4.5 asks for
   `lambda_per_trait`.
2. **Per-trait lambda_k everywhere.** The joint covariance becomes `Sigma_P (x) R + Sigma_E (x) I`, which
   has no Kronecker form; estimating Sigma_P and Sigma_E needs a two-matrix REML over an nK system. A spec
   non-goal ("joint estimation across traits ... separate lambda per off-diagonal").
3. **Hybrid (chosen).** Per-trait lambda_k on continuous-family columns (continuous, count, ordinal,
   proportion, zi magnitude), used by the default prediction path, which is already per column
   (`henderson_bm_predict` per column; `max_iter = 0`, `predict_method = "per_column"`). One block
   lambda_bar wherever a single common R is required: the Sigma M-step, the opt-in exact conditional and
   EM refinement, and the discrete liability columns (binary, zi gate, ordinal, OVR synthetic columns).

Why liability columns take lambda_bar rather than their own estimate: their entries are plug-in posterior
means `E[L | y]` under an N(0, 1) prior, so a lambda profiled on them measures the shrinkage of the E-step,
not the phylogenetic signal of the liability. Inheriting the continuous block's lambda_bar is the spec's
"free upgrade" (decision 4) stated honestly. When a block has no continuous-family column, lambda_bar = 1
and the discrete path is bit-identical to today.

Costs of the hybrid: at most K Henderson builds (cached by unique lambda rounded to 1e-3), each O(n), plus
one 1-D optimisation per column and one for the block, each O(n_o) per evaluation after the existing
one-off eigendecomposition. The `Sigma` estimate uses lambda_bar and is therefore a shared-lambda Kronecker
MLE even when the columns' own lambda_k differ; this is recorded on the fit as `lambda_block` next to
`lambda_per_trait` so a reader can see both.

## 5. Guards carried over from the spec (section 6.1)

Optimiser clipped to [0.01, 0.99]. Columns with fewer than 10 observed cells take lambda_bar (or 1 when
lambda_bar is undefined). `lambda_mode = "fixed_1"` bypasses every new branch and must remain bit-identical
to the current output; a test asserts this at 1e-8. `"cv"` and `"bayes"` stay per-column modes: they have
no joint analogue and keep forcing the per-column path.

## 6. What changes where (the build slices)

- `R/joint_mvn_solver.R`: `fit_mvn_bm_inhouse(..., lambda = "fixed_1")` accepting `"fixed_1"`, `"estimate"`
  or a numeric scalar; returns `lambda_per_trait`, `lambda_block`; `fit_joint_solver(lambda=)`;
  `.fit_mvn_bm_rphylopars` passes `model = "lambda"` when lambda is estimated.
- `R/bm_internal.R`, `R/pagel_lambda.R`: `bm_impute_col_with_cov(lambda=)`, cache with a design matrix.
- `R/fit_baseline.R`: `force_per_column` only for `"cv"` / `"bayes"`; lambda threaded to the joint,
  threshold-joint and OVR delegates and to the covariate branch; the "runs at lambda = 1" warning goes.
- `R/fit_pigauto.R`, `R/impute.R`, `R/multi_impute.R`, `R/multi_impute_trees.R`, `R/predict_pigauto.R`:
  default `"estimate"`; `model_config$lambda_mode`, `$lambda_per_trait`, `$lambda_block`.

## Review verdict (Rose, 2026-09-22)

Read: this note; `specs/2026-05-18-pagel-lambda-baseline-design.md`; `R/pagel_lambda.R`;
`R/henderson_s_inv.R`; `R/joint_mvn_solver.R`; `R/fit_baseline.R`;
`R/joint_threshold_baseline.R:290-440`; `tests/testthat/test-lambda-per-type.R`; the lane plan
`~/.claude/plans/eager-weaving-pike.md`. Items marked MEASURED were run this session against the
worktree sources (R 4.6.0, ape, Matrix) in throwaway scratch scripts; nothing in the repo was edited
except this section.

**Overall: PROCEED-WITH-CHANGES.** The algebra is right, the per-trait lambda on continuous columns
is buildable, and the default flip is defensible. Three things must change before S2 is dispatched:
the discrete-liability inherit (B iii) must be cut or re-scoped, the fixed_1 oracle (C) must be
anchored to pre-change output, and three gates (E) must be changed from "different" to "not worse".

### A. `cov2cor(vcv(T(lambda))) = lambda R + (1 - lambda) I` exactly, for any tree. SOUND.

The derivation in section 3 is correct and the code supports it. Every edge on a root-to-MRCA path
has an internal node as its child, so all of them are scaled by lambda at `pagel_lambda.R:252-253`,
giving `vcv(T)[i,j] = lambda A[i,j]`; the terminal compensation at `pagel_lambda.R:256-258` adds
exactly `(1 - lambda) * node_depths[parent]`, so root-to-tip depth returns to `lambda h_p + t +
(1 - lambda) h_p = A[i,i]`. MEASURED on four tree classes at lambda in {0.9, 0.5, 0.3, 0.05, 0.01}:
`max|cov2cor(vcv(T(lambda))) - (lambda R + (1 - lambda) I)| <= 3.3e-16` for an ultrametric
`rcoal(60)`, a non-ultrametric `rtree(60)`, a non-ultrametric tree with exponential edge lengths, and
a tree carrying a zero-length internal edge. The comment at `pagel_lambda.R:224-227` is over-cautious
and this note is right to correct it.

The Henderson `cor_scale` path stays valid. MEASURED at n = 300: `tip_sqrt_d` rebuilt from
`T(lambda)` differs from the original by <= 1.11e-16 at every lambda, and
`henderson_R_inv_apply(b, build_henderson_S_inv(T(lambda)), cor_scale = TRUE)` matches
`solve(lambda R + (1 - lambda) I, b)` to a relative error of 2e-14 to 5e-12, on both ultrametric and
non-ultrametric trees, with **no** conditioning degradation as lambda falls to 0.01 (I expected the
1/lambda blow-up in Q and it does not materialise at these scales). That was the main numerical risk
in section 3 and it is clear.

One precision nit, not a verdict driver: section 4's "at most K builds, each O(n)" is wrong.
`build_henderson_S_inv` opens with `tip_depths <- diag(ape::vcv(tree))` (`henderson_s_inv.R:66`),
which forms the dense n x n covariance. The build is O(n^2) in time and memory, not O(n). It is
still cheap at the scales in play (MEASURED ~8 ms at n = 300), but the note should say O(n^2) so a
later reader sizing this to n = 20000 is not misled.

### B. The hybrid decision.

**(i) "the default prediction path is per-column, so per-trait lambda_k breaks nothing that runs by
default". WEAK.** True as plumbing, not as statistics. The default is indeed per-column
(`joint_mvn_solver.R:417` returns at `max_iter <= 0L`, `predict_method = "per_column"` at
`:294`/:566-573`), so no Kronecker object has to absorb K different R matrices. But for K >= 2 that
per-column path is `henderson_bm_predict` (`joint_mvn_solver.R:202-206`), which assumes a **zero root
state** by construction, as the comment at `joint_mvn_solver.R:315-320` states outright. The per-trait
lambda_k this note proposes comes from `build_pagel_nll_cache`, whose profile **estimates a free GLS
mean** (`mu_hat <- sum_b / sum_a`, `pagel_lambda.R:198-200`). At lambda = 1 the two models nearly
coincide, which is why the existing "identical to ~1e-3" claim has held. They come apart as lambda
falls, which is exactly the regime this lane exists to serve. MEASURED, n = 200, 30% missing, BM DGP
at the stated lambda, comparing `henderson_bm_predict` on `T(lambda)` against
`bm_impute_col(y, R, lambda = lambda)`:

| lambda | column mean | max abs difference in mu | RMSE henderson | RMSE dense |
|---|---|---|---|---|
| 1.0 | 1 | 0.001 | 0.070 | 0.070 |
| 0.7 | 1 | 0.071 | 0.536 | 0.533 |
| 0.3 | 1 | 0.246 | 0.788 | 0.777 |
| 0.1 | 1 | 0.439 | 1.012 | 0.984 |

The zero-root predictor is measurably worse, and the gap grows monotonically as lambda falls. "The
column is z-scored so its mean is zero" does not dispose of this: split masking leaves the *observed
subset* off-centre, and liability columns are not centred at zero after the E-step either. Fix
before S2: either estimate lambda_k with the mean pinned at zero (drop the `mu_hat` profile in the
block cache, so the estimator matches the predictor), or centre each column by its lambda-dependent
GLS mean before the Henderson call and add it back after. Whichever is chosen, add a gate: on a
lambda = 0.3 DGP with a non-zero observed-subset mean, the Henderson path must match
`bm_impute_col(y, R, lambda = lambda_hat)` within a stated tolerance.

**(ii) lambda_bar for the Sigma M-step while columns use lambda_k. SOUND, with one disclosure owed.**
Nothing a default user sees depends on Sigma. `em_iterations` defaults to `0L`
(`fit_baseline.R:172`), `joint_refine_iter` to `0L`, `predict_method` to `"per_column"`, so
`.mvn_sigma_kron_M`'s output is computed and then discarded: it reaches `.mvn_estep_refine` only when
`max_iter > 0` (`joint_mvn_solver.R:452`), `exact_conditional_mvn` only under
`predict_method = "exact"` (`:398`), and `extract_liability_variances` only through the EM wrappers
(`joint_threshold_baseline.R:585`, `ovr_categorical.R:263`). So the inconsistency is real but inert
by default. It is not invisible, though: `pars$phylocov` is returned on the fit and a user can print
it, and under any of the three opt-ins the mismatch becomes load-bearing. Section 4 already records
`lambda_block` next to `lambda_per_trait`, which is the right disclosure; add one sentence to the
roxygen for `predict_method = "exact"` and `joint_refine_iter` saying that those paths use
lambda_block, not the per-trait values.

**(iii) Liability columns inherit lambda_bar. WRONG.** Three separate problems.

First, the identification with the spec is false. Section 4 calls this "the spec's 'free upgrade'
(decision 4) made explicit". Decision 4 (spec line 205) describes discrete traits inheriting the
lambda that the **joint fit over the whole liability matrix** estimates, because those paths are
"deterministic delegates of the joint MVN"; under the spec's architecture that is
`phylopars(model = "lambda")` on all columns at once, liabilities included. This note's lambda_bar is
`argmin` of the summed NLLs over **continuous-family columns only** (note line 36, plan line 99), then
imposed on the liabilities. Those are different estimators. The spec's version borrows from a fit that
contains the liability data; this one borrows from traits that may share no signal with them. Do not
write the spec's name on it.

Second, it is a spec **non-goal**. Spec line 52: "Binary / categorical liability lambda ... we don't
separately tune the liability prior." D-278 overrides decision 1 (the default), not the non-goals.

Third, and decisive: **it breaks committed tests that exist precisely to protect this behaviour.**
`tests/testthat/test-lambda-per-type.R:109` asserts
`expect_equal(bl_est$mu[, migr_col], bl_fixed$mu[, migr_col])` and the same for the categorical
columns and for `se` (lines 124-128); line 139 repeats it for `lambda_mode = "bayes"`. Those tests
were written to lock in the August `arc/lambda-per-type` fix, whose NEWS entry (`NEWS.md:100-123`)
records that the previous attempt to let lambda_mode touch the discrete path cost 19 pp of
Trophic.Level accuracy. The hybrid inverts all of it. The plan's gate G5 ("discrete columns' mu
**differ** between fixed_1 and estimate") is the direct negation of a committed assertion, and gate G8
(`devtools::test()` FAIL 0) would then be satisfied by rewriting the guard to agree with the new
behaviour. That is a guard deletion wearing a green gate, and a PR reviewer should reject it.

Recommendation, in order of preference:
1. **Cut it from this lane.** Keep discrete liability columns at lambda = 1, exactly as today.
   D-278's motivation is the continuous weak-signal case (sla; the four-arm sim's continuous floor);
   nothing in the trigger evidence is about binary or categorical traits. This makes the discrete
   path bit-identical for free, keeps `test-lambda-per-type.R` green as written, and removes the only
   part of the design that contradicts the approved spec.
2. If it is kept, make it the spec's version, not a proxy: compute lambda_bar from the summed NLLs
   over the **full** liability block, liability columns included (they are continuous pseudo-data
   after the E-step, so `build_pagel_nll_cache` accepts them unchanged), put it behind an explicit
   opt-in argument that defaults off, and replace G5 with a quality gate (below).

### C. "lambda_mode = fixed_1 stays bit-identical", and the 1e-8 oracle. WEAK.

Achievability is fine. `transform_tree_pagel(tree, 1.0)` is the exact identity in IEEE arithmetic
(`1.0 * x == x`, and `t + (1 - 1) * h_p == t + 0 == t`), so even an implementation that always routes
through the transform returns the same doubles; and the dispatcher's fixed_1 branches
(`fit_baseline.R:195-201`, `:422`) are untouched by the planned edits. The claim is reachable.

The **oracle is wrong**, in two ways. (1) 1e-8 is not bit-identity. It is roughly eight orders of
magnitude looser than double precision and would pass a changed numerical route (a re-ordered
Cholesky, an added nugget, a rebuilt tree object) while the note and the plan say "bit-identical" in
three places (note lines 93-94; plan lines 115, 128). Either assert `identical()` / tolerance 0 and
keep the word, or keep 1e-8 and stop calling it bit-identical. (2) More seriously, gate G1 is
`test_file("tests/testthat/test-joint-lambda.R")` EXPECT FAIL 0, on a test file the same agent writes
in the same slice, with **no stated pre-change reference**. A self-consistent run cannot detect a
regression. Before S2 touches `R/joint_mvn_solver.R`, generate and commit an `.rds` snapshot from
`origin/main` (ab02e31) of `fit_baseline(..., lambda_mode = "fixed_1")$mu`/`$se` on a fixed mixed-type
fixture and of `fit_mvn_bm_inhouse()$anc_recon`/`$anc_var`/`$pars$phylocov`, and assert against that
file. `test-lambda-per-type.R:69` already does the right thing for the dispatcher; the new file needs
the same discipline at the solver.

### D. Prior-work sweep receipt. WEAK.

Not vacuous: three of the four lines cite an actual command with an actual result (`git status -sb;
git stash list; branch_drift_check.sh` with the clean 0/0 verdict and `.gitignore:60-61`;
`search_notes(...)` plus `grep -in lambda memory/AGENT_LOG.md` plus two zero-hit files). The
twin-repo line asserts a negative with no command, and the external-prior-art line is a reasoned
exemption rather than a query, which is acceptable for a no-novelty lane.

The failure is coverage, not form. The sweep looked at git and at the brain and never swept **this
repo's own prior work on this exact feature**. Missing from the receipt, all of it decision-relevant:
`tests/testthat/test-lambda-per-type.R` (the committed guard the hybrid breaks, section B iii),
`NEWS.md:100-123` (the 19 pp Trophic.Level regression from the last time lambda touched the discrete
path), `specs/2026-05-18-cv-lambda-selection-design.md` and
`specs/2026-05-18-pagel-lambda-eigendecomp-speedup-design.md` (two sibling approved specs, one of
which owns the cache this lane builds on), and `bayes_lambda_for_col` / `cv_lambda_for_col`
(`pagel_lambda.R:58,113`), whose existence is what makes plan Question 2 a live question. A sweep that
found the eigendecomp cache but not the test file that forbids the design is a sweep that stopped at
`git grep lambda` in `R/`. Re-run it over `tests/`, `specs/` and `NEWS.md` and record the hits.

### E. What the gate list would let through.

1. **G5 is a difference test, not a quality test.** "discrete columns' mu differ between fixed_1 and
   estimate" passes on any change, including a change that costs another 19 pp. Replace with: on
   AVONET300, categorical accuracy and binary accuracy under `estimate` are >= their `fixed_1` values
   minus a stated tolerance, plus the same on the lambda = 0.3 DGP. If recommendation B(iii).1 is
   taken, G5 becomes its opposite (identical, tolerance 0) and the existing test file already is the
   gate.
2. **G1 cannot fail** as specified. See C: no pre-change reference.
3. **G2n is good and should be kept**, but is mis-worded. "a test that ... FAILS" would surface as
   FAIL 1 and contradict G1/G8's FAIL 0. Express it as an assertion that the estimator separates the
   hypotheses, e.g. `expect_false(abs(lambda_hat_on_lambda1_DGP - 0.3) < 0.05)`, so the suite stays
   green while the oracle demonstrably discriminates.
4. **G6 is a spelling test that contradicts the plan's own Question 2.** It greps for
   `lambda_mode = c("estimate", "fixed_1"` and expects 4 hits. The live enums are
   `c("fixed_1", "estimate", "cv", "bayes")` at `fit_baseline.R:171`, `fit_pigauto.R:339`,
   `impute.R:333`, and `R/multi_impute.R` has no `lambda_mode` argument at all. Satisfying G6 means
   reordering to a two-element enum and dropping "cv" and "bayes", which breaks the `switch` at
   `fit_baseline.R:195-201` and `test-lambda-per-type.R:139`, while plan Question 2 says both modes
   stay. Replace the grep with a behavioural check: `formals(impute)$lambda_mode[[2]]` evaluates to
   `"estimate"` at each of the four entry points. Also flip the hard-coded fallback at
   `multi_impute_trees.R:559` (`baseline_arg("lambda_mode", "fixed_1")`), or multi-tree fits stay at
   lambda = 1 under the new default and nothing in the ledger would notice.
5. **G3 has no replication and a lopsided band.** `lambda_hat in [0.15, 0.5]` on a lambda = 0.3 DGP,
   apparently one seed, against G2's 20 seeds within 0.05. State seeds and reps, and make the band
   symmetric or justify the asymmetry.
6. **G12's coverage band is two-sided on a one-sided guarantee.** Split conformal guarantees
   coverage >= 0.95; a run at 0.97 is over-coverage, not a defect, yet the 0.95-0.96 band marks it
   failed. Use >= 0.94 and report the observed value. And "target 'most'" sitting next to a hard
   number (<= 0.972) is unfalsifiable; keep the number, drop the word.
7. **G10's CHECK is a token the script prints about itself.** `AVONET_OK` is only as good as the
   script that emits it, and a two-sided "< 5%" passes a 4.9% regression. Record the two RMSE numbers
   in the gate EVIDENCE line, not just the token.
8. **Two spec-named risks have no gate at all.** Spec 4.5 says the stored lambda must let
   `predict.pigauto_fit()` rebuild the same baseline: nothing asserts that the prediction-time
   baseline equals the fit-time one under the stored lambda. Spec 6.3 says to confirm
   `multi_impute_trees()` re-estimates lambda per tree rather than caching the MCC tree's: nothing
   checks it, and `multi_impute_trees.R:559` is exactly where that would go wrong. Add both.
9. **No runtime gate.** Per-trait lambda adds K one-dimensional optimisations, K Henderson builds
   (O(n^2) each, see A) and, on the dense per-column route, K eigendecompositions at O(n_o^3). S8b is
   estimated at 3-4 h with a 5 h stop rule, but that estimate is carried over from a run at lambda = 1.
   Have S8a record seconds-per-fit and compare against the committed run before S8b is launched; a
   1.6x slowdown eats the whole margin.
10. **Spec 6.2's callsite audit is not carried into any gate.** S3 adds a `lambda` argument to
    `bm_impute_col_with_cov`, whose only callsite (`fit_baseline.R:747`) passes its first three
    arguments positionally. Adding the parameter anywhere but after `lrt_threshold` breaks it
    silently. Cheap gate: the callsite passes `lambda` by name.

11. **`docs/` is git-ignored in this repo** (`.gitignore:44`), so this note, the planned
    `docs/dev-log/lambda-default/benchmark.md` and `prerun.md`, and the S12 after-task report under
    `docs/dev-log/after-task/` are all invisible to git and will not appear in the draft PR. G14
    (`tools/check-after-task.R` passes) can therefore go green while the PR carries no reviewable
    evidence at all, which is the opposite of what the repo's after-task discipline is for. Decide
    before S12 whether these artifacts move to a tracked path or get force-added.

### Required before S2 is dispatched

- Decide B(iii): cut the discrete-liability inherit, or re-scope it to the spec's version behind an
  opt-in that defaults off. Either way, say in this note what happens to
  `tests/testthat/test-lambda-per-type.R:109` and `:139`.
- Pin the lambda_k estimator and the default predictor to the same mean model (B i), and gate it.
- Commit the pre-change `.rds` reference from origin/main before any solver edit, and either assert
  tolerance 0 or stop saying "bit-identical" (C).
- Rewrite G5, G6, G2n's wording, G12's coverage bound; add the predict-rebuild, per-tree-lambda,
  runtime and callsite gates (E).
- Re-run the prior-work sweep over `tests/`, `specs/` and `NEWS.md` and append the hits to the plan's
  receipt (D).

Nothing above touches the core of the lane. The solver lambda, the per-trait lambda_k on
continuous-family columns, the covariate GLS-in-profile path, and the default flip are all sound and
worth building.

## 7. Decision after review (Ada, 2026-09-22, acting on Rose's verdict)

- **B(iii) cut.** Discrete liability columns (binary, zi gate, ordinal liability, OVR synthetic columns)
  stay at lambda = 1, exactly as today. The committed guard `tests/testthat/test-lambda-per-type.R`
  (lines 24-27 and the bayes twin) stays unchanged and becomes the gate. Section 4's "free upgrade"
  sentence is withdrawn: the spec's decision 4 describes a lambda estimated over the whole liability
  block, which this lane does not build. The joint solver's `lambda_cols` argument names the columns that
  estimate lambda; every other column is fit at lambda = 1.
- **B(i) fixed.** The per-trait estimator and the default predictor share one mean model: the GLS
  phylogenetic mean mu_k(lambda_k) from the profile cache is subtracted from the observed cells before
  the Henderson solve on T(lambda_k) and added back to the predictions. Gate G2b requires the solver's
  per-column prediction to match `bm_impute_col(y, R, lambda = lambda_hat)` within 1e-3 on a DGP with a
  non-zero observed-subset mean.
- **C fixed.** The fixed_1 reference is `tests/testthat/fixtures/lambda_fixed1_reference_ab02e31.rds`,
  generated from a `git archive origin/main` copy before any solver edit (solver anc_recon / anc_var /
  phylocov on a K = 3 fixture; dispatcher mu / se / path on a 3-continuous + binary + categorical
  fixture). Tolerance 1e-12. The phrase "bit-identical" is withdrawn in favour of "identical within 1e-12".
- **A nit accepted.** `build_henderson_S_inv` is O(n^2) in time and memory because it forms `ape::vcv`;
  section 3's "O(n)" refers to the solve, not the build.
- **B(ii) disclosure.** `predict_method = "exact"` and `joint_refine_iter > 0` use `lambda_block`, not
  the per-trait values; the roxygen for those arguments says so (S5).
- **D re-swept** over `tests/`, `specs/`, `NEWS.md`: 16 test files mention lambda, three sibling lambda
  specs exist (baseline, cv-selection, eigendecomp-speedup), NEWS lines 100-123 record the 19 pp
  Trophic.Level loss when lambda last touched the discrete path. The decision above is consistent with all
  of them.
- **E adopted** in the ledger: G5 is the unchanged per-type guard; G5b callsite by name; G6 behavioural
  (formals) plus the `multi_impute_trees` fallback flip; G7b predict-rebuild; G7c per-tree lambda; G10
  records the numbers; G11 records seconds-per-fit before S8b; G12 coverage one-sided >= 0.94 and the
  word "most" dropped; G15 force-adds the evidence documents past the `docs/` ignore rule.
