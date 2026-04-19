# Phase 6 EM: iterated Σ-aware prior for the threshold-joint baseline

**Status:** Draft (awaiting user review)
**Author:** Claude Opus 4.7 (working with Shinichi Nakagawa)
**Date:** 2026-04-19

---

## Problem

pigauto's Level-C threshold-joint baseline (shipped in v0.9.0, `R/joint_threshold_baseline.R` + `R/ovr_categorical.R`) handles binary, ordinal, and categorical traits by:

1. Applying `estep_liability()` with a **plug-in prior `N(0, 1)`** to each discrete cell, converting observations into a continuous liability column.
2. Running `Rphylopars::phylopars()` on the combined liability + continuous matrix to produce posterior means and `Σ` (cross-trait covariance).
3. Decoding the liability posteriors back to the trait's native type.

The plug-in prior treats each discrete trait's liability as independent at step (1). But real datasets have **cross-trait correlation on the liability scale** — e.g., body size and beak depth are correlated, AND so are discrete categories that depend on them (foraging guild, diet class). The plug-in ignores that correlation precisely at the step where it would help imputation quality.

This mismatch is the Phase 6 EM target.

## Goal

Add an **opt-in iteration loop** around the threshold-joint path that feeds the Σ learned by `phylopars()` in iteration `k` back as the per-trait prior SD for iteration `k + 1`. Under correlated-discrete data, the prior tightens around the correct liability and discrete-trait imputation accuracy improves. Under independent-discrete data, the iteration converges at the plug-in answer in one step.

## Non-goals

- **Full off-diagonal conditioning** (use Σ's off-diagonals to condition each trait's liability on observed others). Deferred to a later Phase 7 — doubles the math complexity and is not justified until we validate diagonal-only EM in the field.
- **Changing v0.9.0 defaults.** The new behaviour is opt-in via `em_iterations = 0L` default. Default-on comes in a later release after field validation.
- **New trait types.** Phase 6 uses the existing liability encoders in `R/liability.R`.
- **Re-running the GNN per iter.** Only the baseline iterates. The GNN consumes the converged baseline once.

## Design decisions

### 1. Scope: binary + ordinal + OVR categorical ("Scope A")

All three discrete-trait paths iterate. OVR categorical is the hardest case, handled as §4 below. Continuous / count / proportion / zi_magnitude traits pass through Phase 6 unchanged — they're not threshold-encoded, so there's no liability E-step to iterate.

### 2. Iteration control

Fixed max `em_iterations = 5L` with early-stop on relative Frobenius norm:

```
delta = norm(Sigma_k - Sigma_{k-1}, "F") / norm(Sigma_{k-1}, "F")
if (delta < em_tol) break   # em_tol default = 1e-3
```

If `phylopars()` fails at any iteration (singular matrix, convergence failure), fall back to the previous iteration's Σ and stop. Emit a one-line warning. Do **not** error out — the Phase 6 path is strictly additive; regression to plug-in is an acceptable degenerate case.

### 3. How Σ enters the E-step: diagonal only

Each liability column's prior at iteration `k + 1` is:

```
sd_prior_j = sqrt(Sigma_k[j, j])
mu_prior_j = 0                 # symmetric around zero; threshold is the decision boundary
```

Off-diagonal `Σ[j, j']` is **not** used to condition `j` on `j'`. That's Phase 7.

Why diagonal-only suffices for the first iteration wave: `phylopars()`'s internal fit already captures cross-trait covariance in its own posterior — the liability column is just a better posterior mean at the cell level once `Σ[j, j]` reflects the trait's evolutionary variance. Diagonal-only is the minimum-diff change that exercises the feedback loop.

### 4. OVR categorical: per-trait diagonal variance vector across K classes

For a K-class categorical trait, pigauto runs K independent OVR binary fits (each "class k vs rest"). Phase 6 iterates **once per categorical trait**, with a length-K **variance vector** `σ_cat[k]` (diagonal only, consistent with §3) shared across the K fits:

```
σ_cat <- rep(1, K)                                    # plug-in start
for iter in 1..em_iterations:
  for each k in 1..K:
     base_k = fit_one_ovr_class(..., sd_prior = sqrt(σ_cat[k]))
  σ_cat_new[k] <- extract_scalar_variance(base_k)     # k = 1..K
  if (||σ_cat_new - σ_cat|| / ||σ_cat|| < em_tol) break
  σ_cat <- σ_cat_new
```

We deliberately avoid a K × K `Σ_cat` because obtaining its off-diagonals would require a joint phylopars fit on all K one-hot columns, which is rank-(K − 1) and fails with singular-matrix errors — the whole reason OVR exists. Stays consistent with §3.

Categorical traits do **not** share σ_cat across each other. Each categorical trait gets its own length-K iteration. Mixed-binary-and-categorical datasets have one shared binary/ordinal `diag(Σ)` plus per-categorical `σ_cat`s.

### 5. Single exit point: converged Σ → one final phylopars fit

After iteration stops, run `phylopars()` **one more time** with the converged Σ as prior, producing the posterior means + SEs that downstream (GNN, conformal) consume. The intermediate iteration fits are discarded. This keeps the post-Phase-6 plumbing identical to v0.9.0 — only the converged state propagates.

### 6. API

One new argument on `fit_baseline()` and `impute()`:

```r
fit_baseline(data, tree, splits = NULL, ..., em_iterations = 0L, em_tol = 1e-3)
impute(traits, tree, ..., em_iterations = 0L, em_tol = 1e-3)
```

- `em_iterations = 0L` (default): current behaviour. Zero code path change for existing users. No `em_state` attached.
- `em_iterations = 1L`: enters the EM wrapper but runs exactly one plug-in fit. Baseline output byte-identical to `0L`; `em_state` is attached for debugging (iterations_run = 1, converged = NA — no previous Σ to compare). Rarely useful in practice; documented for completeness.
- `em_iterations >= 2L`: the EM loop kicks in. Iter 1 uses the plug-in prior, iter 2+ uses Σ from the previous fit, up to `em_iterations` or until `em_tol` convergence. **This is the Phase 6 path.**
- `em_tol = 1e-3`: relative-Frobenius convergence threshold. Controls early-stop when `delta = ||Σ_k - Σ_{k-1}||_F / ||Σ_{k-1}||_F < em_tol`.

Both args are pass-through from `impute()` to `fit_baseline()`. `fit_pigauto()` does not gain the arg — EM only affects the baseline.

### 7. Result object: no change to `pigauto_fit` shape

The converged Σ is internal to the baseline fit. We store it in `fit$baseline$em_state` as:

```r
list(
  iterations_run  = integer(1),  # how many iters actually ran before convergence
  converged       = logical(1),  # TRUE if delta < em_tol, FALSE if hit max
  final_delta     = numeric(1),  # last Frobenius relative norm
  Sigma_path      = list()       # per-iter Σ snapshots for debugging (opt-in)
)
```

Not written to `pigauto_fit$model_config` (which is user-visible). Kept on the baseline object alone so downstream prediction code doesn't branch on whether EM ran.

## Implementation sketch

New internal wrapper in `R/joint_threshold_baseline.R`:

```r
# Iterate the plug-in E-step → phylopars → Σ → prior loop.
# Returns the final baseline (same shape as fit_joint_threshold_baseline()).
fit_joint_threshold_baseline_em <- function(data, tree, splits, graph,
                                             em_iterations, em_tol) {

  Sigma_prior <- diag(1, length(all_liability_cols))   # start at plug-in
  prev_base   <- NULL

  for (iter in seq_len(em_iterations)) {
    # Use Sigma_prior's diagonal as per-trait SD
    L_mat <- build_liability_matrix(data, splits,
                                     sd_prior = sqrt(diag(Sigma_prior)))
    base <- fit_joint_threshold_baseline_once(data, tree, splits, graph, L_mat)

    Sigma_new <- extract_covariance(base$phylopars_fit)  # K x K from phylopars
    delta <- frobenius_rel_diff(Sigma_new, Sigma_prior)

    if (!is.null(prev_base) && delta < em_tol) break
    Sigma_prior <- Sigma_new
    prev_base   <- base
  }

  base$em_state <- list(iterations_run = iter, converged = delta < em_tol,
                         final_delta = delta)
  base
}
```

Dispatch in `fit_baseline()`:

```r
if (em_iterations >= 1L && use_threshold_joint) {
  base <- fit_joint_threshold_baseline_em(data, tree, splits, graph,
                                           em_iterations, em_tol)
} else {
  base <- fit_joint_threshold_baseline(data, tree, splits, graph)  # current path
}
```

OVR categorical loop lives in `R/ovr_categorical.R`:

```r
fit_ovr_categorical_fits_em <- function(data, cat_trait, tree, splits, graph,
                                         em_iterations, em_tol) {
  K <- length(levels(data$traits[[cat_trait]]))
  Sigma_cat <- diag(1, K)

  for (iter in seq_len(em_iterations)) {
    per_k_fits <- lapply(seq_len(K), function(k) {
      fit_one_ovr_class(data, cat_trait, k, tree, splits, graph,
                         sd_prior = sqrt(Sigma_cat[k, k]))
    })
    Sigma_new <- assemble_ovr_covariance(per_k_fits)  # K x K
    if (frobenius_rel_diff(Sigma_new, Sigma_cat) < em_tol) break
    Sigma_cat <- Sigma_new
  }

  list(fits = per_k_fits, Sigma_cat = Sigma_cat,
       em_state = list(iterations_run = iter, converged = TRUE))
}
```

Helpers:

- `estep_liability_binary()` / `_ordinal()` / `_categorical()` in `R/liability.R` already accept `sd_prior`. No signature change.
- `build_liability_matrix()` in `R/joint_threshold_baseline.R` gains an optional `sd_prior` vector arg (length = number of liability columns). Default `NULL` preserves current `1` behaviour.
- `extract_covariance(phylopars_fit)` and `frobenius_rel_diff(A, B)` are new ≤10-line helpers.

## Testing

1. **Unit tests** (`tests/testthat/test-phase6-em.R`, new):
   - `em_iterations = 0L` → byte-identical to v0.9.0 (regression guard).
   - `em_iterations = 5L` on a tiny synthetic correlated-binary dataset → converges within 5 iters AND produces different Σ than the plug-in first iter.
   - `em_iterations = 5L, em_tol = 1.0` → stops at iter 1 (delta always < 1.0).
   - `phylopars()` failure mid-iter → falls back to previous Σ, emits warning, returns a valid baseline.
   - OVR categorical with K = 3 traits → `em_state$iterations_run` recorded per-trait.

2. **Smoke benchmark** (`script/bench_phase6_em.R`, new):
   - **Datasets:** AVONET 300 (real, mixed-type) + synthetic correlated-binary (n = 300, 4 traits with known Σ of rho = 0.6).
   - **Comparison:** `em_iterations = 0L` (baseline) vs `em_iterations = 5L`.
   - **Metrics:** per-trait accuracy, ECE (expected calibration error) for binary/categorical, wall time, iterations_run.
   - **Expected result:** synthetic correlated-binary should show +2–5pp accuracy; AVONET 300 should show ±2pp (real data is noisier, cross-trait correlation may be weak). Wall time should be `iterations_run × baseline wall`, typically 2–3× slower.

3. **Existing tests**: `tests/testthat/test-joint-threshold-baseline.R` and `test-ovr-categorical.R` stay green because `em_iterations = 0L` is the default — behaviour is byte-identical.

## Documentation

1. `R/joint_threshold_baseline.R` + `R/ovr_categorical.R` roxygen gain `@param em_iterations` and `@param em_tol` entries.
2. `R/fit_baseline.R` and `R/impute.R` roxygen gain the same.
3. `NEWS.md` entry in the v0.9.2 (or v0.10.0) dev section:
   > **Phase 6 EM for threshold-joint baseline.** New opt-in arg `em_iterations = 0L`. When ≥ 1, the threshold-joint path (binary + ordinal + OVR categorical) iterates: feeds `Σ` from the previous `phylopars()` fit back as per-trait prior SD, up to `em_iterations` times (default 5 when enabled), early-stopping when `‖Σ_k − Σ_{k−1}‖_F / ‖Σ_{k−1}‖_F < em_tol` (default 1e-3). Closes the final open item from the v0.9.0 "Deferred to future releases" list.
4. No new pkgdown article. The `em_iterations` knob is a power-user refinement; no tutorial needed for v0.9.2. If field usage shows it needs one, add later.
5. CLAUDE.md internal note: update the "Phase 6 EM will replace this with an iterated, Σ-aware prior" line under the threshold-joint baseline section to reflect that Phase 6 shipped.

## Backward compatibility

- `em_iterations = 0L` default → zero change for any existing user code. Byte-identical baseline output.
- Saved `pigauto_fit` objects from v0.9.0/v0.9.1 load unchanged; they have no `em_state` attribute, which is fine (it's optional metadata).
- `fit_baseline()` signature gains two args with defaults; not a breaking change.

## Out of scope (deferred to Phase 7)

- **Full off-diagonal conditioning.** Using `Σ_k[j, j']` to condition trait `j`'s liability on trait `j'`'s posterior mean during the E-step. Mathematically richer, roughly 2× the code, empirically unknown benefit. Revisit after Phase 6 field data.
- **EM for continuous-only Level-C (`R/joint_mvn_baseline.R`).** The continuous-only path already uses phylopars' full Σ — no plug-in prior to iterate on. No Phase 6 work needed there.
- **`em_iterations` default change from 0 to 5.** Shipped in v0.9.2 / v0.10.0 as opt-in. A later release can flip after validation data accumulates.
- **GNN-in-the-loop EM.** The GNN consumes the converged baseline once. No iteration of GNN training alongside baseline EM. Combinatorial blow-up risk.

## Success criteria

- [ ] `em_iterations = 5L` on a synthetic correlated-binary dataset (rho = 0.6, n = 300) converges within 5 iters and shows ≥ 2pp accuracy lift vs `em_iterations = 0L`.
- [ ] `em_iterations = 0L` produces byte-identical baseline output vs v0.9.1 (regression test passes).
- [ ] AVONET 300 smoke bench: `em_iterations = 5L` completes without error, within 3× the baseline wall time.
- [ ] `em_state` is populated on the baseline fit with `iterations_run`, `converged`, `final_delta`.
- [ ] Test suite stays green (778 → ~785 tests with the new Phase 6 block).
- [ ] `R CMD check` stays 0 errors / 0 warnings / 1 note.
