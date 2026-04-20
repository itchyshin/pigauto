# Phase 7 EM: off-diagonal Σ conditioning in the threshold E-step

**Status:** Draft (awaiting Phase 6 merge)
**Author:** Claude Opus 4.7 (working with Shinichi Nakagawa)
**Date:** 2026-04-19

---

## Problem

Phase 6 EM (v0.9.1.9000) feeds the **diagonal** of the phylopars Σ back as
per-trait prior SD on each iteration. The off-diagonal entries — which
encode the evolutionary covariance between different discrete-trait
liabilities — are ignored:

```
# Phase 6 today (diagonal only)
sd_prior_j = sqrt(Sigma_k[j, j])
mu_prior_j = 0

# Phase 7 target (off-diagonal conditioning)
mu_prior_j = Sigma_k[j, -j] %*% solve(Sigma_k[-j, -j]) %*% L_k[-j]   # per-cell
sd_prior_j = sqrt(Sigma_k[j, j] - Sigma_k[j, -j] %*% solve(Sigma_k[-j, -j]) %*% Sigma_k[-j, j])
```

Where `L_k[-j]` is the current posterior mean liability on the other
traits for the same species. When two binary traits are correlated on
the liability scale (e.g., `Σ[1, 2] = 0.6`), observing trait 1 should
shift trait 2's prior mean AWAY from zero, not just tighten its SD.

Phase 6's diagonal-only prior leaves that signal on the table. Phase 7
captures it.

## Goal

Add a second opt-in flag (`em_offdiag = FALSE` default) that, when
combined with `em_iterations >= 2L`, uses the full conditional-MVN prior
for each liability cell at iteration k+1. Default behaviour unchanged;
Phase 6's diagonal-only path remains the primary opt-in.

## Non-goals

- **Changing Phase 6 defaults.** `em_iterations = 0L` still disables EM;
  `em_offdiag = FALSE` still uses the Phase 6 diagonal path.
- **New trait types.** Same scope as Phase 6 (binary + ordinal + OVR
  categorical).
- **Σ full inverse caching.** Per-cell conditional-MVN computes a
  `(K-1) × (K-1)` inverse per liability column; we compute those K
  inverses once per iteration, not per cell.
- **Iterating OVR categorical's off-diagonals.** Phase 7 keeps OVR
  categorical on a length-K diagonal σ_cat (Phase 6 behaviour). Going
  off-diagonal on OVR means the K classes of one trait would condition
  each other, but OVR already decorrelates via the "is_class_k vs rest"
  one-hot — off-diagonal would re-couple and reintroduce the rank-(K−1)
  instability that OVR was built to dodge. Binary + ordinal get the
  Phase 7 upgrade; OVR stays Phase 6.

## Design decisions

### 1. Scope: binary + ordinal (off-diagonal); OVR stays Phase 6

The Phase 7 conditional-MVN prior fires for binary and ordinal liability
columns in the threshold-joint baseline. OVR categorical is unchanged —
still uses the Phase 6 diagonal σ_cat loop.

### 2. Per-cell conditional-MVN prior

For each observation cell `(i, j)` where `j` is a binary/ordinal liability
column at iteration `k + 1`:

```
# Inputs at iter k+1
Σ     := Σ_k          # K x K, from previous phylopars fit
L_hat := L_k           # n x K, posterior liability means from previous fit

# Per-cell prior for col j, species i
mu_i_j   <- t(Σ[j, -j]) %*% solve(Σ[-j, -j]) %*% L_hat[i, -j]
var_i_j  <- Σ[j, j] - t(Σ[j, -j]) %*% solve(Σ[-j, -j]) %*% Σ[-j, j]
sd_i_j   <- sqrt(max(var_i_j, eps))
```

Then the per-cell E-step is called with `(mu_prior = mu_i_j, sd_prior = sd_i_j)`
instead of Phase 6's `(mu_prior = 0, sd_prior = sqrt(Σ[j, j]))`.

### 3. Cost management

Naive per-cell `solve(Σ[-j, -j])` is K cell-level matrix inversions. To
keep cost reasonable we precompute the `K` (K−1)×(K−1) conditional
matrices once per iteration and reuse them across all `n × K` cells.
Plus the **missing-pattern quirk**: `L_hat[i, -j]` may have NAs when
species `i` has no observed value on some other trait. For those cells
we restrict the condition to the observed subset of `-j` — a further
inner loop per cell, but K is typically small (2–10 for real datasets)
so this is tractable.

### 4. New argument + API

```r
fit_baseline(..., em_iterations = 0L, em_tol = 1e-3, em_offdiag = FALSE)
impute(..., em_iterations = 0L, em_tol = 1e-3, em_offdiag = FALSE)
```

- `em_offdiag = FALSE` (default): Phase 6 diagonal-only path (unchanged).
- `em_offdiag = TRUE` AND `em_iterations >= 2L`: Phase 7 full conditional prior.
- `em_offdiag = TRUE` AND `em_iterations <= 1L`: no effect (there's no
  previous Σ to condition on at iter 1). Emit a one-line message at
  verbose level.

### 5. Iteration control: reuse Phase 6's

- `em_iterations` max iters (same semantics).
- `em_tol` convergence check on relative Frobenius of **full** Σ (not
  just diag) since off-diagonals now matter.
- Same phylopars-failure fallback: revert to previous iter's Σ and stop.

### 6. Single exit point: converged Σ → one final phylopars fit

Same as Phase 6: after iteration stops, run phylopars one more time with
the converged Σ's conditional-MVN prior applied to each cell, producing
the downstream baseline.

### 7. Result object: reuse `em_state`

`em_state` gains one more field:

```r
em_state = list(
  iterations_run = integer(1),
  converged      = logical(1),
  final_delta    = numeric(1),
  em_offdiag     = logical(1)     # NEW: records whether Phase 7 fired
)
```

OVR's `em_state` stays as in Phase 6 (adds `sigma_cat` but no `em_offdiag`
— OVR doesn't support off-diagonal in Phase 7).

## Implementation sketch

New internal helper in `R/joint_threshold_baseline.R`:

```r
# Given Σ (K x K) and current L_hat (n x K posterior means), build a
# per-(i, j) list of (mu_prior, sd_prior) scalars for every cell.
# Returns an n x K matrix of mu_prior and an n x K matrix of sd_prior.
#' @keywords internal
#' @noRd
build_conditional_prior <- function(Sigma, L_hat, eps = 1e-8) {
  K <- ncol(Sigma)
  n <- nrow(L_hat)
  mu_mat <- matrix(0,  nrow = n, ncol = K)
  sd_mat <- matrix(1,  nrow = n, ncol = K)

  for (j in seq_len(K)) {
    sigma_jj  <- Sigma[j, j]
    sigma_jnj <- Sigma[j, -j, drop = FALSE]       # 1 x (K-1)
    sigma_njn <- Sigma[-j, -j, drop = FALSE]       # (K-1) x (K-1)

    # Σ_jj|minus_j  (scalar)
    inv_njn  <- solve(sigma_njn + diag(eps, K - 1))
    var_j    <- sigma_jj - as.numeric(sigma_jnj %*% inv_njn %*% t(sigma_jnj))
    sd_mat[, j] <- sqrt(max(var_j, eps))

    # mu_j|minus_j  for each species: σ_jnj %*% inv_njn %*% L_hat[i, -j]
    # Some L_hat[i, -j] entries may be NA — fall back to the unconditional
    # prior (mu = 0, sd = sqrt(Σ_jj)) for those cells.
    for (i in seq_len(n)) {
      li <- L_hat[i, -j]
      na <- is.na(li)
      if (all(na)) {
        mu_mat[i, j] <- 0
        sd_mat[i, j] <- sqrt(sigma_jj)
        next
      }
      if (any(na)) {
        # Restrict to observed subset
        obs    <- which(!na)
        inv_obs <- solve(sigma_njn[obs, obs, drop = FALSE] +
                          diag(eps, length(obs)))
        mu_mat[i, j] <- as.numeric(
          sigma_jnj[, obs, drop = FALSE] %*% inv_obs %*% li[obs])
        sd_mat[i, j] <- sqrt(max(
          sigma_jj - as.numeric(sigma_jnj[, obs, drop = FALSE] %*%
                                   inv_obs %*%
                                   t(sigma_jnj[, obs, drop = FALSE])),
          eps))
      } else {
        mu_mat[i, j] <- as.numeric(sigma_jnj %*% inv_njn %*% li)
      }
    }
  }

  list(mu_prior = mu_mat, sd_prior = sd_mat)
}
```

Wire into `fit_joint_threshold_baseline_em()`:

```r
if (em_offdiag && !is.null(prev_liab_means) && !is.null(prev_Sigma)) {
  # Phase 7: use conditional prior
  cp <- build_conditional_prior(prev_Sigma, prev_liab_means)
  built <- build_liability_matrix(data, splits = splits,
                                   soft_aggregate = soft_aggregate,
                                   mu_prior_mat = cp$mu_prior,
                                   sd_prior_mat = cp$sd_prior)
} else {
  # Phase 6: diagonal-only
  built <- build_liability_matrix(data, splits = splits,
                                   soft_aggregate = soft_aggregate,
                                   sd_prior_vec = sd_prior_vec)
}
```

This requires `build_liability_matrix()` to also accept per-cell
`mu_prior_mat` / `sd_prior_mat` (n × K matrices) as an alternative to the
length-K `sd_prior_vec`.

## Testing

1. **Unit tests** (`tests/testthat/test-phase7-em.R`, new):
   - `em_offdiag = FALSE` (default) byte-identical to Phase 6 at every
     `em_iterations` value (regression guard).
   - `em_offdiag = TRUE, em_iterations = 1L` → emits a message that
     off-diagonal has no effect at iter 1; falls through to Phase 6.
   - `em_offdiag = TRUE, em_iterations >= 2L` on synthetic correlated-
     binary (ρ = 0.6, n = 200) → `em_state$em_offdiag == TRUE`,
     converges, produces different liability posteriors than Phase 6.
   - Missing-pattern handling: cell with all-NA other traits falls back
     to unconditional prior.

2. **Smoke benchmark** (`script/bench_phase7_em.R`, new):
   - AVONET 300 real-data + synthetic correlated-binary sim.
   - Phase 6 diagonal vs Phase 7 off-diag, `em_iterations = 5L`.
   - Expected: Phase 7 gives a small additional lift over Phase 6 on
     synthetic correlated data; AVONET may or may not show lift (real
     correlation structure uncertain).

3. **Existing tests** stay green (`em_offdiag = FALSE` default).

## Documentation

1. Roxygen `@param em_offdiag` on `fit_baseline()` and `impute()`.
2. NEWS.md entry in the v0.9.1.9000 dev section (appended to the Phase 6
   block).
3. No new pkgdown article. Phase 7 is a refinement of Phase 6; users who
   care about EM enable both flags together.

## Backward compatibility

- `em_offdiag = FALSE` default → zero code-path change. Byte-identical to
  Phase 6.
- Saved `pigauto_fit` objects from pre-Phase-7 load unchanged.
- `fit_baseline()` / `impute()` signature gains one arg with default;
  not a breaking change.

## Out of scope (deferred further)

- **OVR categorical off-diagonal.** Would re-couple the K classes, undo
  the OVR rank-(K−1) fix. Stay Phase 6 on OVR.
- **Continuous-only joint_mvn_baseline**. Already uses full Σ via
  phylopars; nothing to iterate.
- **GNN-in-the-loop EM.** Same decision as Phase 6 — baseline only.

## Success criteria

- [ ] Regression: `em_offdiag = FALSE` byte-identical to Phase 6 at
      matching `em_iterations`.
- [ ] Convergence: `em_offdiag = TRUE, em_iterations = 5L` on synthetic
      correlated-binary (ρ = 0.6, n = 200) converges within 5 iters.
- [ ] `em_state$em_offdiag == TRUE` populated when Phase 7 fires.
- [ ] Test suite stays green.
- [ ] `R CMD check` stays 0 errors / 0 warnings / 1 note.
