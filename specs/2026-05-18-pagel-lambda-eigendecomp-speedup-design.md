# Pagel-λ NLL via eigendecomposition cache — performance spec

**Status:** draft for review.
**Depends on:** v0.10's `ml_lambda_for_col()` (current naive O(n³) per evaluation).
**Required by:** `2026-05-18-cv-lambda-selection-design.md` — CV-based λ at production scale is unusable without this speedup.
**Target version:** v0.11.0 (land together with CV spec).
**Branch (target):** `feature/pagel-lambda-perf` off main, merged before CV.

---

## 1. Why this exists

The current ML loop in `ml_lambda_for_col()` (and the proposed CV
loop in `cv_lambda_for_col()`) evaluates the profile NLL across a
λ grid. Each evaluation does:

1. Form `R(λ)_OO = λ * R_OO + (1-λ) * I` (cheap O(n²) memory write)
2. Cholesky-factor: O(n³) FLOPs

At n_obs = 1500 (BIEN bench scale), each Cholesky is ~3×10⁹ FLOPs
≈ 0.3 sec on a laptop. The ML grid + optim does ~16 evaluations →
~5 sec per trait. CV does 20 × 5 = 100 evaluations → ~30 sec per
trait. With 5 traits × 20 imputations = 100 traits per fit, that's
**50 minutes of overhead per CV fit at BIEN scale** — unusable for
the bench, and prohibitive for the user.

This is an **algorithmic** issue, not a language issue. R already
calls LAPACK (Fortran) for the Cholesky — moving to C++ would call
the same LAPACK and gain ≤ 2× from removing R-call overhead. The
real fix is to avoid redoing the O(n³) work at each λ.

## 2. Goal

Cut the per-λ-evaluation cost from O(n³) to O(n²) in pure R via
eigendecomposition caching. Total speedup: **~600× at n = 1500**,
making CV's 100-evaluation grid feasible.

Hard target: BIEN-scale (n_obs ≈ 1500) ML + CV combined complete
in < 10 sec per trait. AVONET-scale (n_obs ≈ 200) in < 0.5 sec
per trait.

## 3. Math

For symmetric R_OO, eigendecompose once:

$$\mathbf{R}_{OO} = \mathbf{U}\, \boldsymbol{\Lambda}\, \mathbf{U}^T$$

where U is orthogonal and Λ is diagonal with entries `λ_i ≥ 0`. Then

$$\mathbf{R}(\lambda)_{OO} = \lambda \mathbf{U}\boldsymbol{\Lambda}\mathbf{U}^T + (1-\lambda)\mathbf{I} = \mathbf{U}(\lambda \boldsymbol{\Lambda} + (1-\lambda)\mathbf{I})\mathbf{U}^T$$

The middle term `λΛ + (1-λ)I` is diagonal with entries
`d_i(λ) = λ λ_i + (1-λ)`. So:

- `log |R(λ)_OO| = sum_i log(d_i(λ))` — **O(n)** per λ
- `R(λ)_OO^{-1} = U diag(1/d_i(λ)) U^T`
- For any vector b: `R(λ)_OO^{-1} b = U diag(1/d_i(λ)) U^T b`
  - Compute `c = U^T b` once: **O(n²)**
  - `R(λ)^{-1} b = U (c / d(λ))` per λ: **O(n²)**

The profile NLL evaluation needs:
- `mu_hat(λ) = sum(R^{-1} y) / sum(R^{-1} 1)` — uses `c_y = U^T y`, `c_1 = U^T 1` cached, evaluates as `sum(c_y / d) / sum(c_1 / d)`: **O(n)** per λ
- `sigma2(λ) = e^T R^{-1} e / (n-1)` — uses `c_e = U^T e(λ)`: **O(n²)** per λ (e depends on mu_hat which depends on λ)

Total per λ: O(n²) with the cache, vs O(n³) without.

## 4. Architecture

### 4.1 Refactor `ml_lambda_for_col()`

Move the inner Cholesky out of the per-λ loop. Replace with:

```r
ml_lambda_for_col <- function(y, R, nugget = 1e-6,
                                lambda_grid = c(0.005, 0.995)) {
  obs <- which(!is.na(y))
  n_o <- length(obs)
  if (n_o < 10L) return(1.0)
  R_oo <- R[obs, obs, drop = FALSE]
  y_o  <- y[obs]
  ones <- rep(1, n_o)

  # ---- ONE-TIME eigendecomp -------------------------------------------
  eig <- eigen(R_oo, symmetric = TRUE)
  U <- eig$vectors           # n_o x n_o
  evals <- pmax(eig$values, 0)  # clip tiny negatives from numerical noise

  # Pre-compute U^T y and U^T 1
  c_y <- as.numeric(crossprod(U, y_o))
  c_1 <- as.numeric(crossprod(U, ones))

  # ---- O(n^2) NLL via cached eigendecomp ------------------------------
  nll <- function(lambda) {
    d <- lambda * evals + (1 - lambda) + nugget
    if (any(d <= 0)) return(.Machine$double.xmax)
    inv_d <- 1 / d
    # mu_hat = (c_y . inv_d) / (c_1 . inv_d)  -- O(n)
    num <- sum(c_y * inv_d)
    den <- sum(c_1 * inv_d)
    mu_hat <- num / den
    # e = y - mu_hat; c_e = U^T e = c_y - mu_hat * c_1
    c_e <- c_y - mu_hat * c_1
    # sigma2 = sum((c_e^2) * inv_d) / (n_o - 1)
    sigma2 <- sum(c_e * c_e * inv_d) / max(n_o - 1L, 1L)
    if (!is.finite(sigma2) || sigma2 <= 0) return(.Machine$double.xmax)
    log_det <- sum(log(d))
    0.5 * ((n_o - 1L) * log(sigma2) + log_det)
  }

  # Coarse grid + local optim as before
  grid <- seq(lambda_grid[1L], lambda_grid[2L], length.out = 11L)
  grid_nll <- vapply(grid, nll, numeric(1L))
  ...
}
```

The eigendecomp is O(n_o³) ONE-TIME (same cost as one Cholesky).
After that, every NLL eval is O(n_o). 11-point grid + 5 optim
refinements + (eventually) 100-point CV grid all share the same cache.

### 4.2 Same refactor for `bm_impute_col(..., lambda = numeric)`

When called with a fixed numeric λ, the prediction also benefits
from the eigendecomp cache IF we're calling it inside an ML/CV
loop. Refactor to expose `bm_impute_col_cached()` that takes
pre-computed `c_y`, `c_1`, `evals`, `U`. The public
`bm_impute_col()` builds the cache once and delegates.

For the prediction (not just NLL), we need:
- mu_pred at missing cells: requires R_mo, the off-diagonal block.
  At eigenbasis: R_mo = R[miss, obs] doesn't directly transform.
  Workaround: form `c_pred = U^T R_mo^T` ONCE outside the loop (O(n²)),
  then per λ: `pred = mu_hat + R_mo (U diag(inv_d) U^T) e
                  = mu_hat + (c_pred^T diag(inv_d)) c_e`
  Per λ: O(n_m * n_o) = O(n²). Same complexity as NLL.

### 4.3 Where the cache lives

Cache is built per `bm_impute_col()` call. For ML/CV which calls
`bm_impute_col()` many times across a λ grid:

- Option A: refactor `bm_impute_col` to expose a "build cache, then
  eval at λ" two-stage API.
- Option B: have `ml_lambda_for_col` and `cv_lambda_for_col` build
  their own cache inline (don't share with `bm_impute_col`).

Option B is simpler. The cache is internal to the ML/CV functions;
`bm_impute_col` stays the same single-call API.

## 5. Testing

### 5.1 Correctness regression tests

For 100 random (y, R, λ) inputs, the eigendecomp-cached NLL must
match the dense-Cholesky NLL to within 1e-8 absolute. Same for
`bm_impute_col_cached` predictions vs `bm_impute_col` predictions.

```r
test_that("[pagel-perf] cached NLL matches dense Cholesky NLL to 1e-8", {
  for (seed in 1:100) {
    set.seed(seed)
    tree <- ape::rcoal(rpois(1, lambda = 60) + 20)
    R <- stats::cov2cor(ape::vcv.phylo(tree))
    R <- R[tree$tip.label, tree$tip.label]
    n <- nrow(R)
    y <- rnorm(n); y[sample(n, n %/% 4)] <- NA
    for (lam in c(0.05, 0.5, 0.95)) {
      ref <- nll_dense_cholesky(y, R, lambda = lam)
      out <- nll_eigendecomp_cached(y, R, lambda = lam)
      expect_equal(out, ref, tolerance = 1e-8)
    }
  }
})
```

### 5.2 Speed regression tests

At n=1500, the ML loop must complete in < 1 sec wall-time. Use
`microbenchmark` against a reference implementation.

### 5.3 Existing tests still pass

All v0.10 `test-pagel-lambda.R` tests must still pass — the cache
is an internal speedup, not a behaviour change.

## 6. Risks

### 6.1 Eigendecomp of large dense matrices

`eigen()` on a 1500x1500 dense matrix is ~1-2 sec on a laptop. Same
ballpark as one Cholesky. Fine.

For n > 10000, eigendecomp gets expensive. At that scale, the Hadfield
sparse precision path (already in `R/joint_mvn_solver.R`) is the
right approach — different code path. Document that the eigendecomp
cache is appropriate for n ≤ 10000.

### 6.2 Numerical accuracy of eigendecomp vs Cholesky

`eigen()` returns eigenvalues with O(eps * cond(R)) error. For
well-conditioned R (most phylogenetic correlation matrices), this is
fine. For near-singular R, both eigendecomp and Cholesky get
unstable.

Mitigation: clip eigenvalues to >= 1e-12; warn if condition number
> 1e10.

### 6.3 Code duplication risk

We have two NLL implementations in flight (dense Cholesky reference
in `ml_lambda_for_col` v0.10, eigendecomp-cached version in v0.11).
Keep the v0.10 dense version as a private fallback for the n_o < 50
edge case where eigendecomp isn't worth it.

## 7. Rollout

### 7.1 v0.11.0 (ships with CV spec)

- Refactor `ml_lambda_for_col` to use eigendecomp cache
- Same refactor for `cv_lambda_for_col` (new)
- Correctness regression tests
- Speed regression tests
- NEWS entry: "X× speedup on lambda_mode='estimate' and 'cv'"

## 8. Open questions

1. **Cache the eigendecomp at the `fit_baseline` level?** If multiple
   traits share the same R_OO (because their observed-cell sets are
   identical), one eigendecomp could serve all of them. Trade-off:
   code complexity. Probably not worth it — each trait usually has
   a different observed set.
2. **Use sparse SVD for n > 5000?** Standard `eigen()` is dense.
   `RSpectra::eigs_sym()` is sparse. Out of scope for v0.11.

## 9. Success criteria

- Hard: ML loop at n=1500 in < 1 sec (was ~5 sec)
- Hard: CV loop at n=1500 in < 10 sec (would have been ~30 sec)
- Correctness: cached NLL matches dense NLL to 1e-8
- All v0.10 tests still pass
- R CMD check clean
