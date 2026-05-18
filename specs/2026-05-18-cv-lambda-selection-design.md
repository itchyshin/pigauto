# Cross-validation Pagel-λ selection — design spec

**Status:** draft for review.
**Depends on:** v0.10's `lambda_mode = "estimate"` (ML path) which is shipped but underperforms on weak-phylo-signal traits.
**Target version:** v0.11.0.
**Branch (target):** `feature/cv-lambda-selection` off main post-v0.10.0.

---

## 1. Why this exists

v0.10's `ml_lambda_for_col()` does point ML estimation of Pagel's λ. On the BIEN h2h CI bench it works correctly on 4 of 5 traits but **regresses sla**: ML picks λ̂ ≈ 0.005 because the profile log-likelihood is flat in the [0, 0.3] interval and finite-sample noise dominates, while RMSE-optimal λ on validation cells is ≈ 0.7–0.9.

Concretely:

| trait | ML λ̂ | val z-RMSE @ ML | val z-RMSE @ RMSE-opt λ | gap |
|-------|-------|-----------------|--------------------------|-----|
| height_m | 0.95 | 0.732 | 0.732 (λ=0.9) | 0 |
| leaf_area | 0.83 | 0.91 | 0.89 (λ=0.99) | 0.02 |
| sla | 0.005 | 0.98 | **0.89 (λ=0.9)** | **0.10** |
| seed_mass | 0.99 | 0.60 | 0.60 | 0 |
| wood_density | 0.79 | 0.79 | 0.77 (λ=0.99) | 0.02 |

The sla gap is structural: ML is the wrong criterion when the
likelihood is plateau-y. BACE's MCMCglmm wins on sla (z-RMSE 0.69)
because Bayesian model averaging integrates over the posterior; an
ML point estimate cannot replicate that with a single λ.

## 2. Goal

Add `lambda_mode = "cv"` to pigauto. Per-trait λ chosen by k-fold
cross-validation of imputation RMSE on the observed cells. Recovers
the RMSE-optimal interior λ even when ML collapses to a boundary.

Hard target: BIEN h2h `sla` z-RMSE drops from 0.98 (current ML) to
≤ 0.80, closing ≥ 80% of the gap. Other traits stay within 0.02 of
their current best. `BIEN_RMSE_TOTAL` drops to ≤ 6.85.

## 3. Architecture

### 3.1 New helper: `cv_lambda_for_col()`

```r
cv_lambda_for_col(y, R, nugget = 1e-6,
                    lambda_grid = seq(0.01, 0.99, by = 0.05),
                    k = 5L, seed = 1L)
```

Returns a scalar λ_hat in `[0, 1]`. Algorithm:

1. Split the observed cells of `y` into `k` folds (seeded).
2. For each candidate λ in `lambda_grid`:
   a. For each fold f: mask the fold's cells as NA; run
      `bm_impute_col(y_masked, R, lambda = λ)`; record squared error
      on the held-out cells.
   b. Sum squared errors across folds → cv_rmse(λ).
3. Return `argmin λ cv_rmse(λ)`.

Cost per trait: `length(lambda_grid) * k` ≈ 20 * 5 = 100 Cholesky
solves on `n_obs * (1 - 1/k)` cells. At n_obs = 1500, ~75 sec per
trait without the eigendecomp speedup (see sibling spec
`2026-05-18-pagel-lambda-eigendecomp-speedup-design.md`); ~0.1 sec
WITH it.

### 3.2 Dispatcher

```r
fit_baseline(..., lambda_mode = c("fixed_1", "estimate", "cv"))
fit_pigauto(..., lambda_mode = c("fixed_1", "estimate", "cv"))
```

`"cv"` uses `cv_lambda_for_col()`. Falls back to `"estimate"` (ML) if
n_obs < 50 per trait (folds become too small).

### 3.3 Store CV results on the fit

```r
model_config$lambda_method <- "cv"
model_config$lambda_per_trait <- named numeric  # per-trait λ̂
model_config$lambda_cv_curves <- list_of_data_frames  # for diagnostics
```

### 3.4 What CV does that ML can't

ML maximises observed-cell likelihood. CV maximises out-of-fold
prediction RMSE. They differ when:

- The likelihood surface has plateaus (weak signal)
- The training distribution doesn't match the test distribution
- Model misspecification matters more than ML can detect

Pagel's λ on weak-signal traits hits all three. CV is the standard
fix in the phylogenetic comparative methods literature (e.g.
`phylolm` recommends it for `model = "lambda"` when the fit looks
suspicious).

## 4. Non-goals

- Bayesian MAP via prior on λ. Equivalent power, more code, harder
  default-prior choices. Out of scope.
- Multi-trait joint CV. We could share folds across traits, but the
  benefit is marginal and code complexity grows.
- Adaptive grid refinement. The 20-point grid is fine.

## 5. Testing

### 5.1 Unit tests (`tests/testthat/test-pagel-lambda-cv.R`)

1. **CV with k=5 on pure-BM data returns λ near 1.** Sanity check.
2. **CV with k=5 on iid white-noise data returns λ near 0.** Sanity check.
3. **CV with k=5 on a hand-crafted "ML-boundary" trait returns the RMSE-optimal interior λ.** This is the structural test that distinguishes CV from ML — we construct data where ML picks 0 but CV picks ≈ 0.5.
4. **CV is reproducible given a seed.** Same data + same seed → same λ̂.
5. **CV with n_obs < 50 falls back to ML gracefully.** Edge case.

### 5.2 BIEN h2h acceptance bench

Same setup as v0.9.4 / v0.10. Compare `lambda_mode = "cv"` against
`lambda_mode = "fixed_1"`:

- sla z-RMSE drops from 0.91 → ≤ 0.80 (CV picks interior λ ≈ 0.7-0.9)
- height_m / leaf_area / seed_mass / wood_density stay within 0.02
- BIEN_RMSE_TOTAL drops 7.236 → ≤ 6.85

If the sla target isn't hit, the spec is wrong and we don't ship.

### 5.3 AVONET smoke regression

Strong-phylo-signal data. CV should pick λ near 1 for all 7
AVONET continuous traits. RMSE should stay within 5% of v0.10
fixed_1 baseline.

## 6. Risks

### 6.1 CV is slow at scale

20-point grid × 5 folds × O(n³) Cholesky per fold = 100× the cost
of a single fixed_1 fit. At n=2000, the v0.10 fit is ~10 min;
adding CV would push it to ~16 hours per BIEN fit — unusable.

**Mitigation:** the sibling spec `2026-05-18-pagel-lambda-eigendecomp-speedup-design.md` proposes ~600× speedup for the inner Cholesky loop via eigendecomposition caching. CV becomes ~1 min per trait at n=2000 with it. **This spec depends on that one.** Land both together in v0.11.

### 6.2 CV variance across folds

5 folds at n_obs ≈ 1500 means ~300 held-out cells per fold. CV RMSE
variance can be substantial. Solution: use `k = 10` if n_obs > 100;
fall back to `k = 5` otherwise.

### 6.3 CV picks 0 if data is genuinely iid

This is correct behaviour — for iid data the grand-mean predictor
IS optimal. CV will pick λ near 0 and the prediction collapses to
grand mean. The test in §5.1 verifies this.

### 6.4 CV depends on the seed for fold assignment

We seed deterministically (default seed = 1L). Different seeds give
slightly different λ̂. Document in roxygen.

## 7. Rollout

### 7.1 v0.11.0

- Add `cv_lambda_for_col()` in `R/pagel_lambda.R`
- Add `"cv"` to the `lambda_mode` enum
- Dispatcher wires CV through `fit_baseline` → `fit_pigauto` → `impute` → `multi_impute`
- Tests (§5.1)
- BIEN h2h acceptance bench (§5.2)
- AVONET smoke regression (§5.3)
- NEWS entry describing the option

### 7.2 v0.12 (later, if CV is the right default)

Consider flipping default from `"fixed_1"` to `"cv"`. Only after
v0.11 release cycle confirms CV is reliably better across the BACE
h2h bench.

## 8. Open questions

1. **CV k:** 5 or 10? 5 is faster, 10 is lower-variance. Default 5,
   user-configurable.
2. **Grid resolution:** `seq(0.01, 0.99, by = 0.05)` = 20 points.
   Coarser (10) is faster; finer (40) is more precise. 20 is the
   compromise.
3. **Should CV's pick be reported in the `pigauto_fit` summary?**
   Recommendation: yes, surface as `summary(fit)$lambda_per_trait`
   for diagnostics.

## 9. Success criteria

- Hard: BIEN h2h `sla` z-RMSE ≤ 0.80 (closes ≥ 80% of the gap)
- Soft: `BIEN_RMSE_TOTAL` ≤ 6.85
- Back-compat: `lambda_mode = "fixed_1"` continues to be the default
- AVONET regression check: < 5% RMSE change with `lambda_mode = "cv"`
- R CMD check clean
