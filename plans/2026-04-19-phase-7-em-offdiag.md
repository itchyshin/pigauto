# Phase 7 EM (off-diagonal) Implementation Plan

> Execute inline, commit per task. Stacks on top of `feature/phase-6-em`.

**Goal:** Add `em_offdiag = FALSE` opt-in. When `em_offdiag = TRUE` AND `em_iterations >= 2L`, each liability cell's prior at iteration `k+1` is the conditional-MVN `(mu, sd)` given the posterior liability of other traits at iteration `k`. Default behaviour = Phase 6.

**Files touched:**

| Path | Change |
|---|---|
| `R/joint_threshold_baseline.R` | `build_liability_matrix` gains `mu_prior_mat` + `sd_prior_mat`; new `build_conditional_prior()` helper; `fit_joint_threshold_baseline_em` gains `em_offdiag` arg |
| `R/fit_baseline.R` | Pass `em_offdiag` through to the EM wrapper |
| `R/impute.R` | Accept + pass-through `em_offdiag` |
| `tests/testthat/test-phase7-em.R` | Regression guard + correlated-binary convergence + missing-pattern fallback |
| `script/bench_phase6_em.R` | New smoke bench comparing Phase 6 vs Phase 7 |
| `NEWS.md` | Append Phase 7 to the v0.9.1.9000 section |

## Task 1: `build_liability_matrix` accepts per-cell mu/sd matrices

- Add `mu_prior_mat = NULL` and `sd_prior_mat = NULL` args (n × K_liab matrices).
- When non-NULL, use `mu_prior_mat[i, j]` and `sd_prior_mat[i, j]` in each `estep_liability*()` call (instead of the legacy scalar or vec-indexed values).
- Precedence: `mu_prior_mat`/`sd_prior_mat` > `sd_prior_vec` > default `sd=1`.
- Regression: `mu_prior_mat = NULL` + `sd_prior_mat = NULL` preserves Phase 6 behaviour byte-for-byte.

## Task 2: `build_conditional_prior(Sigma, L_hat)` helper

New function producing `(n × K) mu_prior_mat` and `(n × K) sd_prior_mat` from a K × K Σ and the n × K posterior-mean liability matrix from the previous iter. Handles:

- All-observed cells: full conditional-MVN formula.
- Any-NA cells: restrict the condition to the observed subset of other traits; recompute a smaller `solve()`.
- All-NA cells: fall back to unconditional prior `(mu = 0, sd = sqrt(diag(Σ)))`.

## Task 3: `fit_joint_threshold_baseline_em` gains `em_offdiag`

- New arg `em_offdiag = FALSE`.
- When `em_offdiag = TRUE`:
  - Iter 1 behaves exactly like Phase 6 (no previous Σ or L yet).
  - Iter 2+: call `build_conditional_prior(Sigma_k, L_k)` and pass `mu_prior_mat` + `sd_prior_mat` to `build_liability_matrix`.
- `em_state` gains `em_offdiag` field.
- Convergence criterion switches from `rel_frobenius(diag(Σ))` to `rel_frobenius(Σ)` (vectorised over all entries) when `em_offdiag = TRUE`.

## Task 4: Plumb through `fit_baseline()` and `impute()`

- Add `em_offdiag = FALSE` to both signatures. Pass through.
- Roxygen `@param em_offdiag`.

## Task 5: Tests in `tests/testthat/test-phase7-em.R`

1. Regression: `em_offdiag = FALSE` byte-identical to Phase 6 at all `em_iterations`.
2. `em_offdiag = TRUE, em_iterations = 1L` falls through to Phase 6 (no previous Σ).
3. `em_offdiag = TRUE, em_iterations = 5L` on synthetic correlated-binary (ρ = 0.6, n = 200): runs without error, `em_state$em_offdiag == TRUE`, `em_state$iterations_run >= 2L`.
4. `build_conditional_prior`: sanity on small inputs + missing-pattern fallback.

## Task 6: Smoke bench in `script/bench_phase6_em.R`

- AVONET 300 mixed-type + synthetic correlated-binary sim.
- Three rows: plug-in (em=0), Phase 6 diagonal (em=5, offdiag=F), Phase 7 offdiag (em=5, offdiag=T).
- Per-trait accuracy, RMSE, wall time.

## Task 7: NEWS + check + push + PR

- Append Phase 7 paragraph to the v0.9.1.9000 NEWS section.
- Full test suite green.
- `devtools::check()` 0 errors / 0 warnings / ≤1 note.
- Push + `gh pr create`.
