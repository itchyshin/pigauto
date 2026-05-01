# Phase F: LP corner for ordinal traits — design

Date: 2026-05-01.
Status: design / proposed (not yet implemented).
Closes: AVONET Migration regression (−10.9 ± 3.1 pp across 3 seeds at
N_IMP=20, see `useful/MEMO_2026-05-01_multiseed_n20_and_default_flip.md`).

## Problem

For K=3 ordinal traits, the threshold-joint baseline (Phase B3, commit
`a541dbd`) produces systematically worse predictions than label
propagation (LP). The existing per-trait ordinal path selection (Opus
#6, commit `f43edb6`) chooses between two paths — `threshold_joint` and
`bm_mvn` — based on val MSE, but does **not** consider LP. As a result
the calibrated gate cannot route around bad ordinal-baseline behaviour
to LP.

Concrete symptom: AVONET Migration (3-level ordinal: 1 = sedentary,
2 = partial, 3 = full migrant) accuracy regresses by ~11 pp vs the
per-class-mode baseline at N_IMP=20, consistent across 3 seeds.

## Goal

Add LP as a third option in the existing per-trait ordinal path
selection. After this change, an ordinal trait's baseline becomes
the lower-val-MSE of {`threshold_joint`, `bm_mvn`, `lp`} per trait.

## Approach options considered

### Option A (recommended) — extend existing per-trait path selection

Mirror the existing Opus #6 pattern in `R/fit_baseline.R` lines
283–349:

- For each populated ordinal trait, compute LP-via-OVR on K classes
  (one binary LP per class), aggregate to expected-value-on-integer-
  class via `E[class] = sum_k p_k * k`, then z-score to match the
  scale of `threshold_joint` and `bm_mvn` outputs.
- Compute val MSE on the z-scored expected value.
- Pick the lowest-MSE of the three paths; record in
  `ordinal_path_chosen`.

**Pros**: minimal code change (~80 LOC); reuses existing LP
infrastructure (`R/baseline_lp.R`); simplex stays 3-D so no
calibrate_gates changes; orthogonal to safety-floor (still works
in the `safety_floor = FALSE` branch).

**Cons**: hard switch per trait (no blending of {BM, LP}). For traits
where threshold_joint and LP are both partly informative, no convex
combination is reachable. Empirically unlikely for K=3 — the bisect
showed threshold_joint is *worse than column mean*, so blending it
in would always hurt — but the option is not on the table.

### Option B — 4-D safety-floor simplex for ordinal

Make `simplex_grid()` 4-dimensional for ordinal traits: corners
(BM, GNN, MEAN, LP). Calibrator searches the 4-simplex.

**Pros**: full flexibility; calibrator can convex-combine any two of
{BM, LP} at varying weights.

**Cons**: ~7000-point grid at step=0.05 vs current 231; 30× more
calibration compute; touches `simplex_grid()`, `compute_corner_loss()`,
and the gate-application code in `predict_pigauto.R`. Bigger change.

### Option C — hybrid (per-trait baseline source + 3-D simplex unchanged)

Pick the better-of-{BM, LP} as the "baseline source" for ordinal
(Option A's selection), then keep the 3-D simplex with that selected
baseline as the BM corner. This is what Option A already does — it's
the same proposal under a different framing.

## Recommendation

**Option A.** The bisect evidence (`MEMO_2026-04-29_phase6_migration_bisect.md`)
shows threshold_joint is strictly worse than LP for K=3 Migration.
There's no signal that a BM × LP convex combination would help.
Option B's 30× calibration-compute hit is unjustified by current
evidence. If a future trait turns up where BM and LP are
complementary, we can upgrade to Option B then.

## Implementation sketch (Option A)

### Files touched

| file | lines | change |
|---|---|---|
| `R/fit_baseline.R` | ~30 LOC added inside the ordinal-selection block (lines ~313–349) | extend the `ordinal_path_chosen` selection to consider LP |
| `R/baseline_lp.R` | 0 LOC — reuse existing LP machinery | — |
| `tests/testthat/test-ordinal-baseline-paths.R` | NEW, ~80 LOC | regression test that an `ordinal_path_chosen == "lp"` outcome reduces val MSE on a fixture where threshold_joint is bad |
| `script/bench_avonet_phase_f.R` | NEW, ~120 LOC | post-merge AVONET Migration verification |
| `NEWS.md` | +20 lines | new `0.9.1.9010 (dev)` section |
| `DESCRIPTION` | 1 line | bump to `0.9.1.9010` |

### Key code shape

```r
# Inside the ordinal-path-selection block in fit_baseline.R:

# Existing: tj_mse and bm_mse computed.
# New: compute LP MSE.
lp_pred_class_expected <- ...   # E[class | x] from LP, per row
lp_pred_z <- (lp_pred_class_expected - tm$mean) / tm$sd
lp_diff <- lp_pred_z[val_rows_j[finite_t]] - truth_j[finite_t]
lp_mse  <- if (any(is.finite(lp_diff))) mean(lp_diff[is.finite(lp_diff)]^2) else NA_real_

# Pick lowest of three (existing logic generalised):
mses <- c(tj = tj_mse, bm = bm_mse, lp = lp_mse)
mses <- mses[is.finite(mses)]
chosen <- names(which.min(mses))
mu[, col] <- switch(chosen,
                    "tj" = tj_pred,
                    "bm" = bm_res$mu,
                    "lp" = lp_pred_z)
se[, col] <- switch(chosen,
                    "tj" = se[, col],          # unchanged
                    "bm" = bm_res$se,
                    "lp" = lp_se_estimate)     # SE estimate from LP probabilities
ordinal_path_chosen[as.character(col)] <- switch(chosen,
                                                 "tj" = "threshold_joint",
                                                 "bm" = "bm_mvn",
                                                 "lp" = "lp")
```

### LP for ordinal — the math

Run K independent OVR label propagations:
- For each class `k`, build a binary indicator `y_k = (class == k)` for
  observed rows (NA elsewhere).
- LP gives `p_k(s)` for every species `s`.
- Normalise across `k` so probabilities sum to 1: `p_k / sum_k p_k`.
- Expected class on the 1..K integer scale: `E[class] = sum_k p_k * k`.
- Z-score using `tm$mean`, `tm$sd` (the same scaling used by
  `threshold_joint` and `bm_mvn`).

For the SE estimate, use the variance of the K-class distribution:
`Var[class] = sum_k p_k * (k - E[class])^2`, then `sqrt()` and divide
by `tm$sd`.

### Tests

1. `test_that("ordinal path selection chooses LP when threshold_joint is misspecified", { ... })` — fixture with K=3 ordinal trait simulated such that LP wins on val MSE; verify `ordinal_path_chosen == "lp"`.

2. `test_that("ordinal path selection still chooses bm_mvn when LP is bad", { ... })` — fixture where LP is wildly off, BM wins.

3. `test_that("ordinal path selection produces finite predictions for all three options", { ... })` — sanity test that each path's output is finite + on the same z-scale.

4. Existing `test-ordinal*` tests must still pass.

### Bench

`script/bench_avonet_phase_f.R`:
- AVONET full n=1500, miss_frac=0.30, N_IMP=20, seeds 2030/2031/2032.
- Compare Migration accuracy pre-Phase-F vs post-Phase-F.
- Pre-registered success criterion: post-Phase-F Migration accuracy
  ≥ baseline (0.80) on all 3 seeds. Currently 0.71 ± 0.03 — needs to
  recover ~9 pp.

### Wall-time impact

Ordinal path selection runs once per ordinal trait at fit time. LP
adds K independent LP calls per ordinal trait. AVONET has 1 ordinal
(Migration, K=3) so ~3 LP calls; ~10 s extra at n=1500. Negligible.

## Risk

- **Risk that LP is also bad**: possible. The path selection picks the
  best of three on val MSE; if all three are bad, we'd still choose
  the least-bad. Worst-case Phase F is no-worse-than-current.
- **Generalisation beyond Migration**: K=4 ordinal traits in other
  benchmarks may show different relative orderings. The selection is
  per-trait so it should adapt automatically; bench coverage should
  include K=4 fixtures.

## Open

- Should LP be considered for ordinal at K > 3? The bisect localised the
  problem to small K. For K=10 ordinal the LP K-class-aggregation may
  be poorly calibrated. **Proposal**: enable for all K; let val MSE
  decide. If empirical evidence at K > 5 shows LP regresses, we can
  cap.
- Should we also add LP corner for binary? Binary already uses LP as
  the default baseline (`fit_baseline()` line ~205 routes binary to LP
  via `lp_baseline_binary()`). No change needed.
