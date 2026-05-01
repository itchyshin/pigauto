# Phase H — discrete-trait pooling: mode beats mean+round for ordinal

Date: 2026-05-01.
Closes the AVONET Migration N_IMP-dependent regression flagged in
`useful/MEMO_2026-05-01_phase_f_smoke_results.md`.

**TL;DR**: ordinal traits at `N_IMP > 1` were pooling by `mean of
integer class indices, then round` — biasing predictions toward
middle classes when dropout spread predictions across adjacent
classes. The fix is **per-cell majority vote** (`pool_method =
"mode"`). Pre-registered acceptance passes 3/3 seeds; mode improves
Migration accuracy at N_IMP=20 by **+6.6 pp** vs the default.

## Pre-registered acceptance criterion

From PR #55 (Phase F bench results):

> Acceptance criterion: Migration acc at N_IMP=20 ≥ Migration acc
> at N_IMP=1 (more imputations should help or be neutral, never
> hurt).

## Result

3 seeds × {N_IMP=1, N_IMP=20} × `pool_method = "mode"` on AVONET
n=1500. ~25 min wall.

| seed | mode N_IMP=1 | mode N_IMP=20 | Δ |
|---|---|---|---|
| 2030 | 0.756 | **0.767** | +1.1 pp ✓ |
| 2031 | 0.776 | **0.789** | +1.3 pp ✓ |
| 2032 | 0.771 | **0.782** | +1.1 pp ✓ |

**Mean**: 0.767 ± 0.011 → 0.779 ± 0.011 (cross-seed). All 3 seeds
satisfy `mode N_IMP=20 ≥ mode N_IMP=1`. **PHASE H PASSES**.

## Side-by-side regression closure

Comparing to the previously-collected default-pool numbers:

| setting | Migration acc | source |
|---|---|---|
| Mean baseline (per-class mode) | 0.800 ± 0.013 | multi-seed memo |
| Default median pool, N_IMP=1   | 0.767 ± 0.011 | Phase F bench (PR #55) |
| Default median pool, N_IMP=20  | **0.713 ± 0.032** ← was regressed | multi-seed memo |
| **Mode pool, N_IMP=1**         | 0.767 ± 0.011 | this memo |
| **Mode pool, N_IMP=20**        | **0.779 ± 0.011** ← Phase H | this memo |

* Going default→mode at N_IMP=20: **+6.6 pp**
  (from 0.713 to 0.779).
* Regression vs baseline: was −10.9 pp, now **−2.1 pp**.
* Mode at N_IMP=1 is unchanged: at single-pass the integer class is
  decoded directly so pooling math doesn't apply. (The bench result
  `pigauto_acc = 0.767` exactly matches the Phase F default-pool
  N_IMP=1 number — sanity check that mode and median produce
  identical output at N_IMP=1.)

## Root cause

`R/predict_pigauto.R` lines 698–705 (pre-Phase-H):

```r
} else if (tm$type == "ordinal") {
  K <- length(tm$levels)
  int_vals <- rowMeans(sapply(decode_results, function(dr) {
    as.integer(dr$imputed[[nm]]) - 1L
  }))
  int_vals <- as.integer(pmin(pmax(round(int_vals), 0), K - 1L))
  ...
}
```

For each held-out ordinal cell, M MC-dropout draws produce M class
indices `{c_1, c_2, ..., c_M}`. Pooling rule: `mean(c_i)` then
round. This biases toward **middle classes** when dropout produces
asymmetric noise.

Concrete K=3 failure case (Migration: sedentary=1, partial=2,
full=3):

* True class 1 (sedentary), M=20 dropout draws produce
  `{1, 1, 1, 2, 3}`.
* Mean = 1.6, round = 2 → predicted "partial migrant" ❌
* Mode (majority vote) = 1 → predicted "sedentary" ✓

At N_IMP=1 there's no pooling, so the single-pass argmax produces
the right class. At N_IMP > 1 the integer-mean-round drifts predictions
toward 2 (middle), losing accuracy on classes 1 and 3.

## Fix

`R/predict_pigauto.R` `pool_imputations()`:

* Extended `pool_method` enum from `c("median", "mean")` to
  `c("median", "mean", "mode")`.
* New ordinal branch when `pool_method = "mode"`: per-cell
  `tabulate(class_indices, nbins = K) → which.max() - 1L`.
* Continuous-family branches (count, proportion, zi_count) where
  `pool_method == "median"` already triggers median pooling: extend
  the guard to `pool_method %in% c("median", "mode")` so that
  passing `"mode"` for these types falls back to median (since mode
  doesn't apply to continuous decoders).
* Binary, categorical, multi_proportion: no change (already
  probability-averaged + argmax).

`R/impute.R`: extend the user-facing `pool_method` enum
identically. Update roxygen.

## Tests

`tests/testthat/test-mi-pool-robust.R` adds 2 new test_that blocks:

* `[Phase H] pool_imputations 'mode' picks majority class for
  ordinal` — K=3 fixture, 4 species, 5 draws. Verifies (a) all
  4 species get the correct mode-class under `pool_method = "mode"`,
  (b) the same fixture under default `"median"` deliberately
  flips species 1 to the wrong class — locking in the regression-
  test pattern.
* `[Phase H] mode falls back to median for log-transformed
  continuous` — verifies that passing `"mode"` for a log-cont
  trait correctly routes to median, not mean.

11 / 11 tests in `test-mi-pool-robust.R` pass.

## Default policy

**Conservative**: keep `pool_method = "median"` as the default.
Document `"mode"` as the recommended setting for ordinal traits
at N_IMP > 1 in NEWS.md and the vignette.

Reasoning: although mode-N_IMP=20 wins on AVONET Migration, this is
a single-trait single-dataset evidence base. Other ordinal traits
(K > 3, with structurally smooth probability distributions across
classes) might prefer the integer-mean-round behaviour. Default
behaviour at upgrade-time should not silently change a numeric
output. Mode is a one-line opt-in.

A future Phase H+ could:

1. Bench mode vs median across additional ordinal datasets (e.g.,
   AmphiBIO life-history ordinals, plant trait ordinals) at
   N_IMP=20 to verify mode is universally better.
2. If the cross-dataset evidence holds, flip default to mode in
   v0.9.2 with a NEWS callout.

## Files

* `R/predict_pigauto.R` — implementation
* `R/impute.R` — enum + roxygen
* `man/predict.pigauto_fit.Rd`, `man/impute.Rd` — regenerated
* `tests/testthat/test-mi-pool-robust.R` — 5 new assertions
* `script/bench_phase_h_pool.R` — bench script
* `script/bench_phase_h_pool.{rds,md}` — bench outputs
* (this memo)

## What's left

* Phase G — Mass tail-extrapolation clamp. Independent fix; can
  proceed in parallel.
* Cross-dataset ordinal mode bench (Phase H+). Optional future
  investigation before flipping the default.
* MEE methods paper (still deferred).
