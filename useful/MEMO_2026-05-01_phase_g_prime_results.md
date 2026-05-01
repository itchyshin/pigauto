# Phase G' (PMM) results — honest reading: tradeoff, not strict win

Date: 2026-05-01.
Implements `match_observed = "pmm"` (Predictive Mean Matching) as
the principled alternative to Phase G's `clamp_outliers` heuristic.
**Bench shows PMM is a tradeoff, not a strict improvement.**  Ship
opt-in; do NOT flip default.

## Pre-registered criteria + bench result

| criterion | result |
|---|---|
| 1. PMM Mass seed-2030 RMSE ≤ clamp Mass seed-2030 RMSE | **PASS dramatically** (1,275 vs 6,277; PMM 5x better than clamp) |
| 2. PMM doesn't regress > 5 % vs none on any (dataset, seed, trait) | **FAIL on 5 cells** (worst: AVONET Mass seed-2032, +805 %) |
| 3. PMM imputed values 100 % in observed range | **PASS** (programmatic guarantee held on every cell) |

**Verdict**: PHASE G' PARTIAL.  Ship PMM as opt-in.  Default unchanged.

## The mixed result on AVONET Mass

| seed | none RMSE | clamp RMSE | pmm RMSE | PMM Δ vs none |
|---|---|---|---|---|
| **2030** (Casuarius blow-up) | **29,930** |  6,277 | **1,275** | **−95.7 %** |
| 2031                         |    309    |    303 |   **259** |  −16.2 %    |
| **2032** (already accurate)  |    **374**|    660 | **3,386** | **+805 %**  |

**Mechanism behind the seed-2032 regression**: pigauto's BM + GNN +
blend pipeline produces RMSE 374 on seed-2032 Mass — about 14×
below the column-mean baseline (5,269).  These predictions are
already excellent.  PMM forces each missing cell's predicted value
to snap to the K=5 nearest observed donor.  When predictions are
already on-target, donor-mismatch noise (the difference between a
prediction and its closest observed-truth value) inflates RMSE.

In other words, PMM trades **typical-case accuracy** for
**tail-extrapolation safety**.  It is a downside-protection tool,
not a free lunch.

## Other AVONET log-cont traits (Beak, Tarsus, Wing)

PMM was within 5 % of `none` for most cells across all 3 seeds, with
small regressions of 5–15 % on a few cells.  No PMM gains either —
these traits don't have the extreme-tail mode that Mass exhibits.
The takeaway: **PMM helps only where predictions are unreliable**.

## Synthetic heavy-tail bench

3 reps, n = 200, 5 % outlier species at 10–50× normal mass.
PMM was ~5 % worse than `none` on all 3 reps.  The outliers in this
synthetic weren't extreme enough to trigger expm1 amplification (the
"none" predictions stayed in observed range — `clamp` was identical
to `none` on this dataset, meaning the clamp threshold was never
breached).  PMM imposed its 5 % cost without delivering the
tail-safety upside.

## When to use PMM (revised recommendation)

The original framing — "PMM should be the default for at-risk types"
— was wrong.  Honest framing:

| context | tool |
|---|---|
| Routine continuous imputation; predictions look reasonable | **default `none`** — accept small expm1 risk |
| Suspect extrapolation (rare clade, isolated tip, log-cont with extreme observed range) | **`match_observed = "pmm"`** — 5–15 % accuracy cost in typical cells for up to 95 % safety on outliers |
| Multi-imputation context (Rubin's rules via `pool_mi()`) | **`match_observed = "pmm"`** — between-draw donor variance is properly calibrated, matching mice convention |
| Want a deterministic point estimate with cheap tail safety | **`clamp_outliers = TRUE`** — captures big blow-ups without donor stochasticity |

## Why the original "PMM should be default" intuition was wrong

The mice convention is to use PMM for continuous imputation in
**multi-imputation contexts** specifically, where between-draw
variance preservation matters more than point-estimate accuracy.
mice users typically run M = 5–20 imputations and use Rubin's rules
to propagate uncertainty.  PMM's donor-mismatch noise is a feature,
not a bug, in that workflow: it gives Rubin's rules properly
calibrated standard errors.

Pigauto's typical workflow is more variable.  Some users want
single-imputation point estimates (where PMM's noise hurts).
Others want MI for downstream comparative analysis (where PMM's
noise is properly calibrated).  Default choice can't satisfy both.

**Conclusion**: the bench correctly flagged the tradeoff.  Default
stays at `match_observed = "none"`; users opt in based on workflow.

## Implementation summary

Ship Phase G' as a **complete and well-tested optional feature**:

* `R/pmm.R` — 4 helpers (`pmm_impute_one_trait`,
  `apply_pmm_to_decoded`, `recover_X_orig`, `pmm_is_eligible`)
* `R/predict_pigauto.R`, `R/impute.R` — `match_observed` and `pmm_K`
  args wired through.
* `R/fit_pigauto.R` — fit object retains `X_scaled` so PMM can
  recover original-units observed values for the donor pool.
* `tests/testthat/test-pmm.R` — **65 assertions across 5 layers**,
  all pass.
* `script/bench_phase_g_prime_pmm.R` — AVONET 4 log-cont × 3 seeds +
  synthetic heavy-tail × 3 reps × {none, clamp, pmm}.
* DESCRIPTION + NEWS bumped to `0.9.1.9012`.

## Files

* `R/pmm.R` (NEW)
* `R/predict_pigauto.R`, `R/impute.R`, `R/fit_pigauto.R` (modified)
* `man/predict.pigauto_fit.Rd`, `man/impute.Rd` (regenerated)
* `tests/testthat/test-pmm.R` (NEW)
* `script/bench_phase_g_prime_pmm.{R,md,rds}`
* `useful/MEMO_2026-05-01_phase_g_prime_results.md` (this memo)
* `DESCRIPTION` + `NEWS.md`

## What's NOT in this PR

* Default flip — bench evidence does not justify (seed-2032 Mass
  shows real downside risk).
* Cross-dataset evidence beyond AVONET + synthetic — could be added
  later if the PMM tradeoff needs more characterisation.
* Phase H+ (cross-K ordinal mode bench) — independent; queued for
  next.

## Open questions for future work

1. **Why does seed-2032 Mass behave so differently from seed-2030?**
   Both have N_IMP=20 and the same fitting pipeline.  Seed-2032 has
   no Casuarius-class outlier in the test set; the 30 % MCAR mask
   produced a "well-behaved" set of held-out cells.  Worth a quick
   diag: dump the per-cell residuals to see if PMM regresses
   uniformly or only on a few cells.

2. **Is K = 5 the right donor pool size?**  Mice's default of 5 is
   convention.  Larger K dilutes donor-mismatch noise (averages over
   more donors); smaller K (e.g. K = 1) is more deterministic but
   variance-reduced.  Could sweep K ∈ {1, 3, 5, 10, 20} on AVONET
   to find the AVONET-optimal value.  Probably future Phase G''.

3. **Should PMM use covariates in the donor-distance metric?**
   Currently uses only the prediction value as a 1-D distance.
   Could include trait correlations (mass-correlated traits like
   Beak.Length) for finer matching.  Big change; not for this PR.
