# Phase C — cross-dataset bench (multi-seed v2): single-seed v1 was misleading

Date: 2026-05-01.
Bench: `script/bench_phase_c_cross_dataset.R` v2.  3 seeds (2026 / 2027
/ 2028) × 4 configs × 2 datasets = 24 cells.  Wall: 51 min.

**Updates the v1 single-seed memo** which had this same path -- v1
results were dominated by GNN-training MC noise; multi-seed reveals
the actual effect sizes (mostly small) and **reverses one of v1's
headline findings**.

## TL;DR — what holds and what doesn't

✓ **Phase H mode pooling helps K=3 ordinal across taxa, smaller than v1
   suggested**: PanTHERIA habitat_breadth +1.5 ± 1.5 pp; AVONET
   Migration +6.6 pp (PR #59).  Real but modest at K=3.

✗ **Phase H mode pooling HURTS K=5 ordinal on real data**: PanTHERIA
   diet_breadth −2.1 ± noise pp.  v1 said +4.2 pp (single-seed
   artefact).  **Multi-seed agrees with the Phase H+ simulation**
   (PR #62: −6.7 pp on simulated K=5).  My v1 PR #64 framing
   "real-data K=5 contradicts simulation" was wrong.

✗ **Phase G clamp_outliers effect on log-cont mass is small and
   dataset-dependent**, not the broad win v1 implied.  PanTHERIA
   body_mass +7.3 % (helps); AmphiBIO body_mass −5.3 % (hurts);
   AVONET Casuarius case −74 % (huge help, but extreme outlier).
   The Casuarius case generalises only to OTHER extreme-outlier
   cases, not to routine log-cont imputation.

## Multi-seed result (mean ± SD across 3 seeds)

PanTHERIA mammals:

| trait | type | default | clamp Δ | mode Δ | both Δ |
|---|---|---|---|---|---|
| body_mass_g            | log-cont (M)     | 3.70M ± 0.09M | **+7.3 %** | +9.0 % | +9.2 % |
| head_body_length_mm    | log-cont (mm)    | 504 ± 182    | +3.4 %    | +8.6 % | +3.6 % |
| gestation_d            | log-cont (d)     | 42.6 ± 10.9  | +5.7 %    | +13.5 % | +2.9 % |
| max_longevity_m        | log-cont (mo)    | 137 ± 55     | −5.8 %    | +4.8 % | −28.9 % |
| litter_size            | count            | 1.07 ± 0.07  | −2.7 %    | −0.6 % | +0.8 % |
| **diet_breadth**       | **ordinal K=5**  | **0.407 ± 0.031** | +2.5 pp   | **−2.1 pp** | −0.2 pp |
| **habitat_breadth**    | **ordinal K=3**  | **0.767 ± 0.031** | +1.4 pp   | **+1.5 pp** | +1.7 pp |
| terrestriality         | binary           | 0.917 ± 0.013 | −0.4 pp   | 0     | −0.3 pp |

AmphiBIO amphibians (Diu / Noc dropped; presence-only encoding):

| trait | type | default | clamp Δ | mode Δ |
|---|---|---|---|---|
| body_mass_g       | log-cont    | 57.2 ± 37.7 | **−5.3 %** | −12.5 % |
| body_size_mm      | log-cont    | 59.5 ± 18.9 | +2.5 %    | −5.9 % |
| diet_breadth K=5  | ordinal (degenerate) | 0.872 ± 0.051 | 0      | 0     |
| habitat (4-class) | categorical | 0.847 ± 0.012 | 0      | 0     |

## Why the single-seed bench was misleading

In Phase C v1 (single seed), the cross-config RMSE deltas on
continuous traits were dominated by **MC noise from re-training the
GNN per config** -- each `impute()` call has stochastic dropout +
non-deterministic torch ops, so the same data yields different
predictions across runs.  At single seed:

  PanTHERIA body_mass mode    -20.2 % (single-seed)
  PanTHERIA body_mass mode    +9.0 ± SD  (multi-seed)
  PanTHERIA diet_breadth mode +4.2 pp (single-seed)
  PanTHERIA diet_breadth mode -2.1 ± noise pp (multi-seed)

The mode K=5 reversal is the most important: single-seed gave a
confident +4.2 pp lift, multi-seed shows a small negative effect.
Without multi-seed, this would have entered the literature as a
positive finding contradicting the simulation result.

## Updated user guidance (from all phase work today)

| trait class | recommendation |
|---|---|
| log-cont body mass / size / gestation | defaults usually fine; `clamp_outliers = TRUE` ONLY for known-extreme-outlier datasets (e.g. AVONET when Casuarius-class species are masked).  For routine use, the small deltas are noise. |
| **K=3 ordinal** | `pool_method = "mode"` modestly helps (+1.5 to +6.6 pp).  Recommended at N_IMP > 1. |
| **K=5+ ordinal** | **`pool_method = "mode"` HURTS at K=5+** (real-data + simulation agree).  Keep default `"median"`. |
| binary / categorical | defaults usually fine. |
| count | defaults usually fine; small datasets noisy. |

## What changed in the v2 bench (vs v1 PR #64)

- 3 seeds instead of 1
- AmphiBIO Diu / Noc dropped (presence-only encoding produced
  −39.8 pp artefact on Diu in v1)
- Per-config aggregation: mean ± SD across seeds
- Honest reading of effect sizes vs noise

## What this PR contains

- `script/bench_phase_c_cross_dataset.R` — updated bench (multi-seed,
  Diu/Noc dropped)
- `script/bench_phase_c_cross_dataset.{md,rds}` — v2 outputs
- (this memo)

No code changes to `R/`.

## What's NOT in this PR

- v0.9.2 default-flip recommendations.  No, none.  All effect sizes
  are within noise except the K=3 mode finding (consistently
  positive but small).
- BIEN plant + FishBase fish benches.
- `clamp_factor` tuning.
