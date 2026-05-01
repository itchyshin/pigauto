# Phase H+ — cross-K ordinal mode pooling: INCONCLUSIVE

Date: 2026-05-01.
Bench: `script/bench_phase_h_plus_cross.R` on simulated K-class
ordinal traits at λ = 0.6 phylogenetic signal.  Wall: ~14 min.

**TL;DR**: Mode pooling beats median-mean+round on 7 / 9 cells (78%,
below the pre-registered 80% threshold) and exhibits a worst-case
−23.3 pp regression at K=5 rep=2.  **Pre-registered criteria FAIL.**
Default for ordinal pool stays `median`.

This is a **negative result that's the right call**: Phase H (PR #59)
correctly closed the AVONET Migration K=3 regression by adding
`pool_method = "mode"` as opt-in, but cross-K evidence does not
support flipping the default.  Mode is K-3-friendly but unreliable
across the K-spectrum.

## Pre-registered acceptance criteria

> 1. Mode acc ≥ median acc on at least 80 % of cells at N_IMP = 20.
> 2. Mode never regresses by > 2 pp on any (K, rep) combo.

## Result

Per-cell:

| K | rep | median | mode | Δ (pp) |
|---|---|---|---|---|
| 3 | 1 | 0.667 | 0.717 | +5.0 ✓ |
| 3 | 2 | 0.367 | 0.367 | 0    |
| 3 | 3 | 0.517 | 0.517 | 0    |
| 5 | 1 | 0.400 | 0.433 | +3.3 ✓ |
| **5** | **2** | **0.550** | **0.317** | **−23.3** ✗ |
| 5 | 3 | 0.467 | 0.467 | 0    |
| 7 | 1 | 0.350 | 0.350 | 0    |
| 7 | 2 | 0.250 | 0.233 | −1.7 |
| 7 | 3 | 0.183 | 0.300 | +11.7 ✓ |

Aggregated by K:

| K | median | mode | mode Δ |
|---|---|---|---|
| **3** | **0.517** | **0.533** | **+1.6 pp** ✓ |
| **5** | **0.472** | **0.406** | **−6.7 pp** ✗ |
| 7 | 0.261 | 0.294 | +3.3 pp |

**Verdict**: Phase H+ INCONCLUSIVE.  Mode wins or ties on 7/9 cells
(target ≥ 80%, achieved 78%); worst per-cell regression −23.3 pp
violates the −2 pp guard.  Default `pool_method = "median"` for
ordinal traits **stays unchanged**.

## Pattern: Phase H+ is qualitatively similar to Phase G' (PMM)

Both PRs found the same thing: **the architectural fix is real but
context-dependent**.

| | Phase H (PR #59) | Phase H+ (this) | Phase G' (PR #61) |
|---|---|---|---|
| When mode/PMM helps | AVONET Migration K=3 N_IMP=20 (+6.6 pp) | Simulated K=3 (+1.6 pp), K=7 (+3.3 pp) | AVONET Mass seed-2030 (−95.7% RMSE) |
| When it hurts | --- | **K=5 (−6.7 pp)** | **AVONET Mass seed-2032 (+805% RMSE)** |
| Recommendation | ship as opt-in, no default flip | no default flip | ship as opt-in (Phase G'' may scope down) |

## Why mode regresses at K=5 (mechanism guess)

K=5 is in the awkward middle: enough classes for mode pooling to
"spread thin" (each class gets ~M/K = 4 votes per cell, low signal
per class), but not so many classes that the integer-mean+round
default has many adjacent-class flip opportunities.

K=3 (Migration regime): mode benefits because mean+round biases
toward the middle class.  K=7+: mode benefits because mean-and-round
spans many classes; mode preserves the modal-class signal.  K=5:
neither benefit applies; mode just adds variance from M/K-thin per-
class voting.

This is an empirical observation, not a proof.  Could be a fixture
artifact (seed=2032 effect) or genuinely K=5-specific.  Not chasing
this rabbit hole; the substantive finding is that cross-K mode is
NOT universally better.

## What this means for users

Per `?fit_pigauto` (post Phase H, PR #59):

* **K=3 ordinal traits with `n_imputations > 1`**: pass
  `pool_method = "mode"` explicitly — this closes the −10.9 pp
  Migration regression on AVONET (Phase H bench evidence).
* **K=5 ordinal**: stick with default `median` (mean+round).
  Phase H+ shows mode loses on K=5.
* **K=7+ ordinal**: try `pool_method = "mode"` — Phase H+ shows it
  wins or ties.
* **Single-imputation (n_imputations = 1)**: pool_method has no
  effect.

Future work: a per-trait auto-selection (like the per-trait
ordinal-baseline-path-selection in `fit_baseline.R`) could pick
mode vs median based on val MSE per ordinal trait.  Out of scope
for this PR.

## What this PR contains

* `script/bench_phase_h_plus_cross.R` (already on main from earlier;
  no diff)
* `script/bench_phase_h_plus_cross.md` (NEW — bench output)
* `useful/MEMO_2026-05-01_phase_h_plus_results.md` (this memo)

No R/ code changes.  No NEWS entry.  Pure evidence-only PR.

## Files

* `script/bench_phase_h_plus_cross.R` — bench (already on main)
* `script/bench_phase_h_plus_cross.{md,rds}` — bench output
* (this memo)
* upstream: PR #55 — Phase F bench results (Phase H proposal)
* upstream: PR #59 — Phase H (mode pooling implementation)
