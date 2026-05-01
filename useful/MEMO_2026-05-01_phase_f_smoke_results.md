# Phase F smoke bench results — AVONET Migration n=1500

Date: 2026-05-01.
Bench: `script/bench_avonet_phase_f.R` (committed today; ~10 min wall).
Prior commit: `71b4feb` (Phase F merged via PR #54).

**TL;DR**: Phase F is mechanically correct but **does not close the
AVONET Migration regression**. LP is never picked on this fixture
(BM-via-MVN wins val MSE on all 3 seeds). Phase F adds optionality
without changing observed behaviour. Most of the previously-reported
−10.9 pp regression at N_IMP=20 is plausibly an **MC-dropout
discrete-trait pooling artifact**, not a baseline misspecification.

## Results

| seed | pigauto acc | baseline acc | lift (pp) | path chosen |
|---|---|---|---|---|
| 2030 | 0.756 | 0.789 | −3.3 | bm_mvn |
| 2031 | 0.776 | 0.798 | −2.2 | bm_mvn |
| 2032 | 0.771 | 0.813 | −4.2 | bm_mvn |

**Mean**: pigauto 0.767, baseline 0.800, lift = **−3.3 ± 1.0 pp**.

Pre-registered Phase F target was Migration acc ≥ 0.78 on ≥ 2 of 3
seeds. **0 / 3** seeds reached this target. Phase F **inconclusive**.

## What Phase F achieved

* The selection logic now considers three options
  ({threshold_joint, bm_mvn, lp}) per ordinal trait, picks
  lowest val MSE. ✓
* On AVONET-300 (the existing regression bottle test):
  Migration still picks `bm_mvn`. ✓ No regression introduced.
* On AVONET-1500 (this bench): Migration still picks `bm_mvn`
  on all 3 seeds. LP loses on val MSE.
* Architecture is additive — if a future trait shows LP
  beats `bm_mvn`, the selection will route to it automatically.

## What Phase F did NOT achieve

* **LP did not displace `bm_mvn`** for AVONET Migration on either
  n=300 or n=1500. The −3.3 pp regression at N_IMP=1 remains.
* The regression is **smaller** than the multi-seed memo reported
  for N_IMP=20 (−10.9 ± 3.1 pp). This means most of the
  multi-seed regression is **N_IMP-dependent** — it gets worse
  with more imputations, not better.

## New finding: N_IMP-dependent regression

| Setting       | pigauto Migration acc |
|---|---|
| N_IMP=1 (this bench)   | 0.767 ± 0.011 |
| N_IMP=20 (multi-seed memo) | 0.713 ± 0.032 |
| Mean baseline (per-class mode) | 0.800 ± 0.013 |

The pigauto acc *drops by ~5 pp* when going from N_IMP=1 to N_IMP=20.
Pure baseline misspecification (which Phase F targeted) cannot
explain this — the baseline doesn't change with N_IMP. The
implication is that **MC-dropout pooling for discrete traits is
introducing predictive noise that hurts rather than helps**.

Hypothesis: at higher N_IMP, the GNN is trained with dropout active
on each forward pass. For a 3-class ordinal, dropout-induced
prediction noise can flip enough cells to the wrong adjacent class
to materially hurt accuracy after pooling. The current pool method
for discrete traits is `mean` of probability vectors — adequate for
calibrated probabilities, suboptimal for hard accuracy.

## Proposed Phase H — discrete-trait pooling fix (separate sprint)

**Investigation queue** (~1 sprint):

1. Inspect `R/predict_pigauto.R` discrete-trait pooling code path —
   confirm the current behaviour for ordinal traits at
   N_IMP > 1.
2. Bench three pooling alternatives:
   * `mean` (current default for discrete probabilities)
   * `median` of per-cell probability vectors
   * `mode` — the most frequent class label across draws
3. Targeted Migration accuracy test at N_IMP ∈ {1, 5, 10, 20}
   on AVONET-1500 with each pooling method.
4. Acceptance criterion: Migration acc at N_IMP=20 with the new
   pooling method ≥ Migration acc at N_IMP=1 (i.e., more imputations
   should help or be neutral, never hurt).

If Phase H closes the gap to baseline, ship it as v0.9.2.
Otherwise the regression at N_IMP=1 (−3.3 pp) is a separate
real-but-small misspecification on Migration that needs its own
investigation (Phase I — possibly a Migration-specific feature
engineering or covariate-aware baseline).

## Should Phase F be reverted?

No. Phase F is purely additive (LP becomes an option, val MSE picks
best of three). It does not change observed behaviour on AVONET
because LP loses there, but it gives the package optionality for
future traits and other datasets where LP might win. The cost is
~30–60 s extra per ordinal trait at fit time (K LPs) which is
negligible relative to typical pigauto wall time.

The one genuine harm from Phase F is **a misleading commit message**
suggesting it closes the Migration regression. PR #54 was conservative
enough to ship the architecture without claiming the closure (no
NEWS/DESCRIPTION bump until bench confirms). This memo replaces that
implicit claim with the honest verdict.

## Files

* `script/bench_avonet_phase_f.R` — the smoke bench
* `script/bench_avonet_phase_f.{rds,md}` — outputs
* upstream: `useful/MEMO_2026-05-01_multiseed_n20_and_default_flip.md` — original multi-seed evidence
* upstream: `specs/2026-05-01-phase-f-lp-corner-ordinal-design.md` — Phase F design memo (now partially superseded by these results)
