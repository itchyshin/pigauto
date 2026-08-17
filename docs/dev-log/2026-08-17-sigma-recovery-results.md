# Known-Σ recovery simulation — results

2026-08-17 · Claude lane · design pre-registered at `2026-08-17-sigma-recovery-design.md`
(committed `75263ea`, before any estimator comparison ran). 12 cells × 50 reps × 6 arms,
**0 failures**, ~72 min local. Runner + log + rds:
`dev/sigma_recovery_sim.{R,log,rds}` on `arc/inhouse-sigma-convergence`.

Arms: **E0** current in-house (`single_pass`, `max_iter = 0`) · **E0i1** same + 1 refinement
iteration · **E1** Fisher-ML port · **E1i1**/**E1i3** Fisher-ML + 1/3 refinement iterations ·
**E2** `Rphylopars::phylopars` REML (yardstick).

## The finding: Σ is estimated, then discarded

**E0 and E1 produce numerically identical predictions in all 12 cells** — paired RMSE
difference exactly `0.0000` everywhere — despite Σ̂ estimates that differ substantially
(mean Frobenius error 0.315 vs 0.193 at λ=0.2). That is only possible if the Σ estimate
never enters the prediction. Confirmed in code: `fit_mvn_bm_inhouse()` short-circuits at
`if (K == 1L || max_iter <= 0L)`, returning **per-column BM predictions** and reporting Σ̂
as a by-product.

So the shipped joint baseline is not joint at prediction time. It fits a cross-trait
covariance and then predicts each trait independently. **This is the mechanism behind the
AVONET continuous gap** (`2026-08-16-continuous-gap-diagnosis.md`): not a bad Σ, but an
unused one.

It also settles the slice's original question: **improving the Σ estimator cannot help
while `max_iter = 0`.** The Fisher-ML port was aimed at the wrong target.

## Does using Σ close the gap? Mostly yes — and it breaks the intervals

Mean prediction RMSE (lower better), and mean SE coverage at nominal 0.95:

| arm | RMSE λ=0.2 | RMSE λ=1.0 | coverage |
|---|---:|---:|---:|
| E0 (shipped) | 0.988 | 0.498 | 0.925 |
| E0i1 | 0.920 | 0.497 | 0.856 |
| E1i1 | 0.857 | 0.480 | 0.831 |
| E1i3 | 0.835 | **0.521** | **0.618** |
| **E2 (phylopars)** | **0.812** | **0.413** | **0.935** |

- At λ=0.2, refinement closes ~85% of the E0→E2 gap (0.988 → 0.835 vs 0.812).
- At λ=1.0 it barely helps, and E1i3 makes things *worse* (0.521 > 0.498).
- **Coverage degrades monotonically with refinement**: 0.925 → 0.856 → 0.618. The refined
  posterior mean moves but `anc_var` is not updated to match, so the intervals become
  badly overconfident. This is a correctness problem, not a tuning knob — and it is
  almost certainly what the 2026-05-17 "EM diverged" note was seeing from a different angle.

## Fisher-ML verdict: does NOT ship

Per the pre-registered rule, with one honest amendment (below):

- **Rule 2 (prediction within 3×MCSE of E2): FAILS.** E1 is bit-identical to E0; both are
  far from E2.
- **Rule 3 (coverage ∈ [0.93, 0.97]): FAILS** for every refined arm; E0/E1 sit marginally
  low at 0.925.
- **Rule 1 (Frobenius ≤ 1.25 × E2): CANNOT BE EVALUATED — my metric was ill-posed.**
  E2's `phylocov` is on a different scale/parameterisation from the in-house `phylocov`;
  its Frobenius error sits at ~1.08–1.15 in every cell regardless of design, which is a
  units artifact, not accuracy. A tree-height normalisation did not fix it. **Ignore the
  E2 Frobenius column entirely.** E0-vs-E1 Frobenius *is* comparable (same code path, same
  convention) and is reported below on that basis only.

**Where Fisher-ML does and doesn't work** (E0 vs E1, comparable):

| | Frobenius λ=0.2 | Frobenius λ=1.0 | sign acc. λ=1.0 |
|---|---:|---:|---:|
| E0 | 0.315 | **0.353** | **0.993** |
| E1 | **0.193** | 0.477 | 0.912 |

Fisher-ML estimates Σ better at low phylogenetic signal and **worse at λ=1**, where its
off-diagonal sign accuracy drops to 0.91. Exactly the failure the design pre-registered as
its discriminating contrast: the estimator treats species as iid, which is least true when
signal is maximal. The v0.9.2 note that replaced an iid init for biasing off-diagonals was
right, and this measures the same defect from the other end.

## What this means for the dependency question

Shinichi's constraint (2026-08-17) is that pigauto should not depend on Rphylopars — the
in-house solver exists by design (PR #100). This sim makes the dependency-free path
**precise and bounded**, which it was not this morning:

1. **Use the Σ that is already computed** — enable cross-trait refinement in the prediction
   step (the machinery exists; `max_iter` is the switch).
2. **Fix the posterior variance under refinement** so coverage survives. This is the real
   work, and the sim now gives it a pass/fail target (coverage back to ~0.93–0.95 at
   E0i1-level RMSE).
3. **Keep `single_pass` Σ, not Fisher-ML** — E0's Σ is better where it matters (high signal)
   and Fisher-ML brings nothing once (1) is done.

That is a well-scoped numerical slice, not a rewrite, and it does not add a dependency.
`joint_solver = "rphylopars"` stays as the measured yardstick until it lands.

## Regime and caveats

- Simulated matrix-normal BM exactly (the model all arms assume), K=4, MCAR 30%,
  single-obs, n ∈ {100, 300}, λ ∈ {1.0, 0.2}, 3 Σ designs, 50 reps. No GNN, no
  mixed types — this is a solver-level study.
- Coverage here is the **BM analytic SE**, not conformal intervals (which are calibrated
  downstream and would mask this defect).
- The λ=1.0 arms use the untransformed tree; λ=0.2 uses `transform_tree_pagel`.
- E1i3 vs E1i1 shows refinement is not monotonically good — more iterations helped at
  λ=0.2 and hurt at λ=1.0. Any refinement fix needs a stopping rule, not a fixed count.
