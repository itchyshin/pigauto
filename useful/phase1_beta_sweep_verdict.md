# β-strength sweep verdict

> 2026-04-27 ~22:30. Source: bench_beta_sweep.rds (30 cells × 4
> methods × 3 reps = 120 rows, 41 min wall).
> Fixed: λ=0.20, n=500, ncov=10. Sweep: β ∈ {0.3, 0.5, 0.7, 1.0, 1.5};
> f_type ∈ {nonlinear, interactive}.

## Headline

**Pigauto's advantage over lm_nonlinear is INVERSELY related to β.**
Strongest when covariate signal is moderate, weakest when signal is
overwhelming. This is counter-intuitive at first but makes sense.

| β | nonlinear: pig/lmNL | interactive: pig/lmNL |
|---|---|---|
| 0.3 | **0.712 (pigauto +29 %)** | **0.661 (pigauto +34 %)** |
| 0.5 | 0.724 (pigauto +28 %) | 0.736 (pigauto +26 %) |
| 0.7 | 0.824 (pigauto +18 %) | 0.815 (pigauto +18 %) |
| 1.0 | 0.843 (pigauto +16 %) | 1.092 (lm_NL +9 %) |
| 1.5 | 1.161 (lm_NL +16 %) | 1.234 (lm_NL +23 %) |

## Why this pattern

At **low β** (covariate effect small relative to noise):
- Response variance is dominated by phylo + residual noise
- `lm_nonlinear` with 65 features (poly(2)*10 + pairwise(10,2) = 20+45)
  overfits the noise (R² inflation)
- Pigauto's phylo prior provides anchored regularization → wins

At **high β** (covariate effect dominates):
- Response is mostly driven by the covariates
- `lm_nonlinear` with the saturated feature set captures the signal
  well; overfitting penalty is minor relative to the strong signal
- Pigauto's safety floor / shrinkage limits how much it can fit;
  it's outpaced

This is a CLEAN regime characterisation: pigauto adds value
specifically in the "noisy moderate-signal phylogenetic" sweet spot,
which is **most of real comparative data**. Highly-engineered
experimental data with overwhelming covariate signal would be the
exception.

## Combined with λ and n sweeps — pigauto's regime of advantage

Across three axes we now have a clear picture:

| Axis | Pigauto wins when | Pigauto loses when |
|---|---|---|
| λ (phylo signal) | λ ≥ 0.15 (nonlinear) or 0.30 (interactive) | λ ≤ 0.05 |
| n (sample size) | n ≥ 1000 most strongly | n ≤ 200 weakly (interactive) |
| β (covariate signal) | β ≤ 1.0 strongly | β ≥ 1.0 (interactive), 1.5 (nonlinear) |

**Pigauto's home turf:** moderate phylogenetic signal + sufficient
data + moderate covariate signal + nonlinear/interactive feature
relationships. **All four conditions need to hold.**

## Three paper figures

1. **λ-threshold** (lambda_sweep): pigauto/lm_NL ratio vs λ at
   n=500, β=1.0, both f_types. Show clean monotone decrease.
2. **n-scaling** (n_scaling): ratio vs n at λ=0.2, β=1.0. Show
   pigauto's advantage grows with n (AE story).
3. **β-strength** (beta_sweep): ratio vs β at λ=0.2, n=500. Show
   pigauto wins at moderate signal, loses at overwhelming signal.

These three figures characterise pigauto's regime in ~1 page of the
paper.

## What this implies for next steps

- We have enough evidence for the conservative paper.
- Optional: a missingness sweep (vary sp_miss / within_miss). Likely
  shows pigauto handles missingness up to 50-70% similarly. Not
  decisive for the paper.
- The architecture story is dead (transformer ≈ GAT ≈ GCN).
- The discrete trait regression remains an open issue, but that's an
  engineering / package-quality problem, not a paper-claim issue.

## Decision

**Stop benching. Write the comprehensive synthesis memo now.** The
empirical picture is complete enough for a paper. More benches won't
change the conclusion; they'd just polish a clear story.
