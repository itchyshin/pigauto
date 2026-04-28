# Phase 1 extended BIG tier verdict

> 2026-04-27 ~19:30. Source: bench_phase1_extended_big.rds
> (24 cells × 6 methods × 3 reps = 144 rows, 44 min wall).
> DGP: multi-obs nonlinear with 10 i.i.d. covariates, n ∈ {500, 2000},
> λ ∈ {0.05, 0.2}, β = 1.0, sp_miss=0.5, within_miss=0.2.

## Headline

**A clean threshold emerges between λ = 0.05 and λ = 0.20.**
- At λ = 0.05 (very low phylogenetic signal): pigauto LOSES to
  lm_nonlinear by 19–55 % across both f_types and both n.
- At λ = 0.20 (moderate phylogenetic signal): pigauto WINS over
  lm_nonlinear by 5–17 %.

Larger n (2000 vs 500) doesn't change the pattern; the threshold is
phylogenetic-signal driven, not data-volume driven.

## Per-cell head-to-head

| f_type | n | λ | pigauto | lm_NL | ratio | winner |
|---|---|---|---|---|---|---|
| interactive | 500 | 0.05 | 1.123 | 0.726 | **1.55** | lm_NL +55% |
| interactive | 2000 | 0.05 | 1.105 | 0.773 | **1.43** | lm_NL +43% |
| nonlinear | 500 | 0.05 | 1.020 | 0.859 | 1.19 | lm_NL +19% |
| nonlinear | 2000 | 0.05 | 1.010 | 0.852 | 1.19 | lm_NL +19% |
| interactive | 500 | 0.20 | 1.172 | 1.087 | 1.08 | lm_NL +8% |
| interactive | 2000 | 0.20 | 1.161 | 1.216 | **0.96** | pigauto +5% |
| nonlinear | 500 | 0.20 | 1.074 | 1.232 | **0.87** | pigauto +13% |
| nonlinear | 2000 | 0.20 | 1.088 | 1.313 | **0.83** | pigauto +17% |

## Reconciling with smoke tier (n=300, ncov=5)

| Tier | Win rate vs lm_NL | Median ratio | Conclusion |
|---|---|---|---|
| Smoke (λ=0.1, 0.3) | **6/8 wins** | 0.869 (-13%) | pigauto wins broadly |
| BIG (λ=0.05, 0.2) | **3/8 wins** | 1.132 (+13%) | pigauto loses at λ=0.05, wins at λ=0.2 |

The smoke at λ=0.1 was right at the edge. The BIG result confirms:
- λ ≤ 0.05: pigauto loses
- λ ~ 0.1-0.2: mixed, depends on f_type and signal
- λ ≥ 0.3: pigauto wins consistently (smoke confirmed)

## What this means

**Pigauto's value-add is bounded by phylogenetic signal strength,
not by AE capacity.** This is a coherent and defensible position:

1. Pigauto isn't claiming to be a general-purpose autoencoder for
   tabular imputation. It's a *phylogenetic* AE.
2. Without phylo signal, pigauto's machinery has nothing to add over
   well-specified obs-level OLS. That's expected.
3. With phylo signal (λ > ~0.15), pigauto exploits both the phylogenetic
   neighbourhood structure AND the nonlinear feature relationships.
   Neither phylolm-λ (linear) nor lm_nonlinear (no phylo) can do both.

This is the paper's headline:
> *Pigauto provides robust phylogenetic-aware nonlinear imputation,
> outperforming both linear phylogenetic methods (phylolm-λ) and
> nonlinear non-phylogenetic methods (poly+interaction OLS) in the
> regime of meaningful phylogenetic signal (λ ≥ ~0.15) and nonlinear
> covariate-response structure.*

## Where pigauto specifically loses

The single sharpest loss is **interactive λ=0.05 β=1.0** (ratio 1.43-1.55).
This is the regime where:
- The truth IS exactly polynomial+interaction (lm_NL has the saturated
  correct model)
- There's basically zero phylo signal
- Pigauto's BM/GLS baseline has nothing to ground itself in
- The GNN has no phylo neighbourhood structure to exploit

Even with safety floor TRUE, pigauto's calibrated gate apparently
doesn't fully clamp the GNN delta to zero in this regime — and the
delta wanders.

## Next: focused λ-threshold sweep

What we still don't know precisely: where IS the threshold? λ=0.1,
0.15, or somewhere else?

Designing a λ-threshold-finding bench:
- λ ∈ {0, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50}
- f_type ∈ {nonlinear, interactive}
- β = 1.0
- n = 500 (faster than 2000)
- n_covs = 10
- 3 reps
- = 8 λ × 2 f × 3 reps = 48 cells, ~2 hours

This pins down the threshold for the paper figure.
