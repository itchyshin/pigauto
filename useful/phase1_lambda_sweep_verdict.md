# λ-threshold sweep verdict

> 2026-04-27 ~20:40. Source: bench_lambda_sweep.rds (48 cells × 6
> methods × 3 reps = 288 rows, 63 min wall).
> Fixed: n=500, ncov=10, β=1.0, sp_miss=0.5, within_miss=0.2.
> Sweep: λ ∈ {0, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50};
> f_type ∈ {nonlinear, interactive}.

## Headline

**Pigauto wins lm_nonlinear above a clean phylogenetic-signal
threshold:**
- **Nonlinear DGP**: pigauto first wins at λ = 0.15 (ratio 0.907,
  9 % better)
- **Interactive DGP**: pigauto first wins at λ = 0.30 (ratio 0.929,
  7 % better)

Below those thresholds, lm_nonlinear (poly+pairwise OLS, no phylo)
wins because it has either the exactly-correct functional form
(interactive) or close to it (nonlinear), and there's no phylo
signal for pigauto to exploit.

## Full table

| λ | nonlinear: pig/lmNL | interactive: pig/lmNL |
|---|---|---|
| 0.000 | 1.506 (lm_NL +51%) | 1.992 (lm_NL +99%) |
| 0.025 | 1.324 (lm_NL +32%) | 1.775 (lm_NL +78%) |
| 0.050 | 1.226 (lm_NL +23%) | 1.455 (lm_NL +46%) |
| 0.100 | 1.015 (tied) | 1.301 (lm_NL +30%) |
| **0.150** | **0.907 (pigauto +9%)** | 1.159 (lm_NL +16%) |
| 0.200 | 0.837 (pigauto +16%) | 0.986 (tied) |
| 0.300 | 0.825 (pigauto +18%) | **0.929 (pigauto +7%)** |
| 0.500 | 0.549 (pigauto +45%) | 0.815 (pigauto +18%) |

## Why interactive's threshold is higher

`lm_nonlinear` uses `poly(z, 2, raw=TRUE) + pairwise interactions`.

- The **interactive** DGP is *exactly* `Z[,1]*Z[,2] + 0.5*Z[,3]^2 +
  0.3*Z[,1]*Z[,4] + ...` — i.e., the saturated correct model. lm_NL
  has zero misspecification bias.
- The **nonlinear** DGP is `sin(Z[,1]) * exp(0.3*Z[,2]) +
  0.5*tanh(Z[,3]) + ...` — fundamentally non-polynomial. lm_NL is
  always misspecified.

Pigauto needs more phylo signal to overcome a saturated correct OLS
than to overcome a misspecified one. Makes sense.

## Why this is a strong paper result

1. **The pattern is monotone and clean.** Pigauto's RMSE relative to
   lm_nonlinear improves smoothly with λ. No noise; the curve is
   crisp.
2. **Three methods are unambiguously beaten across the table:**
   `column_mean`, `species_mean`, and `lm`. Pigauto beats all three
   at every λ (RMSE comparisons in the source).
3. **Two methods compete differently with phylogenetic signal:**
   - `phylolm-λ`: linear with phylo. Pigauto beats it everywhere.
   - `lm_nonlinear`: nonlinear without phylo. Wins at low λ; loses
     at moderate-to-high λ.
   - **Pigauto's value-add IS the joint exploitation of phylogeny +
     nonlinear feature relationships.** Neither lm_NL nor phylolm-λ
     can do both.
4. **Empirical λ ≥ 0.15 is realistic for many ecological datasets.**
   Real comparative datasets typically show measurable Pagel's λ in
   the 0.3–0.9 range for body-size-related traits. So pigauto's
   regime of advantage covers a meaningful slice of the literature.

## Paper figure suggestion

A two-panel plot, one per f_type:
- x-axis: λ (log-ish or linear, 0 to 0.5)
- y-axis: RMSE
- Lines: lm_nonlinear (no phylo, OLS), phylolm-λ (linear, phylo),
  pigauto (nonlinear, phylo)
- Shaded "pigauto-wins" region above the threshold
- Annotation of thresholds at λ=0.15 and λ=0.30

This is the headline figure for the paper. It tells the entire story
in one plot.

## What we still need to confirm

1. **n-scaling**: does the threshold shift with n? If pigauto wins at
   smaller λ when n is larger, that's an important finding (more data
   → AE wins more easily). Launching now: n ∈ {200, 500, 1000, 2000}
   at λ=0.20 (just above the nonlinear threshold).
2. **β-scaling**: does the threshold shift with covariate-signal
   strength? Probably yes — weaker signal = harder for everyone, but
   pigauto's relative advantage may grow.
3. **Real-data confirmation**: re-run a real ecological dataset with
   non-trivial λ to verify pigauto wins there too.

## Decision

This is enough for an honest paper:
- Position pigauto as a *phylogenetic AE* that wins when phylo signal
  is meaningful (λ ≥ 0.15) and feature relationships are nonlinear.
- Show the λ-threshold figure as the headline.
- Acknowledge the failure mode: at λ ≈ 0 with saturated polynomial
  truth, lm_nonlinear is the right tool.

Tonight: launch n-scaling sweep to confirm the threshold is stable
across n. ~2 hours. Then synthesis.
