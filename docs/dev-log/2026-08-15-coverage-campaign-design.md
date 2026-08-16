# Conformal-coverage campaign design (Tier 1) — ADEMP

Companion to `2026-08-15-regime-map-design.md`; same anchors (**Morris et al. 2019 ADEMP**;
**Williams et al. 2024** 11 items, self-audit at end). **Design only — D-139 gates
execution.** Pipeline under evaluation: `main` ≥ `c655d75`.

## A — Aims

- **Primary:** measure empirical coverage of pigauto's conformal prediction intervals
  against the nominal 95%, as a function of the validation-set size that feeds the
  conformal quantile — the axis the theory says matters and the one the only existing
  number (0.884–0.887, leak-tainted, n_val ≈ 10–15) could not separate from noise.
- **Secondary S1:** quantify the `conformal_split_val` trade: coverage and interval width
  under the default `FALSE` (single-set; documented undercoverage risk from gate/conformal
  double-use) vs `TRUE` (split; exchangeability restored). This turns the roxygen caveat
  into a measured recommendation.
- **Secondary S2:** point-accuracy cost of `conformal_split_val = TRUE` (the documented
  2–26% bench regression) re-measured on the corrected pipeline.

## D — Data-generating mechanism

Reuses the regime-map DGP module (families F1/F2 as implemented in
`regime_prerun_cell.R` → campaign scaffold), restricted to where conformal applies
(continuous family). One level: species; traits BM-on-λ-tree, K = 4, Σ exchangeable
r = 0.7; F2's nonlinear links included so intervals are stressed where the GNN
contributes.

**Factors:**

| Factor | Levels | Count |
|---|---|---|
| n species | 100, 300, 1000 | 3 |
| λ | 0.5, 1.0 | 2 |
| family | F1 linear, F2 nonlinear | 2 |
| `conformal_split_val` | FALSE, TRUE | 2 |
| MCAR m | 0.30 (fixed) | 1 |

Cells: 3 × 2 × 2 × 2 = **24**. The n axis is the *effective* n_val axis: with
`val_frac` defaults, n = 100/300/1000 gives roughly n_val ≈ 10/30/100 per trait —
spanning "quantile is noise" to "quantile is stable".

**Replicates:** coverage is the headline → binomial MCSE. Per-rep coverage is itself an
average over that rep's masked cells (~120 at n=100·m=0.3·4 traits, more above), so the
per-cell MCSE is driven by between-rep variation: target MCSE ≤ 1% at 0.95 →
**n_rep = 100** per cell for n ∈ {100, 300}; **n_rep = 40** at n = 1000 (MCSE ≈ 1.6%,
flagged). Fits: 16 cells × 100 + 8 × 40 = **1,920 fits** — same order as the regime map;
the two campaigns share the scaffold and can share a Totoro session back-to-back.

## E — Estimands

Per cell, evaluated at the MCAR-masked cells (truth stored per rep):

- **cover95**: mean over reps of P(truth ∈ [conformal_lower, conformal_upper]).
- **width**: median interval width on the original scale (S1's other half — coverage
  bought by inflation is visible here).
- **Δ_rmse_split**: paired per-rep RMSE difference (split TRUE − FALSE), same seeds (S2).
- Convergence/failure rate; wall per fit.

## M — Methods

pigauto `impute()` defaults ± `conformal_split_val`; the paired baseline BM-SE interval
(`mu ± 1.96·se`) as the model-dependent comparator the docs tell users *not* to rely on —
included precisely to show *why* (its coverage under F2 misspecification is the
motivating contrast). Excluded: external packages (no comparable distribution-free
interval to compare against; Tier-3 owns package comparisons).

## P — Performance measures

- cover95 with binomial MCSE `√(p(1−p)/n_rep)`; decision rule: a cell is "adequate" iff
  its 95% MCSE interval contains 0.95 or lies above it.
- width (median, IQR); Δ_rmse_split with `sd/√n_rep`.
- Failure rate per cell (item 10b), wall per fit.

## Compute plan

Totoro, same conventions as the regime map (100 workers cap, OPENBLAS=1, torch threads 1,
master seed 20260816, per-(cell,rep) RDS, `sessionInfo()` archived, SAFE 3-stage layout).
**Cost**: same per-fit walls as the regime-map pre-run measures *(that pre-run prices this
campaign too — no separate pre-run needed beyond a 2-fit `split_val = TRUE` check, since
the split path has never been timed)*. Provisional: ~300–700 CPU-h → **3–7 h wall at 100
workers**. Stop rules identical to the regime map.

## Reporting

Coverage-vs-n curves per family × λ × split setting (the "minimum n_val guidance" figure
the docs currently cannot cite); width panel; the split-val recommendation stated with
numbers. Full cell table + MCSEs in supplement. Real-data anchor: AVONET300 conformal
coverage (exists in `validate_avonet_full`, re-run post-#158 in Tier 3).

## Williams 11-item self-audit

| # | Item | Where |
|---|---|---|
| 1 | Aims | §A |
| 2 | DGP | §D |
| 3 | Estimands | §E |
| 4 | Methods | §M |
| 5 | Measures + formulas | §P |
| 6 | Software/versions | §Compute |
| 7 | Seeds | §Compute |
| 8 | Code availability | shared scaffold, in-repo |
| 9 | Real-data example | §Reporting |
| 10 | Full results + failures | §P/§Reporting |
| 11 | MCSE everywhere | §P |

**Sequencing:** after the regime-map pre-run prices per-fit walls, both campaigns come
back for one combined go/no-go; if approved they run back-to-back on Totoro.
