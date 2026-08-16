# λ-attribution campaign — design (ADEMP-lite)

2026-08-16 · Claude lane · follows Morris et al. 2019 (ADEMP) + Williams et al. 2024
reporting items, in the pattern of `2026-08-15-regime-map-design.md`. Runner:
`lambda_cell.R` on Totoro (`~/pigauto_regime_map/`).

## A — Aims

**Primary:** determine how much of the regime map's low-λ GNN gains
(`2026-08-16-regime-map-results.md`: F1 mean Δ −0.104 at λ=0.2 vs −0.001 at λ=1.0) is
obtainable by fitting Pagel's λ in the *baseline* (`lambda_mode = estimate/cv/bayes`)
instead of routing around the λ=1 misspecification with the GNN. The regime-map results
doc pre-committed: "no claim about the low-λ regime should ship until [this] is run."

**Secondary:** (S1) does the F2 nonlinear GNN lift survive a λ-fitted baseline —
i.e. blend−baseline Δ *within* each λ-fitted mode; (S2) first multi-seed measurement of
the ML boundary pathology (λ̂→0, the sla case) and whether `cv`/`bayes` avoid it —
closing the NEWS "v0.11 CV-λ" promise or documenting that it fails.

## D — Data-generating mechanism

Identical to the regime map (paired comparability): `rtree(n)`; `V_λ = λV + (1−λ)diag(V)`;
K=4 traits, Σ exchangeable r=0.7; families **F1 linear** (specificity control) and
**F2 nonlinear** (`sin(2Z₁)`, `Z₂²·sign(Z₂)`); MCAR 30%; single-obs.
**F3/discrete excluded** — λ modes force the per-column BM path and never reach
LP/threshold/OVR, so there is nothing to attribute there.
**No covariates** — λ is silently dropped under covariates (`fit_baseline.R:546–548`),
so a covariate arm would confound λ-fitting with covariate handling.

Grid: family {F1,F2} × λ_DGP {0.2,0.5,0.8,1.0} × n {100,300} ×
`lambda_mode` {fixed_1, estimate, cv, bayes} × 30 reps = **1,920 jobs**.

**Seeding:** seed = f(family, λ, n, rep) — *mode-invariant*, so all four modes see
identical data and masks; mode comparisons are paired within rep. n=1000 dropped
(regime map showed the λ effect at every n; n=1000 tripled cost without changing the
answer's shape).

## E — Estimands

Per job, on the genuinely-missing cells (truth known): per-trait `loss_blend`,
`loss_base` (MSE, decoded scale); `r_cal_gnn`; floor-fired flag; for λ-fitted modes,
λ̂ per column (recomputed diagnostically on all input-observed cells — the in-fit
estimate additionally masks val/test, so this is approximate and labelled as such;
the fit discards its own λ̂ at `fit_baseline.R:550`).

Derived, per (family, λ_DGP, n) cell:
1. **Baseline repair:** `loss_base(mode) − loss_base(fixed_1)`, paired per rep. How much
   of the misspecification does λ-fitting remove?
2. **Residual GNN value:** `loss_blend(mode) − loss_base(mode)` within each mode. Does
   the GNN still earn its place once the baseline fits λ?
3. **Gap closure fraction:** 1 − [blend(fixed_1) − base(mode)]/[blend(fixed_1) −
   base(fixed_1)] — how much of the GNN's apparent low-λ gain the λ-fitted baseline
   reproduces without any neural network.
4. **Boundary rate:** P(λ̂ < 0.05) per mode (the sla pathology), and λ̂ bias vs λ_DGP.

## M — Methods

One method family (pigauto), four baseline configurations. This is an *attribution*
study, not a competitor bench (that is Slice C). `impute(df, tree, epochs = 500,
lambda_mode = mode)` — production defaults otherwise, same as the regime map.

## P — Performance & MCSE

MCSE = sd(paired Δ)/√30 per cell; detection rule |Δ|/MCSE ≥ 3, as pre-registered for
the regime map. At 30 reps the continuous per-cell MCSE from the regime-map analogue is
≈0.015–0.03, so detectable effects are ≥0.05–0.09 — the λ=0.2 F1 gap (0.104) is well
above; the λ=0.8 gap (0.006) is below, and is negligible regardless.

**Interpretation, pre-registered (before seeing numbers):**
- If λ-fitted baselines close ≥~70% of the low-λ gap → the regime-map low-λ rows are
  re-labelled "obtainable by λ-fitting"; the GNN's accuracy claim narrows to
  F2-nonlinear structure; the robustness claim ("graceful under misspecification")
  stands but is bounded by the cheaper remedy.
- If they close little of it (ML collapsing to boundary at low signal, per the sla
  case) → the GNN's misspecification-absorption is not reproducible by point-λ fitting;
  report per-mode, since `bayes` may behave differently from `estimate`.
- Either way F1-λ=1.0 must stay null for every mode (specificity control); a mode that
  "wins" there is a red flag.
- S2 verdict: `cv`/`bayes` judged on boundary rate + loss at λ_DGP ≤ 0.5 vs `estimate`.

## D-139 gate

Analogue: regime map 5,400 jobs ≈ 2.5 h wall at ~100 workers (n≤300 jobs ≈ 100–150 s).
λ modes add per-fit baseline cost (cv/bayes grid-refit; small vs 500-epoch GNN
training). **Estimate: ~1–1.5 h wall at 100 workers.** Pre-run: 8 jobs (2 cells × 4
modes) timed before launch; results appended below. Overrun → stop and re-report.

## Pre-run results (2026-08-16, D-139)

8/8 OK, 0 failures. F2 n=100: 44.0–46.1 s/job; F1 n=300: 58.6–60.6 s/job. λ modes add
no measurable cost over fixed_1 (bayes 58.6 s vs fixed_1 59.2 s at n=300). Revised
estimate: 1,920 jobs ≈ 28 CPU-h → **~20–25 min wall at 100 workers** — under the D-139
30-minute just-run line. Launched at 100 workers alongside mech_cov at 40 (140 ≤ 150,
D-143).
