# Known-Σ recovery simulation — design (Tier 2.2, gating the in-house Σ convergence)

2026-08-17 · Claude lane · ADEMP-lite (Morris 2019 / Williams 2024).
Pre-registered BEFORE any estimator comparison is run.

## A — Aims

Decide whether the `e7ca41c` Fisher-ML estimator closes the in-house solver's measured
accuracy gap (`2026-08-16-continuous-gap-diagnosis.md`: +0.14–1.27 z-RMSE vs converged
phylopars REML) **without** reintroducing an external solver — honouring PR #100's
"fully in-house" decision. Secondary: first-ever SE-calibration measurement for the
joint path's per-tip variances (a Tier 2.2 debt independent of the estimator choice).

## D — Data-generating mechanism

`vec(L) ~ MVN(0, Σ ⊗ R)` exactly (matrix-normal BM — the model all three estimators
assume). `R = cov2cor(vcv(rtree(n)))`, λ-transformed. Factors:

- n ∈ {100, 300} × K = 4
- Σ designs: exchangeable r=0.7 (the campaigns' standard); exchangeable r=0.3;
  heterogeneous (variances 0.5–2, off-diagonals mixed sign) — 3 designs
- tree signal λ ∈ {1.0, 0.2} — 1.0 is where the iid-Σ approximation in Fisher-ML is
  WORST (species maximally correlated); 0.2 where it is nearly exact. This axis is the
  design's discriminating contrast.
- MCAR 30% · 50 reps/cell · **72 cells × 3 estimators**. Solver-level fits (no GNN):
  local, minutes total.

## E/M — Estimators and estimands

- **E0** current in-house single-pass (control — the measured-deficient baseline)
- **E1** Fisher-ML (`e7ca41c` port: pairwise-complete init → Cholesky-optim observed-data
  NLL, iid across species → per-column BM init → EM cell refinement)
- **E2** `Rphylopars::phylopars` REML (yardstick only, per #100 not a shipping path)

Per (estimator, cell): Σ̂ error — relative Frobenius `‖Σ̂−Σ‖_F/‖Σ‖_F`, off-diagonal
sign accuracy, mean element bias; tip prediction RMSE on masked cells (z scale);
**SE calibration** — empirical coverage of `anc_recon ± 1.96·√anc_var` on masked cells.
MCSE = sd/√50 throughout.

## P — Pre-registered decision rule

E1 ships iff, in EVERY cell:
1. Frobenius rel. error ≤ 1.25 × E2's, and
2. prediction RMSE within 3×MCSE of E2's, and
3. SE coverage ∈ [0.93, 0.97].

E0 is expected to fail (that is the point); if E0 *passes* everywhere the AVONET gap is
not Σ-driven and the diagnosis needs revisiting. If E1 fails **only** at λ=1.0 (the iid
approximation's worst case), the escalation is a corrected matrix-normal EM M-step
(second moments included) — designed then, not now. Real-data gate afterwards: the
5-seed AVONET bench must reach Mass ≤1.30, Beak ≤0.60, Tarsus ≤0.87, Wing ≤0.45.

## Compute

Local Mac, sequential, minutes (solver-only). No D-139 gate needed (<30 min by
construction; measured at pre-run).
