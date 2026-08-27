# GNN evidence campaign — design companion

**Date:** 2026-08-26  
**Pre-registration:** [`2026-08-26-gnn-evidence-preregistration.md`](2026-08-26-gnn-evidence-preregistration.md)  
**Regime-map lineage:** `af2e80a` → [`2026-08-15-regime-map-design.md`](2026-08-15-regime-map-design.md) (issue #135)

Short companion to the pre-registration. The prereg is authoritative for gates,
estimands, and approval; this note records design rationale and driver pointers.

---

## Why a new campaign under `6fddd79`

The 2026-08-16 regime map (`5677b2a`, pipeline `c655d75`) established the paired
estimand and F1/F2/F3 discrimination. Lambda-attribution (`a64671e`) resolved
low-λ interpretation. This campaign re-runs the grid on the **0.11 candidate**
(`6fddd79`) so evidence aligns with the release surface under review (PR #174
product lane — not modified here).

Changes from the original regime-map design:

| Aspect | Regime map (`af2e80a`) | This campaign |
|---|---|---|
| Candidate | `c655d75` | `6fddd79` |
| λ levels | 0.2, 0.5, 0.8, 1.0 | 0.2, 0.5, 1.0 |
| Mechanisms | MCAR only | Phase A MCAR; Phase B adds 3 |
| Reps | 60 (n≤300), 30 (n=1000) | 30 explore; 60 confirm (G4 cells) |
| `lambda_mode` | not varied | primary `fixed_1`; sensitivity `bayes` at low λ |
| Pre-run | 6-fit D-139 test | 12-fit sentinel (complete) |

---

## Driver path

| Stage | Script | Status |
|---|---|---|
| Sentinel pre-run (G0) | `script/gnn_evidence_sentinel_prerun.R` | **Complete** — 12/12 fits |
| Phase A campaign | `script/gnn_evidence_campaign.R` | **Pending** — write after prereg commit + G0 Phase A approval |
| Phase B extension | TBD (mechanism generators) | Blocked on Phase A + separate approval |

Sentinel driver pins `CANDIDATE_SHA <- "6fddd79"`, uses regime-map DGP helpers
(F1/F2/F3), records paired Δ and all three `r_cal_*` weights. Phase A driver
will extend the factorial grid (λ × missingness axes) while preserving the same
receipt schema.

---

## Phased rollout rationale

1. **Phase A (MCAR + λ + missingness)** — 81 cells, 2,430 fits. Answers the
   core "where does the GNN buy?" question without mechanism-generator risk.
   Budget-conscious: ~1–1.5 h @ 100 workers (measured pre-run anchor).

2. **Phase B (mechanisms)** — adds phylo-MAR, genuine covariate-MAR, MNAR.
   Gated separately because mechanism DGPs need validation (prior sweep had
   degenerate trait-MAR) and because **no MAR sentence is permitted before
   Phase B completes**.

---

## Priority cell: F2 @ λ = 1

The historically defensible niche:

- At λ = 1.0 the BM baseline is correctly specified on the DGP tree.
- F2's nonlinear cross-trait links (`sin(2·Z₁)`, `Z₂²·sign(Z₂)`) are structurally
  beyond a linear joint-MVN baseline.
- Regime map: F2 mean |Δ|/MCSE = 3.29 vs F1 = 1.51 at λ = 1.0 (`5677b2a`).
- Lambda-attribution: F2 @ λ = 1, n = 300, `gnn_res = −0.0169`, |Δ|/MCSE = 6.5
  on top of `lambda_mode = "bayes"` baseline (`a64671e`).

Protocol: 30-seed explore (embedded in Phase A) → G4 screen → 60-seed confirm.

---

## Compute anchors

From sentinel pre-run (`LOOP/checkpoint.md`):

```
n=100:  38.2 s/fit
n=300:  49.8 s/fit
n=1000: 187.7 s/fit
```

Phase A: 81 cells × 30 seeds = 2,430 fits → ~1.0–1.5 h wall @ 100 workers.

---

## Related artifacts

- Handover: `docs/dev-log/handover/2026-08-26-cursor-handover.md`
- Sentinel checkpoint: `LOOP/checkpoint.md`
- Sentinel results (local): `script/returned_gnn_sentinel/`
- Totoro results: `~/pigauto_gnn_sentinel_prerun/{results,logs}/`
