# Checkpoint — GNN evidence Phase A campaign

**Date:** 2026-08-26  
**Lane:** `evidence/gnn-sentinel-prerun`  
**Worktree:** `~/local-scratch/lanes/pigauto-gnn-sentinel-prerun`  
**Candidate SHA:** `6fddd79`

## STATE

**G0 Phase A APPROVED** (Shinichi, 2026-08-26).  
**Overnight autonomy:** approved until ~05:00 2026-08-27.

| Stage | Status |
|---|---|
| Sentinel pre-run (12 fits) | **DONE** — 0 failures, wall ~3.2 min |
| Pre-registration | **DONE** — `docs/dev-log/2026-08-26-gnn-evidence-preregistration.md` |
| Phase A primary (`fixed_1`) | **DONE** — 2430/2430 RDS, 0 failures, wall 4736 s (~79 min), G8 PASS |
| Phase A analysis | **DONE** — `docs/dev-log/2026-08-26-gnn-evidence-phase-a-results.md` |
| Bayes sensitivity (λ∈{0.2,0.5}) | **LAUNCHING** — 1620 fits, `results_bayes/` |

## PHASE A RESULTS (primary arm)

| Metric | Value |
|---|---|
| Fits | 2430 / 2430 |
| Failure rate | 0.0% (G8 PASS) |
| Wall (Totoro) | 4736 s (~79 min, < 2 h ceiling) |
| G1–G3, G5–G8 | PASS |
| F1 @ λ=1 specificity | PASS — no G4 explore passes |
| F2 @ λ=1 G4 explore | **5 / 9 cells PASS** — 60-seed confirm warranted for passing cells |
| F3 | Descriptive only (per prereg G7) |

**G4 explore passes (F2 @ λ=1):** n=300 and n=1000 at all miss rates (10/30/50%).  
**Near-miss:** n=100 @ 10% (z=2.98), n=1000 @ 50% (rel improve 1.16%).

## ARTIFACTS

| Path | Contents |
|---|---|
| `script/returned_gnn_campaign/results/` | 2430 job RDS + summaries (local pull) |
| `script/returned_gnn_campaign/results/gnn_campaign_cell_summary.csv` | 81-cell Δ/MCSE/gate table |
| `script/returned_gnn_campaign/results/gnn_campaign_gates.csv` | Full r_cal_gnn per latent col (11340 rows) |
| `docs/dev-log/2026-08-26-gnn-evidence-phase-a-results.md` | Phase A results report |
| Totoro `~/pigauto_gnn_evidence_campaign/results_bayes/` | Bayes sensitivity (in flight) |

## NEXT GATES (for Shinichi @ 05:00)

1. **F2 @ λ=1 60-seed confirm** — 5 cells passed G4 explore; run separate confirmatory arm before any public GNN-positive prose.
2. **Bayes sensitivity** — poll `gnn_evidence_sensitivity_bayes_totoro.sh status`; collect to `results_bayes/`; compare low-λ `gnn_res` per prereg §4.
3. **Phase B (MAR/MNAR)** — requires new G0; do not launch from this checkpoint.
4. **PR #174 / product lane** — out of scope; no edits.

## RESUME @ 05:00

```
PLATFORM: cursor | ON BRANCH: evidence/gnn-sentinel-prerun | LANE: gnn-evidence
Phase A primary DONE (0/2430 fail, 79 min). F2@λ=1: 5/9 G4 explore passes.
Check bayes sensitivity status → collect results_bayes → decide 60-seed F2 confirm cells.
Poll: ssh totoro 'bash ~/pigauto_gnn_evidence_campaign/script/gnn_evidence_sensitivity_bayes_totoro.sh status'
```
