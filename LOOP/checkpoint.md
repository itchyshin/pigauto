# Checkpoint — GNN evidence sentinel pre-run

**Date:** 2026-08-26  
**Lane:** `evidence/gnn-sentinel-prerun` @ `dd66b33`  
**Worktree:** `~/local-scratch/lanes/pigauto-gnn-sentinel-prerun`  
**Candidate SHA:** `6fddd79`

## STATE
S4 **DONE** — driver committed (`script/gnn_evidence_sentinel_prerun.R`).  
S5 **DONE** — Totoro 12/12 fits, 0 failures, wall ~3.2 min (PID 69761).  
S6 **DONE** — timing table + revised estimate below.

## TIMING (Totoro, 12 parallel, OPENBLAS=1)

| job | family | n | seed | fit_sec | delta | floor |
|-----|--------|---|------|---------|-------|-------|
| 0 | F1 | 100 | 101 | 38.0 | +0.059 | yes |
| 1 | F2 | 100 | 101 | 37.6 | +0.002 | yes |
| 2 | F3 | 100 | 101 | 39.0 | −0.034 | yes |
| 3 | F1 | 300 | 101 | 49.6 | +0.015 | yes |
| 4 | F2 | 300 | 101 | 49.3 | +0.024 | yes |
| 5 | F3 | 300 | 101 | 50.1 | +0.037 | yes |
| 6 | F1 | 1000 | 101 | 183.0 | −0.008 | no |
| 7 | F2 | 1000 | 101 | 185.8 | −0.011 | no |
| 8 | F3 | 1000 | 101 | 194.3 | 0.000 | no |
| 9 | F2 | 300 | 112 | 50.0 | −0.024 | no |
| 10 | F2 | 300 | 123 | 50.0 | −0.004 | no |
| 11 | F2 | 300 | 134 | 49.9 | +0.003 | no |

Per-n mean fit_sec: n=100 **38.2s**, n=300 **49.8s**, n=1000 **187.7s**.

## REVISED FULL-CAMPAIGN ESTIMATE

- **108-cell × 30-seed grid (3240 fits):** ~1.5 h wall @100 workers (regime-map anchor 5400/2.5h); ~50 min by measured CPU/n alone (likely optimistic).
- **81-cell MCAR-only × 30 seeds (2430 fits):** ~1.1 h @ anchor throughput.
- **Handover 4–8 h band:** still plausible once λ × missingness × mechanism axes expand; n=1000 dominates at ~3 min/fit.

## GATE AUDIT (G1–G8)

| Gate | Status | Evidence |
|------|--------|----------|
| G1 provenance | PASS | All 12 RDS name candidate/driver/host/seed/family/n/miss |
| G2 paired isolation | PASS | blend_loss, baseline_loss, paired_delta, r_cal_* recorded |
| G3 fallback | PASS | floor_fired retained; r_cal_gnn=0 reachable (job 8 F3 n=1000) |
| G4 positive claim | N/A | Pre-run only; no claim gate |
| G5 no-benefit | PASS | Mixed deltas retained (F2 n=1000 negative) |
| G6 missingness | PASS | MCAR only; labeled MCAR |
| G7 trait boundary | PASS | F1/F2/F3 reported separately |
| G8 stop | PASS | 0/12 failures; all fields present; wall 3.2 min < approval |

## ARTIFACTS

- Totoro: `~/pigauto_gnn_sentinel_prerun/{results,logs}/`
- Local: `script/returned_gnn_sentinel/`

## NEXT

Return deliverables to parent. No push/merge. Optional: commit F3 trait_types fix (`dd66b33+`).
