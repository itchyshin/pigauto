# Checkpoint — GNN evidence Phase A campaign

**Date:** 2026-08-26  
**Lane:** `evidence/gnn-sentinel-prerun`  
**Worktree:** `~/local-scratch/lanes/pigauto-gnn-sentinel-prerun`  
**Candidate SHA:** `6fddd79`  
**Branch HEAD:** `709ea58` (F2 confirm scripts)

## STATE

**G0 Phase A APPROVED** (Shinichi, 2026-08-26).  
**Overnight autonomy:** approved until ~05:00 2026-08-27.

| Stage | Status |
|---|---|
| Sentinel pre-run (12 fits) | **DONE** — 0 failures, wall ~3.2 min |
| Pre-registration | **DONE** — `docs/dev-log/2026-08-26-gnn-evidence-preregistration.md` |
| Phase A primary (`fixed_1`) | **DONE** — 2430/2430 RDS, 0 failures, wall 4736 s (~79 min), G8 PASS |
| Phase A analysis | **DONE** — `docs/dev-log/2026-08-26-gnn-evidence-phase-a-results.md` |
| Bayes sensitivity (λ∈{0.2,0.5}) | **DONE** — 1620/1620 RDS, 0 failures, wall 2302 s (~38 min), `results_bayes/` |
| F2 @ λ=1 60-seed confirm | **IN FLIGHT** — 300 fits, `results_confirm/`, launcher PID ~441886 |

## PHASE A RESULTS (primary arm)

| Metric | Value |
|---|---|
| Fits | 2430 / 2430 |
| Failure rate | 0.0% (G8 PASS) |
| Wall (Totoro) | 4736 s (~79 min, < 2 h ceiling) |
| G1–G3, G5–G8 | PASS |
| F1 @ λ=1 specificity | PASS — no G4 explore passes |
| F2 @ λ=1 G4 explore | **5 / 9 cells PASS** — 60-seed confirm launched |
| F3 | Descriptive only (per prereg G7) |

**G4 explore passes (F2 @ λ=1, 30-seed):** n=300 @ 10/30/50%; n=1000 @ 10/30%.  
**Excluded from confirm:** n=1000 @ 50% (rel improve 1.16% < 2% G4 threshold).  
**Near-miss:** n=100 @ 10% (z=2.98), n=1000 @ 50% (rel improve 1.16%).

## F2 CONFIRM SPEC (prereg §3.4)

| confirm_cell | explore_cell_id | n | miss | seeds |
|---:|---:|---:|---:|---:|
| 0 | 42 | 300 | 10% | 60 (60000-block) |
| 1 | 43 | 300 | 30% | 60 |
| 2 | 44 | 300 | 50% | 60 |
| 3 | 51 | 1000 | 10% | 60 |
| 4 | 52 | 1000 | 30% | 60 |

**Total:** 5 cells × 60 seeds = **300 fits**, `lambda_mode=fixed_1`.

## ARTIFACTS

| Path | Contents |
|---|---|
| `script/returned_gnn_campaign/results/` | 2430 primary job RDS + summaries (local pull) |
| Totoro `~/pigauto_gnn_evidence_campaign/results_bayes/` | 1620 bayes sensitivity RDS (complete) |
| Totoro `~/pigauto_gnn_evidence_campaign/results_confirm/` | F2 confirm RDS (in flight) |
| `docs/dev-log/2026-08-26-gnn-evidence-phase-a-results.md` | Phase A results report |

## TOTORO JOBS (2026-08-26)

| Job | Launcher | Fits | Status | Wall |
|---|---|---:|---|---:|
| Bayes sensitivity | `gnn_evidence_sensitivity_bayes_totoro.sh` | 1620 | **DONE** 0 fail | 2302 s |
| F2 confirm | `gnn_evidence_f2_confirm_totoro.sh` | 300 | **RUNNING** | ETA ~10–15 min |

**Poll commands:**

```bash
ssh totoro 'bash ~/pigauto_gnn_evidence_campaign/script/gnn_evidence_f2_confirm_totoro.sh status'
ssh totoro 'bash ~/pigauto_gnn_evidence_campaign/script/gnn_evidence_sensitivity_bayes_totoro.sh status'
```

**Pull when complete:**

```bash
cd ~/local-scratch/lanes/pigauto-gnn-sentinel-prerun
bash script/rsync_gnn_evidence_campaign.sh pull
```

## NEXT GATES

1. **F2 confirm complete** → collect `results_confirm/`, run G4 on 60-seed cells.
2. **Bayes analysis** → collect `results_bayes/`, compare low-λ `gnn_res` per prereg §4.
3. **Phase B (MAR/MNAR)** — requires new G0; do not launch.
4. **PR #174 / product lane** — out of scope.

## RESUME

```
PLATFORM: cursor | ON BRANCH: evidence/gnn-sentinel-prerun @ 709ea58 | LANE: gnn-evidence
Phase A DONE. Bayes sensitivity DONE (1620/1620). F2 confirm IN FLIGHT (poll status).
After F2 confirm: pull → G4 confirm audit → bayes low-λ comparison doc.
```
