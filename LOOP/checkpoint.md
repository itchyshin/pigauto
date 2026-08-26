# Checkpoint — GNN evidence Phase A campaign

**Date:** 2026-08-26  
**Lane:** `evidence/gnn-sentinel-prerun`  
**Worktree:** `~/local-scratch/lanes/pigauto-gnn-sentinel-prerun`  
**Candidate SHA:** `6fddd79`

## STATE

**G0 Phase A APPROVED** (Shinichi, 2026-08-26).

| Stage | Status |
|---|---|
| Sentinel pre-run (12 fits) | **DONE** — 0 failures, wall ~3.2 min |
| Pre-registration | **DONE** — `docs/dev-log/2026-08-26-gnn-evidence-preregistration.md` @ `b7e597a` |
| Phase A driver | **DONE** — `script/gnn_evidence_campaign.R` |
| Totoro launcher | **DONE** — `script/gnn_evidence_campaign_totoro.sh` |
| Local smoke (1 fit) | **DONE** — job 0 OK, fit_sec ~200s (laptop) |
| Totoro launch | **IN FLIGHT** — see below |

## CAMPAIGN SCOPE (Phase A)

- **81 cells:** F1/F2/F3 × n∈{100,300,1000} × λ∈{0.2,0.5,1.0} × miss∈{10%,30%,50%} × MCAR
- **2,430 fits:** 30 paired seeds per cell
- **Primary arm:** `lambda_mode = "fixed_1"`
- **Sensitivity:** `lambda_mode = "bayes"` on λ∈{0.2,0.5} — **not in this launch** (separate arm)
- **Host:** Totoro, ≤100 workers, OPENBLAS=1, OMP=1
- **Wall ceiling:** ≤2 h (G8 stop)

## TOTORO LAUNCH

- **Remote dir:** `~/pigauto_gnn_evidence_campaign`
- **Launcher:** `nohup bash script/gnn_evidence_campaign_totoro.sh > logs/campaign.log 2>&1 &`
- **Poll:**
  ```bash
  ssh totoro 'bash ~/pigauto_gnn_evidence_campaign/script/gnn_evidence_campaign_totoro.sh status'
  ssh totoro 'tail -f ~/pigauto_gnn_evidence_campaign/logs/campaign.log'
  ```
- **Collect locally:**
  ```bash
  bash script/rsync_gnn_evidence_campaign.sh pull
  ```

## GATE AUDIT (sentinel pre-run)

| Gate | Status |
|---|---|
| G1 provenance | PASS |
| G2 paired isolation | PASS |
| G3 fallback | PASS |
| G5 no-benefit retention | PASS |
| G6 MCAR labeling | PASS |
| G7 trait boundary | PASS |
| G8 stop (sentinel) | PASS — 0/12 failures |

Phase A campaign G8 monitored at end of Totoro run (>20% failures or wall >2h).

## ARTIFACTS

| Location | Contents |
|---|---|
| Totoro `~/pigauto_gnn_evidence_campaign/results/` | `gnn_campaign_job_*.rds`, summary CSV |
| Totoro `~/pigauto_gnn_evidence_campaign/logs/` | per-job + campaign.log |
| Local `script/returned_gnn_campaign/` | post-pull mirror |

## NEXT

- Monitor Totoro to completion (G8).
- Post-run: cell-level Δ/MCSE tables, gate distributions.
- Sensitivity arm (`bayes` @ low λ): separate G0 or tagged follow-on.
- F2 @ λ=1 confirm (60 seeds): only after G4 pass — not in Phase A launch.
