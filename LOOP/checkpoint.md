# Checkpoint — GNN evidence campaign (Phase A + Phase B)

**Date:** 2026-08-27 (Phase B launch)  
**Lane:** `evidence/gnn-sentinel-prerun`  
**Worktree:** `~/local-scratch/lanes/pigauto-gnn-sentinel-prerun`  
**Candidate SHA:** `6fddd79`  
**Branch HEAD:** see RESUME block

## STATE

**G0 Phase A APPROVED** (Shinichi, 2026-08-26).  
**G0 Phase B APPROVED** (Shinichi, 2026-08-27, tasks 1+2+3).

| Stage | Status |
|---|---|
| Sentinel pre-run (12 fits) | **DONE** — 0 failures, wall ~3.2 min |
| Pre-registration (Phase A) | **DONE** — `docs/dev-log/2026-08-26-gnn-evidence-preregistration.md` |
| Phase A primary (`fixed_1`) | **DONE** — 2430/2430 RDS, 0 failures, wall 4736 s (~79 min), G8 PASS |
| Phase A analysis | **DONE** — `docs/dev-log/2026-08-26-gnn-evidence-phase-a-results.md` |
| Bayes sensitivity (λ∈{0.2,0.5}) | **DONE** — 1620/1620 RDS, 0 failures |
| F2 @ λ=1 60-seed confirm | **DONE** — 300/300 RDS; **G4 confirm 3/5 PASS** |
| Phase B prereg addendum | **DONE** — `docs/dev-log/2026-08-27-gnn-evidence-phase-b-preregistration.md` |
| Phase B driver + launcher | **DONE** — see ARTIFACTS |
| Phase B primary (`fixed_1`) | **LAUNCHED** — 7290 fits target on Totoro |

## PHASE A RESULTS (primary arm)

| Metric | Value |
|---|---|
| Fits | 2430 / 2430 |
| Failure rate | 0.0% (G8 PASS) |
| Wall (Totoro) | 4736 s (~79 min, < 2 h ceiling) |
| G1–G3, G5–G8 | PASS |
| F1 @ λ=1 specificity | PASS — no G4 explore passes |
| F2 @ λ=1 G4 explore | **5 / 9 cells PASS** (30-seed screen) |
| F2 @ λ=1 G4 confirm | **3 / 5 cells PASS** (60-seed) |
| F3 | Descriptive only (per prereg G7) |
| Bayes low-λ closure (F2) | **0% closed** — gnn_res survives on all 5 fixed_1 G4 cells |

**G4 confirm PASS (manuscript-eligible):** n=300 @ 10%/30%; n=1000 @ 10%.  
**G4 confirm FAIL:** n=300 @ 50% (z=2.39, rel 1.32%); n=1000 @ 30% (z=7.91 but rel 1.84% < 2%).  
**Explore-only (not confirmed):** n=100 @ 10% (z=2.98 explore), n=1000 @ 50% (rel 1.16% explore).

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
| `script/returned_gnn_campaign/results_bayes/` | 1620 bayes RDS + closure tables (local pull) |
| `script/returned_gnn_campaign/results_confirm/` | 300 confirm RDS + G4 cell summary (local pull) |
| `script/collect_gnn_evidence_f2_confirm.R` | G4 confirm collector |
| `script/collect_gnn_evidence_bayes_sensitivity.R` | Bayes closure collector |
| `docs/dev-log/2026-08-26-gnn-evidence-phase-a-results.md` | Phase A + confirm + bayes report (gitignored; local copy) |
| `script/returned_gnn_campaign/PHASE_A_SUMMARY.md` | Durable committed copy of phase-a-results report |
| `script/gnn_evidence_campaign_phase_b.R` | Phase B driver (phylo_MAR / covariate_MAR / MNAR) |
| `script/gnn_evidence_campaign_phase_b_totoro.sh` | Phase B Totoro launcher (7290 fits) |
| `script/rsync_gnn_evidence_phase_b.sh` | Phase B push/pull |
| `script/collect_gnn_evidence_phase_b.R` | Phase B collector + G6 audit |

## TOTORO JOBS

| Job | Launcher | Fits | Status | Wall |
|---|---|---:|---|---:|
| Phase A primary | `gnn_evidence_campaign_totoro.sh` | 2430 | **DONE** 0 fail | 4736 s |
| Bayes sensitivity | `gnn_evidence_sensitivity_bayes_totoro.sh` | 1620 | **DONE** 0 fail | 2302 s |
| F2 confirm | `gnn_evidence_f2_confirm_totoro.sh` | 300 | **DONE** 0 fail | 693 s |
| **Phase B mechanisms** | `gnn_evidence_campaign_phase_b_totoro.sh` | 7290 | **RUNNING** | est ~4 h |

**Phase B poll:**

```bash
ssh totoro 'bash ~/pigauto_gnn_evidence_phase_b/script/gnn_evidence_campaign_phase_b_totoro.sh status'
```

**Phase B pull:**

```bash
cd ~/local-scratch/lanes/pigauto-gnn-sentinel-prerun
bash script/rsync_gnn_evidence_phase_b.sh pull
Rscript script/collect_gnn_evidence_phase_b.R
```

## PHASE B SCOPE

| Quantity | Value |
|---|---|
| Cells | 243 (81 base × 3 mechanisms) |
| Fits | 7,290 |
| Mechanisms | phylo_MAR, covariate_MAR (G6 MAR), MNAR |
| Wall estimate | ~237 min linear (79 min × 3); G8 ceiling ≤ 5 h |

## TOTORO JOBS (2026-08-26, historical)

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

| Track | Status |
|---|---|
| **Manuscript** | **DEFERRED** — scratch at `script/returned_gnn_campaign/MANUSCRIPT_DRAFT.md` (banner: pending AVONET + Phase B). Do not use for submission. |
| **AVONET panel** | **ACTIVE** — real-data corroboration; parallel agent lane |
| **Phase B (MAR/MNAR)** | **RUNNING** — Totoro 7290-fit campaign; no MAR prose until complete |
| **PR #174 / product lane** | Out of scope for evidence lane |

## MANUSCRIPT CLAIM FENCE (Phase A, MCAR only)

Under candidate `6fddd79`, paired same-fit estimand, MCAR missingness, F2 nonlinear DGP at Pagel λ=1.0: the calibrated GNN blend beats the fixed-λ=1 phylogenetic baseline on held-out test cells in **3 of 5** confirmatory cells (60 seeds each, G4 prereg thresholds). Confirmed regimes: n=300 at 10% and 30% missing; n=1000 at 10% missing (2.9–4.8% relative RMSE improvement, |Δ|/MCSE ≥ 4.5). F1 specificity control shows no systematic false-positive GNN wins. Low-λ sensitivity (λ_DGP ∈ {0.2, 0.5}, `lambda_mode = "bayes"`) shows **0% closure** on F2 cells that passed G4 under fixed_1 — incremental GNN value survives a λ-aware baseline. **Not claimed:** n=300 @ 50% or n=1000 @ 30% (confirm fail), n=100 cells, λ<1 primary-arm claims, F3 mixed-type, MAR/MNAR.

## RESUME

```
PLATFORM: cursor | ON BRANCH: evidence/gnn-sentinel-prerun | LANE: gnn-evidence-phase-b
OTHER LANES: codex PR#174 (do not touch)
Phase A COMPLETE. Phase B LAUNCHED (7290 fits, mechanism axis).
Next: poll Totoro → collect → Phase B gate audit. No MAR claims until done.
```
