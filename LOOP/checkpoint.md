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
| Phase B lane 3b (MNAR) | **RUNNING** — 284/2430 RDS, 0 failures @ 2026-08-27 05:45 MDT |
| Phase B lane 3a (phylo/cov MAR) | **RUNNING** — 188/4860 RDS, 0 failures @ 2026-08-27 05:45 MDT |
| AVONET panel | **RUNNING** — 5/15 RDS, **5 failures** (BACE unavailable; Rphylopars OK) |

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
| `script/gnn_evidence_campaign_phase_b.R` | Phase B driver; `PIGAUTO_MECHANISM_ARM` hook |
| `script/gnn_evidence_phase_b_mnar_totoro.sh` | Lane 3b MNAR Totoro launcher (2430 fits) |
| `script/rsync_gnn_evidence_phase_b_mnar.sh` | Lane 3b push/pull |
| `script/collect_gnn_evidence_phase_b_mnar.R` | Lane 3b collector + G6 audit |
| `script/gnn_evidence_campaign_phase_b_totoro.sh` | Monolithic Phase B launcher (7290 fits; optional) |
| `script/rsync_gnn_evidence_phase_b.sh` | Monolithic push/pull |
| `script/collect_gnn_evidence_phase_b.R` | Monolithic collector + G6 audit |

## TOTORO JOBS (coordinator poll — 2026-08-27 05:45 MDT)

**Coordination:** AVONET + Phase B 3a + 3b run **in parallel** on Totoro (384 cores). Worker budget ≈205 (5 + 100 + 100). Load ~91 at poll — within capacity.

| Job | Launcher | Fits done | Failures | Wall | Status |
|---|---|---:|---:|---:|---|
| Phase A primary (legacy) | `gnn_evidence_campaign_totoro.sh` | 2430/2430 | 0 | 4736 s | **DONE** — pulled |
| Bayes sensitivity (legacy) | `gnn_evidence_sensitivity_bayes_totoro.sh` | 1620/1620 | 0 | 2302 s | **DONE** — pulled |
| F2 confirm (legacy) | `gnn_evidence_f2_confirm_totoro.sh` | 300/300 | 0 | 693 s | **DONE** — pulled |
| **Phase B 3a** phylo+cov MAR | `gnn_evidence_phase_b_phylo_cov_mar_totoro.sh` | 188/4860 | 0 | — | **RUNNING** (launched 05:38; ETA ~08:16 MDT) |
| **Phase B 3b** MNAR | `gnn_evidence_phase_b_mnar_totoro.sh` | 284/2430 | 0 | — | **RUNNING** (launched 05:35; ETA ~06:54 MDT) |
| **AVONET panel** | `gnn_evidence_avonet_panel_totoro.sh` | 5/15 | 5 | — | **RUNNING** — investigate failures; BACE=FALSE |

**Poll commands (Totoro monitor agent inactive — poll manually or relaunch):**

```bash
ssh totoro 'bash ~/pigauto_gnn_evidence_phase_b_mnar/script/gnn_evidence_phase_b_mnar_totoro.sh status'
ssh totoro 'bash ~/pigauto_gnn_evidence_phase_b/script/gnn_evidence_phase_b_phylo_cov_mar_totoro.sh status'
ssh totoro 'bash ~/pigauto_gnn_evidence_campaign/script/gnn_evidence_avonet_panel_totoro.sh status'
```

**ETA (Phase A anchor: 4736 s / 2430 fits @ 100 workers):** lane 3b ~79 min from 05:35 → ~06:54 MDT; lane 3a ~158 min from 05:38 → ~08:16 MDT. AVONET pending bootstrap + launch.

**Lane 3b pull:**

```bash
cd ~/local-scratch/lanes/pigauto-gnn-sentinel-prerun
bash script/rsync_gnn_evidence_phase_b_mnar.sh pull
Rscript script/collect_gnn_evidence_phase_b_mnar.R
```

## PHASE B SCOPE (lane 3b MNAR)

| Quantity | Value |
|---|---|
| Cells | 81 (global cell_id 162–242) |
| Fits | 2,430 |
| Mechanism | MNAR only (G6: **not MAR**) |
| MNAR generator | `P(miss) ∝ plogis(1.5·scale(y))`, calibrated to `miss_frac` |
| Wall estimate | ~79 min (Phase A anchor); G8 ceiling ≤ 1.5 h |
| Full Phase B grid | 243 cells / 7,290 fits (3a+3b parallel) |

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
| **Manuscript** | **DEFERRED pending AVONET + Phase B** — scratch `MANUSCRIPT_DRAFT.md` (`1280cb0`). Do not expand or submit. |
| **AVONET panel** | **ACTIVE** — Lane 2; 5/15 RDS, 5 failures (pigauto-only until BACE bootstrap) |
| **Phase B lane 3b (MNAR)** | **ACTIVE** — Lane 3b; 284/2430; ETA ~06:54 MDT |
| **Phase B lane 3a (MAR arms)** | **ACTIVE** — Lane 3a; 188/4860; ETA ~08:16 MDT |
| **PR #174 / product lane** | Out of scope for evidence lane |

## MANUSCRIPT CLAIM FENCE (Phase A, MCAR only)

Under candidate `6fddd79`, paired same-fit estimand, MCAR missingness, F2 nonlinear DGP at Pagel λ=1.0: the calibrated GNN blend beats the fixed-λ=1 phylogenetic baseline on held-out test cells in **3 of 5** confirmatory cells (60 seeds each, G4 prereg thresholds). Confirmed regimes: n=300 at 10% and 30% missing; n=1000 at 10% missing (2.9–4.8% relative RMSE improvement, |Δ|/MCSE ≥ 4.5). F1 specificity control shows no systematic false-positive GNN wins. Low-λ sensitivity (λ_DGP ∈ {0.2, 0.5}, `lambda_mode = "bayes"`) shows **0% closure** on F2 cells that passed G4 under fixed_1 — incremental GNN value survives a λ-aware baseline. **Not claimed:** n=300 @ 50% or n=1000 @ 30% (confirm fail), n=100 cells, λ<1 primary-arm claims, F3 mixed-type, MAR/MNAR.

## RESUME

```
PLATFORM: cursor | ON BRANCH: evidence/gnn-sentinel-prerun @ 1280cb0 | LANE: gnn-evidence
OTHER LANES: codex PR#174 (do not touch)
Phase A COMPLETE (G4 confirm 3/5). Manuscript DEFERRED pending AVONET + Phase B.
Active lanes 2+3: AVONET 5/15 (5 fail); Phase B 3a 188/4860; 3b 284/2430 (0 fail). Unified launcher NOT running (3a+3b only).
Totoro parallel (~205 workers). Next: poll → pull → collect. No manuscript work.
```
