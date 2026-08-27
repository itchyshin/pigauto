# GNN evidence Phase A — Sys Biol manuscript draft

**Date:** 2026-08-27  
**Lane:** `evidence/gnn-sentinel-prerun` @ `5d41452`  
**Status:** Internal draft only — not for public posting  
**Source:** `script/returned_gnn_campaign/PHASE_A_SUMMARY.md`, confirm/bayes CSVs  
**Canonical path:** `script/returned_gnn_campaign/MANUSCRIPT_DRAFT.md` (`docs/` gitignored)

---

## Results (F2 @ λ=1, G4-confirmed cells only)

Under MCAR missingness and a nonlinear continuous DGP (family F2) with Pagel λ fixed at 1.0 on both sides of the comparison, the calibrated GNN blend `(1 − r_cal) × baseline + r_cal × δ` beat the fixed-λ=1 phylogenetic Brownian baseline on held-out test cells in **three** of five confirmatory regimes re-run at 60 seeds each (G4 preregistered thresholds: F2 family, Δ < 0, |Δ|/MCSE ≥ 3, relative RMSE improvement ≥ 2%). The confirmed cells were n = 300 at 10% missing (mean Δ = −0.028, MCSE = 0.003, |Δ|/MCSE = 8.8, 4.8% relative improvement), n = 300 at 30% missing (Δ = −0.018, MCSE = 0.004, |Δ|/MCSE = 4.5, 2.9%), and n = 1000 at 10% missing (Δ = −0.018, MCSE = 0.001, |Δ|/MCSE = 13.4, 3.3%). In all three, the calibrated gate remained materially open (fraction gate open 0.92–1.00; mean r_cal ≈ 0.40–0.47), indicating that the GNN contribution was active rather than falling back to baseline-only prediction. The estimand is paired within-fit (Δ = blend loss − baseline loss on the same held-out cells); primary-arm explore used 30 seeds, with 60-seed confirmation required before any positive claim (prereg §3.4). Two explore-flagged cells did not survive confirmation (n = 300 @ 50% missing: |Δ|/MCSE = 2.4, rel improve 1.3%; n = 1000 @ 30%: |Δ|/MCSE = 7.9 but rel improve 1.8%) and are excluded from this sentence.

---

## Specificity (F1 null control)

Family F1 implements a linear Brownian DGP at λ = 1 — the regime where the phylogenetic baseline is correctly specified and the GNN should not systematically win. Across all nine F1 @ λ=1 cells (n ∈ {100, 300, 1000} × missing ∈ {10%, 30%, 50%}, 30 seeds each), **no cell passed G4 explore**. Mean Δ ranged from −0.004 to +0.030 with MCSE 0.0003–0.011; the largest |Δ|/MCSE was 2.6 (n = 100, 50% missing), well below the G4 threshold of 3. Relative improvements were near zero (|rel improve| ≤ 4.3%, mostly < 1%). Although calibrated gates were often open (frac gate open 0.40–0.87), the blend did not produce preregistered positive GNN wins on held-out loss — consistent with the gate-closed fallback safety when the baseline is already adequate. F1 therefore supports specificity of the F2 signal rather than a global GNN advantage.

---

## Bayes sensitivity (0% F2 closure)

To test whether incremental GNN value survives a stronger baseline, we re-ran all low-λ cells (λ_DGP ∈ {0.2, 0.5}) with `lambda_mode = "bayes"`, defining gnn_res = blend − λ-fitted baseline (baseline Pagel λ estimated per replicate rather than fixed at 1). **Closure** is defined as a fixed_1 arm that passed G4 explore but whose bayes arm fails G4 (|Δ|/MCSE < 3 or rel improve < 2%). Among F2 cells, five passed G4 explore under fixed_1 at each of λ = 0.2 and λ = 0.5; **zero were closed** under bayes (0% closure). At λ = 0.2, six of nine F2 cells passed gnn_res under bayes (including one cell that failed fixed_1 explore); at λ = 0.5, five of nine passed — matching or exceeding the fixed_1 count. F1 and F3 showed no G4 fixed_1 passes and no closures. This sensitivity analysis supports the interpretation that GNN incremental value on F2 is not an artefact of baseline λ misspecification at λ = 1, though we do **not** extend primary positive claims to λ < 1 on the basis of this arm alone.

---

## Limitations and claim fences

**Scope of positive claim.** Manuscript-eligible GNN-positive language is restricted to the three G4-confirmed F2 @ λ=1 cells listed above. Do not pool wins across cells, families, or λ values into a single headline rate.

**Not claimed — sample size.** n = 100 cells were explore-only; none reached G4 confirm. No positive claim at n = 100.

**Not claimed — failed confirm cells.** n = 300 @ 50% and n = 1000 @ 30% failed G4 confirm despite passing explore; cite as sensitivity/context only, not as evidence of benefit.

**Not claimed — missingness mechanism.** Phase A is MCAR only. MAR, MNAR, and missing-not-at-random phylogenetic missingness require Phase B (separate preregistration and G0).

**Not claimed — mixed types.** F3 (eight-trait mixed-type panel) is descriptive only per prereg G7; no G4 passes and no inferential claim.

**Not claimed — method superiority.** These results do not establish pigauto superiority over BACE, phylopars-only pipelines, or other phylogenetic imputers. Comparisons here are internal: gated GNN blend vs fixed-λ baseline within the same fit.

**Not claimed — λ < 1 primary arm.** Bayes sensitivity is supportive context; primary confirmatory arm is fixed_1 at λ_DGP = 1.

**Not claimed — explore-only cells.** Five F2 @ λ=1 explore passes at 30 seeds; only three survived 60-seed confirm. Explore flags are screening, not publication claims.

**Candidate version.** Fits under package candidate at evidence-lane SHA; product release claims await merge outside this lane (PR #174 out of scope).

---

## Figure specifications

### Figure 1 — Primary-arm Δ heatmap (explore, 30 seeds)

**Purpose:** Show full Phase A landscape; positive signal localized to F2 @ λ=1.

**Data:** `script/returned_gnn_campaign/results/gnn_campaign_cell_summary.csv` — `mean_delta` (or `rel_improve`) by cell.

**Layout:** Facet or panel rows = family (F1, F2, F3); columns = phylo_lambda (0.2, 0.5, 1.0); within each panel, heatmap axes = n_species (100, 300, 1000) × miss_frac (10%, 30%, 50%). Color = mean Δ (blue = blend better) or diverging centered at 0.

**Annotations:** Mark G4 explore passes (F2 @ λ=1 only: 5 cells). Do not mark confirm passes here (explore stage).

**Caption fence:** "Explore stage, 30 seeds; G4 confirm required before positive claim. MCAR, fixed_1 arm."

---

### Figure 2 — F2 @ λ=1 confirm: Δ with MCSE bars (5 cells)

**Purpose:** Confirmatory evidence and explicit failures.

**Data:** `script/returned_gnn_campaign/results_confirm/gnn_confirm_cell_summary.csv`.

**Layout:** Five bars (x = cell label: n × miss), y = mean Δ, error bars = ± MCSE. Horizontal reference lines: Δ = 0; optional secondary axis or annotation for |Δ|/MCSE = 3 and rel improve = 2% thresholds.

**Cells (all five re-run):**

| n | miss | mean Δ | MCSE | |Δ|/MCSE | rel improve | G4 confirm |
|---:|---:|---:|---:|---:|---:|---|
| 300 | 10% | −0.028 | 0.003 | 8.8 | 4.8% | PASS |
| 300 | 30% | −0.018 | 0.004 | 4.5 | 2.9% | PASS |
| 300 | 50% | −0.009 | 0.004 | 2.4 | 1.3% | FAIL |
| 1000 | 10% | −0.018 | 0.001 | 13.4 | 3.3% | PASS |
| 1000 | 30% | −0.011 | 0.001 | 7.9 | 1.8% | FAIL |

**Visual encoding:** Fill or outline PASS vs FAIL; manuscript text cites PASS cells only.

**Caption fence:** "60 seeds per cell; G4 confirm thresholds per prereg §3.4."

---

### Figure 3 — Calibrated gate distribution (r_cal)

**Purpose:** Show that confirmed wins occur with gate open (GNN active), not baseline fallback.

**Data:** `script/returned_gnn_campaign/results/gnn_campaign_gates.csv` (11,340 per-latent rows). Filter to F2 @ λ=1 primary arm and/or confirm RDS-derived gates if per-rep distribution needed.

**Layout:** Violin or ridge plot of r_cal by family (F1/F2/F3) or faceted by confirmed F2 cells; secondary panel: frac_gate_open by cell from cell summaries.

**Caption note:** r_cal = 0 is valid fallback; discrete-trait columns initialise near closed; continuous F2 columns show material opening in confirmed regimes.

---

## Methods (stub)

Simulation design, G1–G8 gate audit, G4 explore/confirm rules, paired Δ estimand, seed ladders, and compute layout are specified in the Phase A pre-registration (`docs/dev-log/2026-08-26-gnn-evidence-preregistration.md`, local gitignored copy; prereg committed at evidence-lane G0). The GNN architecture (graph transformer backbone, rate-aware attention, calibrated per-trait gate, conformal intervals, baseline blend) is documented in the package vignette `vignettes/gnn-architecture.Rmd` (pkgdown: `articles/gnn-architecture.html`). Phase A drivers: `script/gnn_evidence_campaign.R`, collectors `script/collect_gnn_evidence_campaign.R`, `script/collect_gnn_evidence_f2_confirm.R`, `script/collect_gnn_evidence_bayes_sensitivity.R`. Primary results: `script/returned_gnn_campaign/PHASE_A_SUMMARY.md`.

**DGP families:** F1 linear BM (null); F2 nonlinear continuous (primary); F3 mixed-type eight-trait panel (descriptive). **Missingness:** MCAR on latent test cells. **Comparison:** blend vs baseline with `lambda_mode = "fixed_1"` (primary) or `"bayes"` (sensitivity). **Metrics:** held-out RMSE-based Δ, MCSE across seeds, relative improvement vs baseline loss.

---

## Internal cross-reference

| Artifact | Path |
|---|---|
| Phase A summary | `script/returned_gnn_campaign/PHASE_A_SUMMARY.md` |
| Confirm cell summary | `script/returned_gnn_campaign/results_confirm/gnn_confirm_cell_summary.csv` |
| Bayes closure | `script/returned_gnn_campaign/results_bayes/gnn_bayes_closure_by_family_lambda.csv` |
| Gate distribution | `script/returned_gnn_campaign/results/gnn_campaign_gates.csv` |
| Checkpoint / claim fence | `LOOP/checkpoint.md` |
