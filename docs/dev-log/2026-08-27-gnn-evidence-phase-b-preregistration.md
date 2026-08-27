# GNN evidence campaign — Phase B pre-registration addendum

**Status:** pre-registered **before** Phase B compute  
**Date:** 2026-08-27 MDT  
**Lane:** `evidence/gnn-sentinel-prerun`  
**Candidate SHA:** `6fddd79` (unchanged from Phase A)  
**Base prereg:** `docs/dev-log/2026-08-26-gnn-evidence-preregistration.md`  
**G0 approval:** Shinichi approved tasks 1+2+3 (= Phase B) 2026-08-27

This addendum registers the missingness-mechanism axis for Phase B. **No MAR
sentence may be written from Phase A (MCAR-only) results.** Phase B must
complete before any mechanism-axis claim.

---

## 1. Scope

Phase B adds three missingness mechanisms to the **same factorial skeleton**
as Phase A primary arm:

| Factor | Levels |
|---|---|
| DGP family | F1, F2, F3 |
| n species | 100, 300, 1000 |
| Pagel λ (DGP) | 0.2, 0.5, 1.0 |
| Missingness rate | 10%, 30%, 50% |
| **Missingness mechanism** | **phylo_MAR, covariate_MAR, MNAR** |
| Replicates | 30 paired seeds per cell |
| `lambda_mode` | `"fixed_1"` (primary arm only) |

Phase A MCAR results are **not re-run** in Phase B.

---

## 2. Mechanism definitions and G6 labeling

| Mechanism code | Generator | G6 label |
|---|---|---|
| `phylo_MAR` | `P(miss) ∝ plogis(c − 1.5·scale(z))`, `z ~ MVN(0, R)` phylogenetic latent; calibrated intercept `c` targets `miss_frac` | **Not MAR** — no observed covariate driver |
| `covariate_MAR` | Driver = first continuous trait (`trait1` / `cont1`), **fully observed**; remaining traits punched with `P(miss) ∝ plogis(1.5·scale(driver))`, calibrated to `miss_frac` | **MAR** — observed covariate drives missingness |
| `MNAR` | Per-cell `P(miss) ∝ plogis(1.5·scale(y))` where `y` is the value being hidden; calibrated to `miss_frac` | **Not MAR** — missingness depends on hidden value |

**G6 enforcement:** only cells with `missing_mechanism == "covariate_MAR"` may be
called MAR in prose. `phylo_MAR` is phylogenetically structured missingness,
not covariate-MAR. `MNAR` is value-dependent missingness.

DEP_STRENGTH = 1.5 matches the mechanism sweep
(`script/sim_bench/mechanism_sweep_methodology.md`, Penone et al. 2014 ICCs).

---

## 3. Cell and fit counts

| Quantity | Value |
|---|---|
| Base cells (Phase A skeleton) | 81 |
| Mechanisms | 3 |
| **Phase B cells (full grid)** | **243** |
| Seeds per cell | 30 |
| **Phase B fits (full grid)** | **7,290** |
| Cumulative (A + B primary) | 324 cells, 9,720 fits |

### Parallel compute split (lanes 3a + 3b)

Phase B launches as **two Totoro jobs** in parallel (same G0, disjoint mechanisms):

| Lane | Mechanisms | Cells | Fits | Totoro dir | Launcher |
|---|---|---:|---:|---|---|
| **3a** | phylo_MAR + covariate_MAR | 162 | 4,860 | `~/pigauto_gnn_evidence_phase_b_mar` | `script/gnn_evidence_phase_b_mar_totoro.sh` |
| **3b** | MNAR only | 81 | 2,430 | `~/pigauto_gnn_evidence_phase_b_mnar` | `script/gnn_evidence_phase_b_mnar_totoro.sh` |

Global `cell_id` is preserved across arms (phylo_MAR 0–80, covariate_MAR 81–161,
MNAR 162–242). Per-arm `job_id` is 0-based within the arm.

No within-arm subset — full factorial on each mechanism axis.

---

## 4. Wall-time estimate (Phase A measured throughput)

Phase A primary (2,430 fits, 100 Totoro workers): **4,736 s (~79 min)**, 0 failures.

| Quantity | Value |
|---|---|
| Phase B fits | 7,290 (= 3× Phase A) |
| Linear extrapolation (monolithic) | 4,736 × 3 = **14,208 s (~237 min, ~3.95 h)** |
| **Lane 3a (MAR arms)** | 4,860 fits → **~158 min** linear |
| **Lane 3b (MNAR arm)** | 2,430 fits → **~79 min** linear (Phase A anchor) |
| Parallel wall (max of 3a, 3b) | **~158 min (~2.6 h)** @ 100 workers each |
| Conservative (overhead anchor) | **~2.5–3.0 h** parallel |
| **G8 approval ceiling (monolithic)** | **≤ 5 h** |
| **G8 approval ceiling (lane 3b MNAR)** | **≤ 1.5 h** |

Mechanism generators add negligible per-fit cost (same `impute()` path as Phase A).

Breakdown by n (2,430 fits each, same mix as Phase A):

| n | Phase B fits | Est. CPU share |
|---|---:|---:|
| 100 | 2,430 | ~15% |
| 300 | 2,430 | ~22% |
| 1000 | 2,430 | ~63% |

---

## 5. Driver and launcher paths

| Artifact | Path |
|---|---|
| Phase B driver | `script/gnn_evidence_campaign_phase_b.R` |
| Totoro launcher (monolithic) | `script/gnn_evidence_campaign_phase_b_totoro.sh` |
| Totoro launcher (lane 3a MAR) | `script/gnn_evidence_phase_b_mar_totoro.sh` |
| Totoro launcher (lane 3b MNAR) | `script/gnn_evidence_phase_b_mnar_totoro.sh` |
| Rsync push/pull (monolithic) | `script/rsync_gnn_evidence_phase_b.sh` |
| Rsync push/pull (lane 3b) | `script/rsync_gnn_evidence_phase_b_mnar.sh` |
| Collector (monolithic) | `script/collect_gnn_evidence_phase_b.R` |
| Collector (lane 3b) | `script/collect_gnn_evidence_phase_b_mnar.R` |
| Totoro results dir (monolithic) | `~/pigauto_gnn_evidence_phase_b/results_phase_b/` |
| Totoro results dir (lane 3b) | `~/pigauto_gnn_evidence_phase_b_mnar/results_phase_b_mnar/` |

Driver env hooks: `PIGAUTO_MECHANISM_ARM` (`all` | `MAR` | `MNAR`),
`PIGAUTO_JOB_PREFIX` (default `gnn_phase_b`; lane 3b uses `gnn_phase_b_mnar`).

Output RDS (lane 3b): `gnn_phase_b_mnar_job_<id>.rds` (job_id 0..2429;
global cell_id 162..242).

Seed block: `30000 + cell_id × 100 + rep` (disjoint from Phase A `10000`-block).

---

## 6. Acceptance gates (unchanged from base prereg)

G1–G8 apply verbatim. Phase B-specific checks:

- **G6:** collector asserts `covariate_MAR` is the only mechanism labeled MAR;
  `phylo_MAR` and `MNAR` never appear as MAR in summaries.
- **G8:** pause if wall > 5 h or > 20% cell failure rate.

**Out of scope (unchanged):** BACE, external rankings, PR #174, push/merge/release,
Bayes sensitivity re-run on Phase B (primary arm only).

---

## 7. Analysis plan (post-collection)

Descriptive by mechanism × family × (n, λ, miss):

1. Mean paired Δ, MCSE, gate-open fraction — same estimand as Phase A.
2. Compare Phase B cells to Phase A MCAR matched cells (same family/n/λ/miss).
3. **No pooled "GNN wins under MAR" headline** — report by mechanism with G6 labels.
4. F2 @ λ = 1 explore-only screen (no new 60-seed confirm unless pre-approved).

---

## 8. G0 Phase B approval record

```
G0 Phase B — pigauto GNN evidence campaign (mechanism axis)

Lane:      evidence/gnn-sentinel-prerun @ evidence worktree
Candidate: 6fddd79
Prereg:    docs/dev-log/2026-08-27-gnn-evidence-phase-b-preregistration.md
           (addendum to 2026-08-26 base prereg)

Scope:
  - 243 cells: F1/F2/F3 × n × λ × miss × {phylo_MAR, covariate_MAR, MNAR}
  - 30 paired seeds per cell → 7,290 fits
  - Primary arm: lambda_mode = fixed_1
  - G6: MAR label only for covariate_MAR

Phase A anchor: 2430 fits, 0 failures, wall 79 min @ 100 workers

Wall estimate @100 workers:
  - Linear: ~4.0 h
  - Conservative: ~4.0–5.0 h
  - G8 ceiling: ≤ 5 h

Out of scope: MCAR re-run, Bayes arm, PR #174, push/merge/release

Approved: Shinichi 2026-08-27 (tasks 1+2+3)
```

---

*Pre-registered 2026-08-27. Phase B compute authorized by G0 approval above.*
