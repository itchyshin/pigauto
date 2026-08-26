# GNN evidence campaign — pre-registration

**Status:** pre-registered **before** any full-campaign compute  
**Date:** 2026-08-26 MDT  
**Lane:** `evidence/gnn-sentinel-prerun`  
**Candidate SHA:** `6fddd79` (pigauto 0.11 trust-and-usability candidate)  
**Handover:** `docs/dev-log/handover/2026-08-26-cursor-handover.md` @ `f808767`  
**Pre-run:** `script/gnn_evidence_sentinel_prerun.R` — 12/12 Totoro fits, 0 failures (G0 complete)

This document pre-registers the factorial grid, phased execution plan, estimands,
acceptance gates, and stop rules for the GNN evidence campaign. **No number from
the 12-fit sentinel pre-run may be cited as campaign evidence** — it measured wall
time and observability only.

---

## 1. Study objective (narrow)

Estimate the **incremental held-out predictive value** of pigauto's calibrated
three-way blend over its **same-fit decoded phylogenetic baseline**, answering:

> When does the calibrated GNN weight (`r_cal_gnn`) open, when does it contribute
> measurable held-out accuracy beyond the baseline from the same fit, and when does
> it safely close to a valid fallback?

This lane does **not** claim blanket GNN superiority, external package rankings,
or that a nonzero gate implies benefit. The package's central value is unified
mixed-type imputation with calibration safety (`r_cal = 0` always valid).

---

## 2. Estimand

**Primary estimand (per cell, per replicate):**

```
Δ = loss_blend − loss_baseline
```

evaluated on **held-out test cells** via the package-convention masked surface
(`.mask_observed_idx = c(val_idx, test_idx)`). Blend and baseline come from the
**same fit** — the baseline is the decoded `baseline$mu` from that fit, not an
external comparator.

- Loss: z-scale RMSE (continuous/count/ordinal/proportion), 0–1 error
  (binary/categorical).
- **Negative Δ** = the calibrated blend beats the same-fit baseline on held-out
  cells.
- Per cell: report mean Δ, MCSE = `sd(Δ)/√n_rep`, and the full paired
  distribution.

**Gate descriptors (secondary, not evidence of benefit):**

- `open` = `r_cal_gnn > 0`
- `materially open` = `r_cal_gnn ≥ 0.10`
- `closed` = `r_cal_gnn = 0` exactly

Record all three calibration weights (`r_cal_bm`, `r_cal_gnn`, `r_cal_mean`) per
fit. Report gate distributions, not only means (historical gates are bimodal).

---

## 3. Factorial grid and phased execution

### 3.1 Design lineage

Grid structure inherits the regime-map campaign (`af2e80a`,
`docs/dev-log/2026-08-15-regime-map-design.md`, issue #135) with these updates:

- Candidate pinned at `6fddd79` (not `c655d75`).
- λ levels: **{0.2, 0.5, 1.0}** (drops λ = 0.8; retains the discrimination axis).
- Phased missingness-mechanism rollout (Phase A = MCAR only).
- Primary baseline arm: `lambda_mode = "fixed_1"`; sensitivity arm documented below.

### 3.2 Phase A — budget-conscious first pass (G0 Phase A approval)

| Factor | Levels |
|---|---|
| DGP family | F1 linear specificity control, F2 nonlinear cross-trait, F3 mixed-type |
| n species | 100, 300, 1000 |
| Pagel λ (DGP) | 0.2, 0.5, 1.0 |
| MCAR missingness | 10%, 30%, 50% |
| Missingness mechanism | **MCAR only** |
| Replicates | **30 paired seeds per cell** |

**Cells:** 3 × 3 × 3 × 3 = **81 cells**  
**Fits:** 81 × 30 = **2,430 fits**

Families (regime-map convention):

- **F1:** 4 continuous traits, exchangeable Σ (r = 0.7). Specificity control —
  expected Δ ≈ 0; a GNN win here is a red flag.
- **F2:** traits 1–2 as F1; trait 3 = `sin(2·Z₁) + 0.5·Z₃`; trait 4 =
  `Z₂²·sign(Z₂) + 0.5·Z₄`. Nonlinear cross-trait links the linear joint baseline
  cannot represent — the historically defensible niche (`5677b2a`).
- **F3:** continuous + binary + 3-class categorical + count. Exercises
  LP/threshold-joint baselines and discrete gate calibration.

Tree: `ape::rtree(n)`. Traits simulated under BM on a Pagel-λ-transformed tree.
Fits: `impute()` defaults (`safety_floor = TRUE`, 500 epochs).

### 3.3 Phase B — mechanism axis (before any MAR sentence)

Phase B **requires separate G0 approval** after Phase A gate audit. No MAR claim
may be written from Phase A alone.

Add three missingness mechanisms (same factorial skeleton as Phase A):

| Mechanism | Label rule |
|---|---|
| phylo-MAR | missingness probability depends on phylogenetic position |
| genuine covariate-MAR | missingness driven by an **observed** covariate in the DGP |
| MNAR | missingness depends on the missing value itself |

**Additional cells:** 81 × 3 mechanisms = **243 cells** (9720 fits at 30 seeds)  
**Cumulative:** 324 cells, 12,960 fits

**G6 enforcement:** call a cell MAR **only** if an observed covariate actually
drives missingness. Processes with no covariate driver are labeled MCAR, never
trait-MAR. (`script/sim_bench/mechanism_sweep_methodology.md` documents a prior
degenerate case where apparent trait-MAR collapsed to MCAR in 3/4 scenarios.)

### 3.4 F2 @ λ = 1 priority cell

Historical evidence (`5677b2a`) identified F2 nonlinear structure at λ = 1.0 as
the defensible niche — where the BM baseline is correctly specified and the GNN
still recovers structure the baseline cannot represent. Lambda-attribution
(`a64671e`) confirmed low-λ gains close under `lambda_mode = "bayes"`; the F2
@ λ = 1 lift survives on top of a λ-fitted baseline.

**Priority protocol for F2 @ λ = 1.0:**

1. **Explore:** 30 seeds across all (n, missingness) cells in this row-block
   (9 cells × 30 = 270 fits; subset of Phase A).
2. **Confirm:** any cell proposed for a public GNN-positive claim must pass G4
   on 30 seeds, then run a **separate 60-seed confirmation** before the claim
   enters manuscript or release prose.

This two-stage gate is mandatory for F2 @ λ = 1 and recommended for any other
cell crossing G4.

---

## 4. Primary vs sensitivity arms

| Arm | `lambda_mode` | Role |
|---|---|---|
| **Primary** | `"fixed_1"` | Main estimand; matches user default and historical regime-map comparator |
| **Sensitivity** | `"bayes"` | **Low-λ cells only** (λ_DGP ∈ {0.2, 0.5}); tests whether incremental GNN value survives a λ-aware baseline |

Rules:

- `lambda_mode = "bayes"` is **not** the primary estimand unless explicitly
  pre-registered and approved as primary in a future amendment.
- Low-λ Phase A rows must report **both** arms side by side.
- A low-λ GNN claim requires primary-arm G4 pass **and** sensitivity-arm
  `gnn_res` (blend − λ-fitted baseline) with `|Δ|/MCSE ≥ 3` — mirroring
  `a64671e` closure logic.
- λ = 1.0 cells run primary arm only (baseline is correctly specified).

---

## 5. Acceptance gates (G1–G8)

Verbatim from handover (`f808767`) unless noted.

### G1 — Provenance

Every result row names: candidate SHA, driver SHA, host, seed, trait family,
signal (λ), missingness rate, missingness mechanism, and outcome.

### G2 — Paired isolation

Every row records: baseline loss, blend loss, paired Δ, and all three calibration
weights (`r_cal_bm`, `r_cal_gnn`, `r_cal_mean`).

### G3 — Fallback

Zero GNN weight remains reachable and yields the valid baseline corner. Retain
calibration-floor interventions and failures. A closed gate is success, not
failure.

### G4 — Positive GNN claim

A cell may be proposed for a public GNN-benefit statement **only if all** of:

1. **Nonlinear structure** — cell is F2 (or F3 sub-component pre-specified in
   the analysis plan), not F1.
2. **Lambda-aware baseline** — for λ_DGP < 1.0, sensitivity arm with
   `lambda_mode = "bayes"` shows residual GNN value (`gnn_res`) passing the
   MCSE rule below.
3. **`|Δ|/MCSE ≥ 3`** on the primary paired estimand (30-seed explore stage).
4. **Relative improvement ≥ 2%** vs same-fit baseline loss.
5. **60-seed confirmation** — separate confirmatory run before any public claim.

### G5 — No-benefit retention

F1 null results, lambda-repair closure, material blend loss (Δ > 0), failed
fits, and floor interventions are **retained and reported by cell**. Absence of
benefit is a valid finding.

### G6 — Missingness labeling

Call a cell MAR only if an observed covariate actually drives missingness;
otherwise label MCAR. Phase A is MCAR-only by construction.

### G7 — Trait-family boundary

Report F1, F2, and F3 independently. Do not pool across families into a broad
"GNN wins" headline. Within F3, report by trait type where sample size permits.

### G8 — Stop rules

Pause and report if:

- **> 20%** of fits fail in any cell,
- any required extraction field is missing from output RDS, or
- wall time exceeds the **approved** estimate for the current phase.

Do not quietly continue past an overrun — re-estimate and seek re-approval.

---

## 6. Compute plan and budget (Totoro)

### 6.1 Pre-run timings (measured, G0 sentinel)

Host: Totoro, 12 parallel jobs, `OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`, candidate `6fddd79`.

| n | mean fit_sec | n fits measured |
|---|---:|---:|
| 100 | 38.2 | 3 |
| 300 | 49.8 | 6 |
| 1000 | 187.7 | 3 |

Pre-run wall: ~3.2 min for 12 fits, 0 failures. All G1–G3, G5–G8 gates passed
(see `LOOP/checkpoint.md`).

### 6.2 Phase A estimate @ 100 workers

| Quantity | Value |
|---|---|
| Cells | 81 |
| Fits | 2,430 |
| Total CPU-seconds (measured weights) | ~223,300 |
| Wall (optimistic, perfect parallelism) | **~37 min** |
| Wall (conservative, regime-map anchor + overhead) | **~1.0–1.5 h** |
| Approval band | **≤ 2 h** (stop at G8 if exceeded) |

Breakdown by n (810 fits each): n = 100 dominates count but not time; n = 1000
dominates wall (~63% of total CPU-seconds).

### 6.3 Phase B incremental estimate (provisional)

243 additional cells × 30 seeds = 7,290 fits. Mechanism generators add negligible
per-fit cost; expect **~3–4 h** wall @ 100 workers (provisional — revise after
Phase A measured throughput).

### 6.4 Execution conventions

- **Host:** Totoro, ≤ 100 PSOCK workers, `OPENBLAS_NUM_THREADS=1`,
  `OMP_NUM_THREADS=1`.
- **Layout:** commit-pinned driver, per-(cell, seed) RDS, resumable; master seed
  derived per cell × seed.
- **Provenance in every RDS:** candidate SHA, driver SHA, host, seed, family,
  n, λ, missingness, mechanism, `lambda_mode`, losses, Δ, `r_cal_*`, floor flag.

---

## 7. Explicit out of scope

The following are **excluded** from this campaign and must not be inferred from
its results:

- **BACE comparison** or any external package head-to-head.
- **Pooled cross-type wins** — no aggregating F1/F2/F3 into a single accuracy
  figure for public claims.
- **Nonzero gate = gain** — an open gate is a descriptive outcome, not evidence
  of benefit (G4 required for claims).
- **Pre-run numbers as claims** — the 12-fit sentinel measured timing and
  pipeline observability only.
- **MAR claims from Phase A** — MCAR only until Phase B completes.
- **Product-lane work** — `codex/pigauto-0-11-trust-usability` / PR #174 is
  PROTECTED; this lane does not modify candidate source.
- **Release, merge, push, or public posting** — evidence lane only.

---

## 8. Historical context (not re-run evidence)

| Artifact | Commit | Retained reading |
|---|---|---|
| Regime-map results | `5677b2a` | F1 null at λ = 1; F2 incremental lift in majority of cells |
| Lambda-attribution | `a64671e` | Low-λ gains close under `lambda_mode = "bayes"`; F2 @ λ = 1 survives |
| Mechanism sweep | `script/sim_bench/mechanism_sweep_methodology.md` | Continuous single-trait n = 500; cannot support general MAR |
| Ablation verdict | `useful/phase1_gnn_ablation_verdict.md` | Historical comparators unusable; paired same-fit design required |

This campaign re-measures under candidate `6fddd79` with the paired estimand and
phased grid above. Historical numbers are context, not substitutes.

---

## 9. Provenance requirements

Every output artifact must record:

```
candidate_sha, driver_sha, host, seed, family, n_species,
phylo_lambda, miss_frac, missingness_mech, lambda_mode,
blend_loss, baseline_loss, paired_delta,
r_cal_bm, r_cal_gnn, r_cal_mean, floor_fired, fit_sec
```

Campaign driver path (to be committed before launch):
`script/gnn_evidence_campaign.R` (Phase A driver; Phase B extension or sibling).

---

## 10. G0 Phase A approval request

> **Copy-paste block for Shinichi — G0 Phase A (full campaign, MCAR axis)**

```
G0 Phase A — pigauto GNN evidence campaign (MCAR axis)

Lane:      evidence/gnn-sentinel-prerun @ evidence worktree
Candidate: 6fddd79
Prereg:    docs/dev-log/2026-08-26-gnn-evidence-preregistration.md (committed BEFORE compute)

Scope:
  - 81 cells: F1/F2/F3 × n∈{100,300,1000} × λ∈{0.2,0.5,1.0} × miss∈{10%,30%,50%} × MCAR
  - 30 paired seeds per cell → 2,430 fits
  - Primary arm: lambda_mode = "fixed_1"
  - Sensitivity: lambda_mode = "bayes" on λ∈{0.2,0.5} cells only
  - F2 @ λ=1: 30-seed explore; 60-seed confirm before any public GNN claim

Pre-run:   12/12 Totoro fits, 0 failures, G1–G8 pass (sentinel complete)
Host:      Totoro, ≤100 workers, OPENBLAS=1, OMP=1

Wall estimate @100 workers:
  - Measured CPU/n: ~37 min (optimistic)
  - Conservative (overhead anchor): ~1.0–1.5 h
  - Approval ceiling: ≤2 h (G8 stop if exceeded)

Out of scope: BACE, external rankings, pooled wins, MAR claims (Phase B),
               PR #174 / product lane, push/merge/release

Approve Phase A launch?  [ ] yes  [ ] no  [ ] revise scope
```

---

*Pre-registered 2026-08-26. No full-campaign compute authorized until G0 Phase A
approval above.*
