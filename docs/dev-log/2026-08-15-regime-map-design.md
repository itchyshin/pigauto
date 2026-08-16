# Regime-map campaign design (issue #135) — ADEMP

Design per **Morris, White & Crowther (2019, Stat Med 38:2074–2102; ADEMP)** with the
**Williams et al. (2024, MEE 15:1926–1939)** 11 transparent-reporting items (self-audit at
the end). Status: **design only — D-139 gates execution** (pre-run test on Totoro, results
shown, Shinichi's approval, then the full run). Drafted 2026-08-15 while the Tier-1 bench
re-runs execute; numbers marked *(pre-run)* are finalized by the pre-run test.

Pipeline under evaluation: `main` @ `c655d75` — post-#158 (leak-free training, two-layer
gate floor) and post-#159. **No number from this campaign is comparable to any pre-#158
measurement.**

## A — Aims

- **Primary:** estimate, per regime cell, the paired accuracy difference between pigauto's
  calibrated blend and its own pure phylogenetic baseline — answering the reviewer question
  "what does the GNN buy, and where?" and fixing the paper's framing (regime-mapped lift
  vs safety-gated architecture).
- **Secondary S1:** characterize gate behaviour per regime — P(gate opens), the calibrated
  `r_cal_gnn` distribution, and how often the #157 floors intervene.
- **Secondary S2:** validate the safety claim at scale: the frequency and size of cells
  where the blend is *worse* than baseline on the delivered surface (the tail the floor
  exists to bound).
- **Secondary S3 (opportunistic, free from the same fits):** empirical conformal coverage
  per cell. The dedicated coverage-vs-n campaign remains separate; S3 is recorded, not
  designed for.

## D — Data-generating mechanisms

Tree: `ape::rtree(n)` (matches the per-type bench convention for continuity). Phylogenetic
signal imposed by simulating traits under BM on a **Pagel-λ transformed** tree.

Latent construction (K = 4 latent traits per dataset):

1. `Z_k ~ BM(σ² = 1)` on the λ-transformed tree, k = 1..4, with cross-trait covariance Σ.
2. Family-specific map from latents to observed traits (below).
3. Missingness: MCAR at rate m over the observed trait matrix (row-preserving for any
   composition). MAR/MNAR are **out of scope** — `bench_missingness_mechanism.R` owns that
   axis.

**Factors:**

| Factor | Levels | Count |
|---|---|---|
| Phylo signal λ | 0.2, 0.5, 0.8, 1.0 | 4 |
| n species | 100, 300, 1000 | 3 |
| MCAR missingness m | 0.10, 0.30, 0.50 | 3 |
| DGP family | F1 linear, F2 nonlinear, F3 mixed-type | 3 |

**Families** (the axis that decides whether the GNN *can* win — under F1 the joint
baseline is near-optimal by construction and the honest expectation is a closed gate):

- **F1 linear-Gaussian:** 4 continuous traits, Σ exchangeable with r = 0.7. The baseline's
  home turf; expected Δ ≈ 0. This family is the *specificity control* — a GNN "win" here
  would itself be suspicious.
- **F2 nonlinear:** traits 1–2 as F1; trait 3 = `sin(2·Z_1) + 0.5·Z_3`; trait 4 =
  `Z_2² · sign(Z_2) + 0.5·Z_4` (nonlinear cross-trait links the linear joint baseline
  cannot represent; the GNN's claimed niche).
- **F3 mixed-type:** continuous (Z_1), binary (threshold Z_2 at 0), 3-class categorical
  (tercile cut of `sin(2·Z_1) + 0.5·Z_3` — nonlinear class boundary), count
  (`rpois(exp(0.5·Z_4 + 0.3·Z_1))`). Exercises LP/threshold-joint baselines + discrete
  gate calibration (0-1 loss path, the April-2026 regression class).

**Cells:** 4 × 3 × 3 × 3 = **108**.

**Replicates** *(pre-run finalizes)*: paired design — blend and baseline come from the
*same fit*, so the headline MCSE is `sd(Δ)/√n_rep` on the paired per-rep difference.
From tonight's bench re-runs, rep-to-rep sd of per-rep z-RMSE ≈ 0.03–0.05; to detect
|Δ| = 0.02 with `3·MCSE < 0.02` needs `MCSE ≤ 0.0067` → `n_rep ≥ (0.045/0.0067)² ≈ 45`.

- n ∈ {100, 300}: **n_rep = 60** (margin over 45).
- n = 1000: **n_rep = 30** (fit cost dominates; MCSE ~1.4× wider, flagged in reporting).

Fits: 72 cells × 60 + 36 cells × 30 = **5,400 fits**.

## E — Estimands

Per cell (all evaluated on **held-out test cells** via the package-convention masked
surface, `.mask_observed_idx = c(val_idx, test_idx)` — the #157 surface decision; truths
stored per replicate):

- **Δ_acc** (primary): per-rep paired difference in test loss, blend − baseline. Loss =
  z-scale RMSE (continuous/count), 0-1 error (binary/categorical). Negative = GNN helps.
- **P(open)**: fraction of trait-columns with calibrated `r_cal_gnn > 0`; plus the
  `r_cal_gnn` distribution (mean, P(≥ 0.1)).
- **floor_rate**: fraction of fits where either #157 floor overrode a gate (captured from
  the `fit` object's gate values vs a re-derived pre-floor candidate, or by message
  capture — implementation detail for the runner; must be programmatic, not log-scraped
  *(pre-run resolves which)*).
- **safety tail**: P(Δ_acc > 0.05 · baseline loss) and max relative excess — the S2
  validation.
- **cover95** (S3): empirical conformal coverage on test cells, continuous-family traits.

## M — Methods

| Method | Why included |
|---|---|
| pigauto blend, `impute()` defaults (safety_floor = TRUE, 500 epochs — bench convention) | the method under evaluation, as users run it |
| pure baseline (same fit's `baseline$mu`, decoded) | the paired comparator; isolates the GNN's marginal contribution exactly |
| grand mean/mode | orients both against the no-information floor |

Excluded, with reasons: Rphylopars / missForest head-to-heads (Tier-3 comparison lane —
this campaign isolates the *internal* GNN contribution, not package rankings); pre-#158
pigauto (leak-tainted; no longer installable state); tuned-hyperparameter variants (the
paper claims default behaviour).

## P — Performance measures

Per cell × method, each with MCSE (Williams item 11):

- mean Δ_acc, MCSE = `sd(Δ)/√n_rep` — headline, `3·MCSE` decision rule against |Δ|=0.02.
- P(open), floor_rate, safety-tail P: binomial MCSE `√(p(1−p)/n_rep)` (≈ 0.065 worst-case
  at n_rep=60 — adequate for the map's qualitative contrast; flagged).
- cover95: binomial MCSE at n_rep=60 ≈ 0.028 — S3 is *descriptive only* at this size.
- convergence/failure rate per cell — failed fits reported, never dropped (item 10b).
- wall time per fit (feeds future D-139 estimates).

## Compute plan (Totoro; D-139 + D-143)

- **Host:** Totoro (384 cores, no queue). **100 PSOCK workers** (≤150 per D-143),
  `OPENBLAS_NUM_THREADS=1`, torch CPU threads = 1 per worker. Master seed 20260815;
  per-(cell, rep) seeds derived once via `sample.int`.
- **Bootstrap** *(pre-run step 0)*: install pigauto from the `c655d75` tarball into the
  Totoro R 4.5.3 library (torch already present); smoke: one `impute()` on n=50.
- **Layout:** SAFE three-stage (`0_prepare / 1_run / 2_summarise`), per-(cell, rep) RDS,
  resumable; `sessionInfo()` saved beside results. Outputs under `~/pigauto_regime_map/`
  on Totoro, rsynced back.
- **Cost estimate** *(provisional — the pre-run test replaces this)*: local MPS fits at
  n=300/500 epochs ≈ 1.6 min. Totoro CPU per-fit unknown → the point of the pre-run.
  Guess-band: n∈{100,300} 2–6 CPU-min, n=1000 15–40 CPU-min → total ≈ 700–2,000 CPU-h →
  **7–20 h wall at 100 workers**. If the pre-run lands at the top of the band, first
  action is reducing n=1000 to n_rep = 20, not exceeding a 24 h wall budget.

**Pre-run test (D-139, must be shown before approval):** 3 cells (one per n, F2, λ=0.5,
m=0.3) × 2 reps on Totoro. Pass = non-empty, non-NA RDS per rep; one fit inspected
(`str`); per-fit wall measured at each n; floor_rate capture mechanism demonstrated
working. Then the finalized rep counts + total estimate come back here for the go/no-go.

**Stop rules:** first-cell output read early — empty/NA aborts the campaign; any cell
with failure rate > 20% pauses the grid and reports; wall overrunning the approved
estimate stops and re-reports (D-139), it does not quietly continue.

## Reporting (Williams items 9–11)

Deliverables: regime-map heatmaps (Δ_acc by λ × n, faceted family × m) + the full
cell-level table with MCSEs in the supplement + gate-behaviour panel (P(open), floor_rate)
— the #135 "gate-firing diagnostics". Worked real-data anchor (item 9): the existing
AVONET300 validation (`script/validate_avonet_full.*`), re-run post-#158 in the Tier-3
lane. Code + per-cell results archived in-repo (`script/regime_map/`), commit-pinned.

## Williams et al. (2024) 11-item self-audit

| # | Item | Where |
|---|---|---|
| 1 | Aims | §A |
| 2 | DGP incl. factors/levels | §D |
| 3 | Estimands, true values stored per rep | §E |
| 4 | Methods + inclusion rationale | §M |
| 5 | Performance measures w/ formulas | §P |
| 6 | Software/versions (`sessionInfo()` archived) | §Compute |
| 7 | Seeds/reproducibility (master seed, derived per cell×rep) | §Compute |
| 8 | Code availability (in-repo, commit-pinned) | §Reporting |
| 9 | Worked real-data example | §Reporting (AVONET300) |
| 10 | Full results incl. convergence failures | §P, §Reporting |
| 11 | MCSE for every aggregate | §P (formulas + targets) |

**Next actions when compute frees:** (1) Totoro bootstrap + the 6-fit pre-run test;
(2) results back here; (3) Shinichi's go/no-go on the finalized estimate; (4) launch.
