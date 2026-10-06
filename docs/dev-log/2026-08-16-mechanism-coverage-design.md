# Mechanism-coverage campaign (B1) — design (ADEMP-lite)

2026-08-16 · Claude lane · Morris et al. 2019 / Williams et al. 2024 pattern.
Runner: `mech_cell.R` on Totoro (`~/pigauto_regime_map/`).

## Why this is NOT a re-run of `bench_missingness_mechanism.R`

The plan originally said "upgrade the April bench". Reading `R/mask_missing.R:85–186`
killed that: in `make_missing_splits()` the mechanism weights ONE sample of cells which
is then split into val and test — so **val and test are exchangeable with each other
even under MNAR**, and that design cannot reproduce the real-data seam. The seam the
fishbase/pantheria undercoverage points at is different: genuine NAs follow the
mechanism, while pigauto's conformal calibration cells are drawn from the **observed
complement** (`fit_pigauto` masks val cells out of observed data). This campaign
reproduces the real usage exactly:

1. simulate complete truth → 2. apply the mechanism as **genuine NAs** →
3. `impute(df, tree)` exactly as a user would (pigauto internally draws val/test from
   observed cells) → 4. score conformal coverage on the genuinely-missing cells.

The April run stays useful as a point-accuracy bench; its SE-based coverage drop
(MCAR 0.949 → MNAR 0.923) is directional but measures a different interval on a
different surface.

## A — Aims

**Primary:** confirm or kill the exchangeability hypothesis for the large-n real-data
undercoverage (fishbase 0.89–0.91 at n=10,654 where the order-statistic ceiling is
~0.998; see `2026-08-16-coverage-results.md` CORRECTION). **Decision rule
(pre-registered): the hypothesis is CONFIRMED iff conformal coverage under MAR/MNAR is
below the MCAR arm by > 3× MCSE, with the MCAR arm itself ≥ ~0.94.** If MCAR also
undercovers at these n, the defect is not (only) exchangeability — stop B2 and file a
new diagnosis lane.

**Secondary:** which mechanism hurts most (value-MNAR vs trait-MAR vs clade-MAR) —
this shapes what a weighted-conformal fix must condition on; realised n_val recorded
per rep.

## D — Data-generating mechanism

DGP: F1 continuous (4 correlated BM traits, Σ r=0.7), **λ=1 tree** — the baseline is
correctly specified, so the mechanism is the only axis varied. Overall missingness
targeted at 30% per column (probabilities rescaled; ≥20 observed per column guard).

Mechanisms (following the April bench's parametrisation):
- **MCAR** (control)
- **MAR_trait**: P(missing in t2–t4) ∝ plogis(2·z(t1)); t1 MCAR (driver observed)
- **MAR_phylo**: two clades (15–35% of tips each) at 7:1 odds vs background
- **MNAR**: P(missing tj) ∝ plogis(2·z(tj)) — depends on the trait's own value

Grid: 4 mechanisms × n {300 (30 reps), 1000 (20 reps)} = **200 jobs**.

## E — Estimands

Per job: conformal coverage of `[conformal_lower, conformal_upper]` at nominal 0.95 on
genuinely-missing continuous cells (pooled + per column); median interval width;
n_val cells; realised missing fraction. Per cell: mean coverage across reps with
MCSE = sd(per-rep coverage)/√n_rep (per-rep coverage averages hundreds of cells, so
across-rep sd captures fit-level variation; expected MCSE ≈ 0.005 at 30 reps).

## M — Methods

pigauto production path only: `impute(df, tree, epochs = 500, seed)`. No competitor
arms — this is a defect-diagnosis campaign, not a bench.

## P — Interpretation, pre-registered

- MCAR ≥ 0.94 and MAR/MNAR ≥ 0.94 → exchangeability hypothesis KILLED at these
  regimes; the real-data undercoverage needs a different explanation (candidates:
  heteroscedastic width vs the global scalar; within-species non-independence in
  multi-obs data; trait distributions far from the sim's). B2 does NOT proceed.
- MCAR ≥ 0.94, any MAR/MNAR arm < MCAR − 3×MCSE → CONFIRMED; B2 (weighted/Mondrian
  conformal) proceeds, conditioning on what the worst mechanism implies.
- MCAR itself < 0.94 at λ=1, n≥300, ~75 val cells/trait → a pipeline defect upstream
  of the mechanism question; stop and diagnose before any fix.

## D-139 gate

Analogue: regime-map n=300 jobs ≈ 100–150 s; n=1000 jobs are the heavy tail (~10–20
min each; 80 of them). **Estimate: ~30–45 min wall at 100 workers**, dominated by the
n=1000 arm. Pre-run: 3 jobs (MCAR n300, MNAR n300, MAR_phylo n1000) timed before
launch; results appended below. Overrun → stop and re-report.

## Pre-run results (2026-08-16, D-139)

3/3 OK, 0 failures. MCAR n=300: 60.4 s (cover 0.896, 1 rep — noise, not a result);
MNAR n=300: 60.5 s (0.942); MAR_phylo n=1000: 237.1 s (0.956). Revised estimate:
200 jobs ≈ 7.3 CPU-h → **~8–12 min wall at 40 workers** — under the D-139 30-minute
just-run line. Launched at 40 workers alongside lambda_attr at 100 (140 ≤ 150, D-143).
