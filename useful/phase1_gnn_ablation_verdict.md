# Phase 1 GNN ablation — verdict (with caveat)

> 2026-04-27 ~14:10. Source: bench_multi_obs_nonlinear_smoke.rds
> (already on disk) + an ablation script that did NOT complete cleanly.

## What I tried

The plan: force `r_cal_gnn = 0` post-fit on the multi-obs nonlinear
cells where pigauto wins, and see if pigauto still wins.

**The override approach failed.** Calling `predict()` on a fit object
with manually-overridden gate values crashed with a tensor shape
mismatch (`linear(): input and weight.T shapes cannot be multiplied
(1473x22 and 23x64)`). All 24 cells errored. `predict()` rebuilds the
data tensor differently than `impute()`'s internal predict call —
something in the data pipeline isn't symmetric. Worth fixing eventually
(it's a real bug for users wanting custom-gate prediction) but not the
critical path here.

## What I learned by re-examining the existing bench data

While debugging, I re-read the bench_multi_obs_nonlinear_smoke results
with `lm` and `lm_nonlinear` in the table — the **obs-level**
baselines we already had — and found a much more important issue.

### The baseline I was comparing against was a straw man

In `bench_multi_obs_nonlinear`, `phylolm_lambda_blup` aggregates to
**species means** before fitting. That throws away all within-species
covariate variation. `lm` and `lm_nonlinear` operate at the
**observation level** with full data. On multi-obs cells with strong
covariate signal, the obs-level methods have access to information
phylolm-λ doesn't.

### When I rank all six methods per cell, pigauto's "win" splits

#### Where pigauto wins (10/14 informative cells, mostly moderate phylo)

`pig < phyL < {lm, lmNL}` — pigauto best, phylolm second, OLS-based
methods worst.
- All β=0 cells (no covariate signal — phylo info is everything)
- Most β=0.5 cells across all f_types
- β=1.0 cells with **higher phylo signal (λ=0.3)**, all f_types

#### Where pigauto LOSES (4/14, low phylo + strong covariate)

`{lm, lmNL} < pig < phyL` — obs-level regression wins, pigauto and
phylolm-λ both lose.

| f_type | λ | β | lm | lm_NL | phylolm-λ | pigauto |
|---|---|---|---|---|---|---|
| linear | 0.1 | 1.0 | **0.987** | 1.003 | 1.162 | 1.037 |
| nonlinear | 0.1 | 1.0 | 0.925 | **0.892** | 1.201 | 1.080 |
| interactive | 0.1 | 1.0 | 1.264 | **0.856** | 1.182 | 1.071 |
| interactive | 0.3 | 1.0 | 1.414 | **1.052** | 1.273 | 1.192 |

These are exactly the cells where I claimed "pigauto's GNN/AE earns
its keep on low-phylo + nonlinear DGPs". But **a no-phylo polynomial
+ interaction OLS beats pigauto on all of them**, by 8–17 %. That
is the regime where literature predicts AE should win, and pigauto's
AE does NOT.

## What this actually means

1. **Pigauto's "designed-to-win" bench was a comparison against the
   wrong control.** The advantage was largely "obs-level model beats
   species-aggregated model". Pigauto IS obs-level, so it benefits;
   so does plain OLS. We're not seeing GNN value; we're seeing
   "anything that uses obs-level data beats anything that doesn't".
2. **The AE/GNN appears to add NOTHING measurable on the cells where
   it should help most.** Low phylo + strong nonlinear covariate
   signal: pigauto is 8–17 % WORSE than `lm_nonlinear`. If the AE
   were earning its keep we'd expect pigauto to be at least
   competitive there.
3. **Where pigauto wins is exactly where pigauto's analytical phylo
   machinery wins** — moderate-to-high phylo signal where the BM/GLS
   baseline genuinely captures the response structure, and OLS
   fails. This isn't an AE story; it's an analytical-baseline story.

## Reconciling with the literature

MIDAS, MIWAE, GAIN show AE wins on tabular data with **larger n
(thousands), no good phylo prior, complex feature relationships**.
Our setting (n=300, multi-obs ~1500 rows, 5 covs) might just be too
small for the AE to find structure that polynomial OLS doesn't.

## Recommended next step

Three options for the planning conversation:

### Option 1 — Fix the bench, then re-test (1–2 days)
- Add `lm_nonlinear` to ALL benches as a baseline (currently only in
  the multi_obs_nonlinear smoke).
- Re-run with reps ≥ 5 for statistical power.
- This will show whether pigauto has ANY regime where it beats
  obs-level OLS-based methods. If yes, that's the paper's regime.
  If no, the AE story is dead.

### Option 2 — Fix the predict() shape bug + re-do real Phase 1 (1 day)
- The intended ablation (force r_cal_gnn = 0 post-fit, re-predict)
  would still be informative even though Option 1's finding is more
  immediate.
- Currently blocked by a real bug in pigauto's predict path.

### Option 3 — Increase n + nonlinearity (3–4 days)
- Test on n=1000, 5000 species with truly complex DGPs (e.g., XOR-like
  feature interactions, deep nonlinearity).
- Possibly compare against MIDAS or a simple PyTorch tabular AE.
- This gives the AE an actual chance to find structure OLS can't.
- Most aligned with "doing justice to the AE story".

### My recommendation

**Option 1 first.** It's cheap, uses data we already have plus a few
re-runs with more reps. The decisive test is: **on the multi_obs_nonlinear
DGP where we thought pigauto won, does pigauto beat lm_nonlinear?**
The current data says no. With more reps we can confirm whether
that's robust.

If pigauto ≥ lm_nonlinear on enough cells with reps, the AE story is
alive (just narrower than thought). If pigauto < lm_nonlinear
across the board, the AE doesn't help on this scale and we either:
- Pursue Option 3 (bigger / harder DGPs) before claiming AE value, or
- Accept that pigauto's value is the multi-obs aggregation + safety
  gating + UQ infrastructure, not the AE.

## What I am NOT recommending right now

- Bisecting the discrete regression — still important but lower
  priority than confirming the AE story.
- Architecture variations — the architecture story is dead and not
  worth more time.
- Writing the paper — premature.

## Stop and think

This is a more pessimistic finding than yesterday's "pigauto wins
14/18". I'd like to discuss before launching anything else.
