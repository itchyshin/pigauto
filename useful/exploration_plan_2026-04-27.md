# pigauto exploration plan — 2026-04-27 PM

> Saved for the planning conversation. Goal: figure out what the
> GNN/autoencoder is actually doing in pigauto, reconcile the
> inconsistencies in past sim results, and produce a defensible
> paper claim.

## Why this plan

After a week of sims, we have:
- One designed-to-win bench where pigauto beats phylolm-λ (multi-obs
  nonlinear smoke, 14/18 cells, 6–10 % at β=1.0)
- Most other benches: pigauto matches phylolm-λ, occasionally regresses
- Architecture is irrelevant (transformer / GAT / GCN within 1.21 %)
- Discrete trait types regressed Apr 17 → Apr 27

But: pigauto IS named after the auto-encoder. Literature (MIDAS,
MIWAE, GAIN, denoising AEs) clearly supports AE-based imputation
when nonlinearity + missingness are present. We expect the GNN to
add value in low-phylo + nonlinear regimes — and we did see this on
multi_obs_nonlinear. So the AE story isn't dead; it's just not yet
isolated from the rest of pigauto's machinery.

## What I haven't done that matters

1. **Never ablated the GNN within pigauto.** Forcing the calibrated
   gate to zero (so `pred = baseline + cov_linear * X`) and comparing
   to full pigauto would tell us exactly what the GNN delta earns.
2. **Never inspected what the GNN learns.** No training loss curves,
   no calibrated-gate distributions, no `delta` magnitude
   distributions, no attention-weight inspection.
3. **Reps too few** (1–3 per cell) — discrete "regression" might be MC
   noise.
4. **Never compared to a non-phylogenetic AE imputer** (MIDAS, MIWAE,
   simple tabular DAE) — we don't know whether pigauto's phylo
   machinery is what wins, or just the AE itself.
5. **No systematic sweep of λ (phylo signal)** down to 0.
6. **No systematic sweep of n_species** — AE methods typically need
   n > 500 to show clear gains.
7. **No systematic sweep of missingness** — literature wins are at
   0.50–0.70; we mostly tested 0.25–0.30.

## Inconsistencies to reconcile

- bench_multi_obs.R (Apr 16) with ONE covariate: pigauto wins
  10–19 %.
- bench_multi_obs_nonlinear (Apr 27) with FIVE covariates and
  nonlinear function: pigauto wins 6–10 %.
- More covariates + nonlinearity should give the AE *more* to learn.
  Why is the win smaller? Need to understand.
- Per-type benches (single-obs): pigauto matches baseline on every
  scenario including nonlinear. Why does the AE not earn anything
  there?

## Five-phase plan

### Phase 1 — Inspect what the GNN does (1–2 days)
- Force `r_cal = 0` on bench_multi_obs_nonlinear smoke cells where
  pigauto currently wins. Does pigauto still win?
- Plot calibrated gate values per trait per regime
- Plot training loss curves
- Plot magnitude of `delta` predictions vs `(truth - baseline)`
  residuals
- **Decisive question:** is the GNN delta what wins, or is it
  the multi-obs aggregation / GLS baseline / cov_linear path?

### Phase 2 — Map pigauto's home turf (2–3 days)
A single systematic sweep on the multi-obs nonlinear DGP:
- n_species ∈ {100, 300, 500, 1000}
- phylo_signal ∈ {0, 0.05, 0.1, 0.3, 0.6}
- n_covs ∈ {2, 5, 10}
- missingness ∈ {0.25, 0.50}
- f_type ∈ {linear, nonlinear, interactive}
- 5 reps for statistical power
Goal: a heatmap showing where pigauto beats phylolm-λ by ≥ X %.
This is the paper's headline figure.

### Phase 3 — Connect to literature (1 day)
Add a non-phylogenetic AE imputer (MIDAS-style tabular DAE, or a
simple PyTorch AE with the same DAE training objective). Compare
on the same cells. Question: does pigauto's phylo machinery add
value over a generic AE?

### Phase 4 — Reconcile inconsistencies (1–2 days)
- Investigate why bench_multi_obs.R (1 cov) outwins
  bench_multi_obs_nonlinear (5 covs).
- Confirm the discrete regression with ≥ 5 reps; if real, git bisect
  Apr 17 → Apr 27 against bench_binary.

### Phase 5 — Conservative paper (1 week)
Position pigauto as: **"a phylogenetic graph-based denoising
autoencoder for trait imputation, with calibrated safety gating,
that adds value in regimes characterised by [whatever Phase 2
identifies]"**.
- The AE/GNN is in the headline (consistent with the package name)
- Architecture is pragmatic (legacy attention is fine)
- Multi-obs aggregation, safety gate, mixed-type API, conformal UQ
  are infrastructure contributions

## Order

Start Phase 1 (cheap, decisive). Stop and think.

If Phase 1 confirms the GNN is doing the work → proceed to Phase 2.
If Phase 1 shows the GNN delta is near zero even on winning cells →
diagnose why, and either fix calibration or accept that pigauto's
value is multi-obs aggregation, not the AE.

Either way: stop after Phase 1, share findings, decide next step.
