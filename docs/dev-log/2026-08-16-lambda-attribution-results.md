# λ-attribution campaign — results

Campaign complete 2026-08-16 on Totoro. **1,920 jobs, 7,680 trait-rows, 0 failures.**
~19 min wall at 100 workers (pre-run estimate was 20–25 min). Design and pre-registered
interpretation: `2026-08-16-lambda-attribution-design.md`. Runner `lambda_cell.R`,
summariser `summarise_lambda_attr.R`, artifacts `~/pigauto_regime_map/results/lambda_attr/`
and `lambda_attr_summary.{md,rds}`.

**Quantities.** All paired within rep (seed is mode-invariant, so the four `lambda_mode`
settings see identical data and masks):

- `repair` = loss_base(mode) − loss_base(fixed_1). Negative = fitting λ improves the baseline.
- `gnn_res` = loss_blend(mode) − loss_base(mode). Negative = the GNN still adds *on top of*
  a λ-fitted baseline.
- `closure` = fraction of the fixed_1 GNN gain that the λ-fitted baseline reproduces
  with no neural network at all.

MCSE = sd(paired Δ)/√30; detection at |Δ|/MCSE ≥ 3.

## The headline: the low-λ gains were the baseline, and λ-fitting takes all of them

The regime map's largest effect was Δ growing ~100× as λ fell. The design pre-registered
that ≥70% closure would re-label those rows. **Measured closure at λ = 0.2 is 1.00–1.17**
— the λ-fitted baseline reproduces *all* of it, in both families, at both n:

| family | λ_DGP | n | mode | repair | MCSE | closure | gnn_res | MCSE |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| F1 | 0.2 | 300 | bayes | −0.1402 | 0.0062 | 1.14 | +0.0002 | 0.0033 |
| F1 | 0.2 | 300 | estimate | −0.1374 | 0.0064 | 1.12 | −0.0022 | 0.0030 |
| F2 | 0.2 | 300 | bayes | −0.2519 | 0.0156 | 1.17 | −0.0084 | 0.0099 |
| F1 | 0.2 | 100 | bayes | −0.1422 | 0.0134 | 1.00 | +0.0007 | 0.0049 |

At F1 λ=0.2 the GNN's residual contribution over the λ-fitted baseline is **statistically
indistinguishable from zero** (+0.0007 ± 0.0049; −0.0022 ± 0.0030). The neural correction
was absorbing a misspecification that a one-parameter baseline fix removes outright.

**Consequence for the manuscript:** the low-λ regime-map rows must not be reported as GNN
accuracy. They are the cost of the λ=1 default, and the remedy already ships
(`lambda_mode`). The robustness reading ("degrades gracefully when the phylogenetic model
is wrong") survives but is now bounded — the cheaper remedy recovers the same ground.

## What survives: the F2 nonlinear lift, on top of a λ-fitted baseline

The pre-registered discrimination came out the same way it did in the regime map, and this
time it cannot be explained by misspecification, because the baseline now fits λ:

| family | λ_DGP | n | mode | gnn_res | MCSE | \|Δ\|/MCSE |
|---|---:|---:|---|---:|---:|---:|
| **F2** | 1.0 | 300 | bayes | −0.0169 | 0.0026 | **6.5** |
| **F2** | 0.5 | 300 | bayes | −0.0626 | 0.0074 | **8.5** |
| **F2** | 0.8 | 100 | bayes | −0.0642 | 0.0108 | **5.9** |
| F1 | 1.0 | 100 | bayes | +0.0055 | 0.0061 | 0.9 (null ✓) |
| F1 | 0.2 | 100 | bayes | +0.0007 | 0.0049 | 0.1 (null ✓) |

F2's nonlinear cross-trait structure (`sin(2Z₁)`, `Z₂²·sign(Z₂)`) is recovered by the GNN
at every λ once the baseline is λ-fitted. F1 — the linear specificity control — stays null
almost everywhere, which is the outcome the design said would validate the comparison.

**One honest exception:** F1 λ=1.0 n=300 shows gnn_res = −0.0028 ± 0.0005 (5.6 MCSE). It is
statistically detectable but ~6× smaller than F2's smallest effect and ~0.3% of baseline
loss. Report it as detectable-but-negligible rather than pretending it is zero.

## The ML boundary pathology is real, and `bayes` is the fix

NEWS (`0.10.0.9000`, "Known limitation") flagged that ML λ can collapse to the boundary on
weak-signal traits (the BIEN `sla` case, λ̂ ≈ 0.005), and deferred the remedy to "a v0.11
target via cross-validation lambda selection". **This is the first multi-seed measurement,
and it says the deferred remedy was the wrong one:**

| λ_DGP | n | P(λ̂ < 0.05): estimate | cv | bayes |
|---:|---:|---:|---:|---:|
| 0.2 | 100 | **0.47** | 0.24 | **0.00** |
| 0.2 | 300 | 0.22 | 0.06 | 0.01 |
| 0.5 | 100 | 0.20 | 0.08 | 0.00 |

`cv` halves the boundary rate but does not remove it, and it is the *worst* of the three on
accuracy at n=100 (F1 λ=0.2: closure 0.93 vs 1.00; repair −0.111 vs −0.142). `bayes`
eliminates the boundary collapse entirely and gives the best repair at low λ. **The NEWS
guidance ("keep `fixed_1` unless verified"; "v0.11 CV-λ") should be updated to recommend
`bayes` for weak-signal continuous traits.**

λ̂ recovery is downward-biased at low signal (mean λ̂ 0.13–0.25 at λ_DGP = 0.2; 0.89–0.99 at
λ_DGP = 1.0) — expected shrinkage behaviour, and the reason the boundary rate matters.

## Cross-reference: this explains the BACE gap measured the same day

`2026-08-16-external-comparison-results.md` records the first working pigauto-vs-BACE
head-to-head: BACE beats pigauto on all four AVONET continuous traits (Tarsus RMSE 11.8 vs
23.3) while losing on categorical. BACE's mechanism is MCMCglmm's posterior over the
phylo/residual variance ratio — *Bayesian averaging over what pigauto calls λ*. This
campaign shows pigauto's own `lambda_mode = "bayes"` does the same job and is off by
default. **Untested and now the obvious next experiment:** re-run the AVONET head-to-head
with `lambda_mode = "bayes"` and see how much of the continuous-trait gap closes.

## Regime, and what this does not establish

- Simulated DGPs only (BM on λ-transformed trees, K=4, Σ exchangeable r=0.7); **MCAR 30%**;
  single-obs; **continuous only**; **no covariates**; 30 reps; n ∈ {100, 300}.
- The no-covariate/continuous-only restriction is forced, not chosen: λ modes set
  `force_per_column` (`fit_baseline.R:224`) and λ is silently dropped under covariates
  (`:546–548`), so any other configuration would confound λ-fitting with the joint-vs-
  per-column and covariate-handling axes. A covariate arm needs the Slice-D fix first.
- λ never reaches binary/categorical (LP), threshold-joint, or OVR paths — F3/discrete was
  excluded for that reason. **Nothing here licenses a λ claim for discrete traits.**
- λ̂ diagnostics are recomputed on all input-observed cells; the in-fit estimate additionally
  masks val/test, so λ̂ values are approximate (the fit discards its own λ̂ at
  `fit_baseline.R:550` — worth exposing, filed as a Slice-D candidate).
- n = 1000 not run (the regime map showed the λ effect at every n; it would have tripled cost
  without changing the shape).

## What this licenses for the paper (Shinichi's decision, not this lane's)

- **Supported:** "the gated GNN recovers nonlinear cross-trait structure that a
  correctly-specified linear phylogenetic baseline cannot represent" — F2, at λ-fitted
  baselines, 5.9–8.5 MCSE.
- **Supported:** the safety-gated architecture framing (F1 null; the gate closes where the
  baseline is right).
- **NOT supported and now positively contradicted:** any reading of the low-λ regime-map
  gains as GNN accuracy.
- **New, and defensible:** pigauto's λ toolkit works, `bayes` is the mode that works, and
  the default should probably change — but that needs a real-data confirmation before it
  ships as a recommendation.
