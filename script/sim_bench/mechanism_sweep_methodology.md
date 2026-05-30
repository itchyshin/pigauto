# Mechanism-axis robustness sweep — methodology and results

**Date:** 2026-05-19 / 2026-05-20
**Drivers:** `script/sim_bench/overnight_2026_05_19_mechanisms.R`,
`script/sim_bench/pigauto_full_sweep.R`
**Data:** `dev/simulation_results_overnight_2026_05_19_mechanisms/results.rds`
(240 reps), `dev/simulation_results_pigauto/results.rds` (600 reps)
**Auto-generated companion:** `useful/MEMO_2026-05-20_mechanism_sweep.md`
**Coverage addendum:** `useful/MEMO_2026-05-21_coverage_remeasurement.md`

---

## 1. Aim

A single, falsifiable question:

> Under realistic non-random missingness, does pigauto's gated GNN
> match or beat the plain Brownian-motion kriging baseline — and does
> it ever do materially *worse*?

This is a safety-floor test, not a capability ceiling test. pigauto's
calibrated gate is designed so that `r_cal = 0` (GNN off, pure BM) is
always a reachable fallback. The claim under test is that the gate
delivers this in practice across evolutionary processes that violate
the BM assumption and across missingness mechanisms that bias the
observed sample.

The estimand is the per-cell mean z-scored RMSE on imputed cells, and
its difference between methods. The comparator of record is
`bm_kriging` — pigauto's own internal BM baseline run standalone — so
the comparison isolates the GNN's contribution from the baseline it
sits on top of.

## 2. Data-generating process

One continuous trait `y` on an `n`-tip coalescent tree
(`ape::rcoal(n)`). `R = cov2cor(vcv.phylo(tree))` is the phylogenetic
correlation matrix; `bm_signal = chol(R)' %*% rnorm(n)` is a unit-variance
Brownian draw. Four process scenarios:

| scenario | `y` = | BM assumption |
|---|---|---|
| `bm_strong` | `sqrt(0.9)·bm_signal + sqrt(0.1)·noise` | correct |
| `ou_strong` | `sqrt(0.9)·ou_signal + sqrt(0.1)·noise`, `ou_signal` from a Pagel-λ=0.3 tree | wrong (process is OU-like) |
| `nonlinear_cov` | `sqrt(0.4)·bm_signal + sqrt(0.4)·sin(2·env) + sqrt(0.2)·noise` | partially wrong; `env ~ N(0,1)` supplied as a covariate |
| `weak_signal` | `sqrt(0.1)·bm_signal + sqrt(0.9)·noise` | technically correct but near-vacuous (10% signal) |

`bm_strong` is the scenario where BM kriging is the (near-)optimal
estimator. The other three each break a different assumption: wrong
evolutionary process (`ou_strong`), an exogenous nonlinear covariate
effect BM cannot represent (`nonlinear_cov`), and a signal-to-noise
ratio so low that any phylogenetic predictor is fragile
(`weak_signal`).

## 3. Missingness conditions

30% of `y` is set missing. The missingness indicator follows a
logistic model `P(miss) = plogis(c + linpred)` with `c` calibrated by
1-D optimisation so the expected missing rate is exactly 0.30. Three
conditions for `linpred` (`DEP_STRENGTH = 1.5`, matching the
intra-class correlations in Penone et al. 2014 and Dan's BACE design):

| condition | `linpred` | depends on |
|---|---|---|
| `phylo_MAR` | `-1.5·scale(z)`, `z ~ MVN(0, R)` | a phylogenetically autocorrelated latent — observed pattern, unobserved cause |
| `trait_MAR` | `-1.5·scale(env)` if a covariate exists, else `0` | the `env` covariate |
| `trait_MNAR` | `-1.5·scale(y)` | the value of `y` itself |

**Design limitation — read before interpreting `trait_MAR`.** In this
single-trait DGP, the only candidate "other observed variable" for a
MAR mechanism is the `env` covariate, and `env` exists *only* for
`nonlinear_cov`. For `bm_strong`, `ou_strong`, and `weak_signal` the
`trait_MAR` condition has nothing to depend on and falls back to
`linpred = 0`, i.e. **pure MCAR**. The conditions genuinely exercised
are therefore: phylo-MAR, MCAR, and value-dependent MNAR for all four
scenarios; plus a genuine covariate-driven MAR for `nonlinear_cov`
only. Dan's multi-trait BACE simulation does not have this limitation
(missingness in `y` can depend on observed `x1`…`x4`); the single-trait
adaptation here trades that realism for a cleaner contrast and direct
comparability across `n` and across the four scenarios.

## 4. Methods compared

| method | description |
|---|---|
| `column_mean` | impute every missing cell with the observed-cell mean. The null model. |
| `bm_kriging` | pigauto's internal `bm_impute_col()` — conditional-MVN BM kriging on `R`. The comparator of record. |
| `pigauto_default` | full `impute()` pipeline: BM baseline + gated attention GNN, `lambda_mode = "fixed_1"`, 300 epochs, patience 5. |
| `pigauto_bayes` | as `pigauto_default` but `lambda_mode = "bayes"`. **240-rep pilot only** — dropped from the 600-rep sweep (see §5). |

## 5. Design and replication

Two sibling runs share the §2–§3 generators verbatim (the full sweep
adds a `+200000` offset to the rep-seed arithmetic so its seeds are
disjoint from the pilot's):

| run | scenarios | conditions | `n` | sims/cell | methods | total reps |
|---|---|---|---|---|---|---|
| 240-rep pilot | 4 | 3 | {200, 500} | 10 | 4 | 240 |
| 600-rep sweep | 4 | 3 | {500} | 50 | 3 | 600 |

The 600-rep sweep is the Dan-comparable run: `n = 500` and 50 sims per
cell match the per-cell replication of `BACE/dev/02_benchmark_simulated_full.R`,
so the two studies can be merged into one figure.

`pigauto_bayes` was dropped from the sweep because the pilot showed it
within one Monte Carlo SE of `pigauto_default` in every one of the 24
pilot cells — it carried no information and doubled pigauto wall time.

Execution is strictly serial: torch on the Apple MPS backend is not
fork-safe, so `mclapply` with more than one worker segfaults. A per-rep
RDS checkpoint (`results_partial.rds`) makes the ~10 h sweep
crash-resilient.

## 6. Performance measure

Per replicate, per method:

```
z_rmse = sqrt(mean((pred - truth)^2)) / sd(y_observed)
```

z-scoring by the observed-cell SD makes the metric comparable across
scenarios with different signal variance, and pins the column-mean
baseline near 1.0 by construction (it is exactly 1.0 in expectation
under MCAR). A method beats column-mean iff `z_rmse < 1`. Per-cell
means are reported with Monte Carlo SE = `sd / sqrt(n_sims)`.

## Addendum: conformal-coverage remeasurement

The original mechanism sweep was an RMSE study. A follow-up coverage
remeasurement (`useful/MEMO_2026-05-21_coverage_remeasurement.md`) used
the same mechanism families at n = 500, 30% missing, 120 replicates, and
measured three different quantities:

- `cov_conformal`: pigauto's actual split-conformal interval coverage
  from `pred$conformal_lower` / `pred$conformal_upper`.
- `cov_drawband`: the older mixed-type sweep's 2.5/97.5% band of 20
  imputation draws.
- `cov_bm`: BM-kriging's analytic Gaussian interval, as a reference.

The distinction matters. `cov_drawband` is not pigauto's conformal
interval; it is a finite-draw MI band and systematically understates
coverage relative to the actual split-conformal interval.

| scenario | condition | cov_conformal | cov_drawband | cov_bm |
|---|---|---|---|---|
| bm_strong | MCAR | 0.937 +/- 0.015 | 0.872 +/- 0.020 | 0.931 +/- 0.007 |
| bm_strong | phylo_MAR | 0.952 +/- 0.010 | 0.890 +/- 0.010 | 0.949 +/- 0.005 |
| bm_strong | trait_MAR | 0.934 +/- 0.018 | 0.853 +/- 0.024 | 0.924 +/- 0.008 |
| bm_strong | trait_MNAR | 0.878 +/- 0.022 | 0.799 +/- 0.021 | 0.928 +/- 0.008 |
| weak_signal | MCAR | 0.966 +/- 0.007 | 0.898 +/- 0.014 | 0.918 +/- 0.010 |
| weak_signal | phylo_MAR | 0.948 +/- 0.012 | 0.876 +/- 0.017 | 0.954 +/- 0.009 |
| weak_signal | trait_MAR | 0.951 +/- 0.008 | 0.893 +/- 0.013 | 0.933 +/- 0.010 |
| weak_signal | trait_MNAR | 0.815 +/- 0.020 | 0.690 +/- 0.022 | 0.890 +/- 0.014 |

Reading:

- Under MCAR, where split-conformal exchangeability is the intended
  reference regime, `cov_conformal` sits near the nominal 0.95 target.
- Under phylo-MAR and the degenerate-MCAR `trait_MAR` cells, coverage is
  also near nominal in this design.
- Under value-dependent MNAR, undercoverage is real rather than a
  measurement artifact: the observed calibration residuals no longer
  represent the missing target cells.

## 7. Results — 600-rep sweep (n = 500, 50 sims/cell)

`Δ = z_rmse(pigauto) − z_rmse(bm_kriging)`. Negative = pigauto better.
Verdict: **PIG** if `Δ < −0.05`, **BM** if `Δ > +0.05`, else tie.
`Δ_se` is the conservative `sqrt(se_pig² + se_bm²)` (a paired SE over
shared `sim_id` would be tighter).

| scenario | condition | column-mean | bm-kriging | pigauto | Δ | Δ_se | verdict |
|---|---|---|---|---|---|---|---|
| bm_strong | phylo_MAR | 1.139 | **0.687** | 0.691 | +0.004 | 0.038 | tie |
| bm_strong | trait_MAR* | 1.014 | 0.589 | **0.584** | −0.005 | 0.031 | tie |
| bm_strong | trait_MNAR | 1.549 | **0.744** | 0.788 | +0.044 | 0.063 | tie |
| nonlinear_cov | phylo_MAR | 1.102 | 0.990 | **0.947** | −0.043 | 0.027 | tie |
| nonlinear_cov | trait_MAR | 1.007 | 1.041 | **0.951** | −0.091 | 0.022 | **PIG** |
| nonlinear_cov | trait_MNAR | 1.508 | **1.302** | 1.328 | +0.025 | 0.046 | tie |
| ou_strong | phylo_MAR | 1.033 | 1.157 | **1.025** | −0.132 | 0.019 | **PIG** |
| ou_strong | trait_MAR* | 1.032 | 1.160 | **1.021** | −0.139 | 0.021 | **PIG** |
| ou_strong | trait_MNAR | 1.542 | 1.510 | **1.502** | −0.008 | 0.030 | tie |
| weak_signal | phylo_MAR | 1.014 | 1.234 | **1.013** | −0.220 | 0.019 | **PIG** |
| weak_signal | trait_MAR* | 1.003 | 1.204 | **1.004** | −0.200 | 0.014 | **PIG** |
| weak_signal | trait_MNAR | 1.565 | 1.667 | **1.566** | −0.101 | 0.021 | **PIG** |

`*` `trait_MAR` = MCAR for these cells (no covariate; see §3).

**Tally: 6 PIG, 6 ties, 0 BM.** Across all 12 cells at n = 500, pigauto
never loses materially to BM kriging. `nonlinear_cov/phylo_MAR`
(Δ = −0.043) is a near-win just inside the tie band.

The 240-rep pilot agrees cell-for-cell in sign and shows two BM-win
cells — `bm_strong/trait_MNAR` (Δ = +0.232) and
`nonlinear_cov/trait_MNAR` (Δ = +0.088) — both at **n = 200 only**. At
n = 500 those same cells are ties (+0.044, +0.025). They were
small-sample artefacts of 10-rep cells, and the pilot's own ±0.23 SE
on the first already flagged them as non-significant.

## 8. Interpretation

The headline is one sentence: **across four evolutionary processes and
three missingness conditions, using pigauto in place of BM kriging is
never materially worse, and in half the cells it is materially
better.** That is exactly the safety-floor guarantee the gated
architecture is designed to provide.

The *reason* pigauto wins differs by scenario, and only one of the
three win-types is the GNN doing genuinely novel work:

- **`weak_signal` (Δ ≈ −0.10 to −0.22).** BM kriging mistakes
  90%-noise for phylogenetic signal and over-shrinks toward
  phylogenetically-close tips, landing *worse than column-mean*.
  pigauto's gate closes (§9) and pigauto reduces to ≈column-mean. The
  win is the **safety floor refusing a bad baseline**, not GNN
  learning.

- **`ou_strong` (Δ ≈ −0.13).** The process is OU-like; BM kriging
  over-commits to deep-time conservation. pigauto beats *both*
  baselines — the gate is mostly closed but the minority of open-gate
  reps supply a real correction. This is "do not over-trust a
  misspecified baseline", with a sliver of genuine GNN contribution.

- **`nonlinear_cov` (Δ ≈ −0.09, genuine `trait_MAR`).** This is the
  scenario where the gate opens most often (§9). The GNN exploits the
  `env` covariate — a `sin(2·env)` term no phylogeny-only method can
  represent — and beats column-mean *and* BM. This is the GNN earning
  its keep.

- **`bm_strong` (tie).** BM is the true process; pigauto ties it. The
  gate is closed in most reps and the blend is harmless in the rest
  (§9). The safety floor holds without a capability cost.

On the missingness axis: `trait_MNAR` uniformly raises difficulty —
the observed sample is value-biased, so every method's z-RMSE rises.
On `bm_strong`, column-mean degrades most (z-RMSE 1.01–1.14 under
phylo-MAR/MCAR → 1.55 under MNAR) while the phylogenetic methods
absorb it far better (bm-kriging 0.59–0.69 → 0.74; pigauto 0.58–0.69 →
0.79). MNAR never flips a verdict — pigauto's standing relative to BM
is stable across phylo-MAR, MCAR, and MNAR.
Note the §3 caveat: genuine three-way mechanism variation is only
tested for `nonlinear_cov`; for the other scenarios the claim is
stability across phylo-MAR, MCAR, and MNAR.

## 9. Calibrated-gate diagnostics

The per-rep mean calibrated gate is **bimodal**: in any given rep it is
either ≈0 (GNN off) or substantially open. The per-cell *mean* gate is
pulled up by the open-gate minority and misrepresents the typical rep,
whose gate is 0. The honest summary is the fraction of reps with the
gate effectively closed (`gate < 0.02`):

| scenario | mean gate | median gate | reps with gate closed |
|---|---|---|---|
| `weak_signal` | 0.009 | 0.000 | 96–100% |
| `ou_strong` | 0.076 | 0.000 | 82–86% |
| `bm_strong` | 0.187 | 0.000 | 54–72% |
| `nonlinear_cov` | 0.222 | 0.000 | 56–58% |

The ordering is mechanistically sensible: the gate opens most often
exactly where the GNN has something to contribute (`nonlinear_cov` has
the covariate), and is essentially always shut where it does not
(`weak_signal` has no learnable structure). Crucially, even in
`bm_strong` — where the gate opens in ~30–45% of reps — pigauto still
ties BM, which means the validation-set gate calibration is admitting
the GNN only when its blended prediction is at least as good as the
baseline. The gate is not a binary on/off switch tuned once; it is
re-calibrated per fit, and the bimodality is the calibration grid
landing at either end depending on whether the held-out residuals
favour the GNN.

## 10. Limitations

1. **`trait_MAR` degeneracy (§3).** Three of four scenarios test MCAR
   under the `trait_MAR` label. The genuine covariate-MAR result rests
   on `nonlinear_cov` alone. A multi-trait DGP (as in Dan's BACE sim)
   would close this gap.
2. **Single trait, continuous only.** The sweep says nothing about
   pigauto's mixed-type machinery (binary, categorical, ordinal,
   count, proportion). It tests the BM-baseline + GNN path for one
   continuous trait.
3. **One `n`, one missing rate in the headline sweep.** n = 500, 30%
   missing. The pilot adds n = 200; neither sweep varies the missing
   rate. Generalisation to n ≫ 500 or to heavier missingness is not
   established here.
4. **Synthetic trees and processes.** `rcoal` trees and exact
   BM/OU/covariate DGPs. Real phylogenies and real traits depart from
   all of these; the companion real-data benches (AVONET, BIEN,
   Delhey) are the external-validity check, not this study.
5. **`Δ_se` is unpaired.** Reported SEs ignore the shared `sim_id`
   pairing between methods within a cell, so they are conservative;
   the true precision on Δ is somewhat better than tabulated.
6. **`pigauto_default` only in the sweep.** `lambda_mode = "bayes"`
   was dropped after the pilot showed equivalence; the sweep does not
   re-confirm that at 50 sims/cell.

## 11. Reproducibility

```r
# 240-rep pilot (n in {200, 500}, 10 sims/cell)
Rscript script/sim_bench/overnight_2026_05_19_mechanisms.R

# 600-rep Dan-comparable sweep (n = 500, 50 sims/cell)
Rscript script/sim_bench/pigauto_full_sweep.R

# Combined memo with per-cell tables (this study's auto-generated companion)
Rscript script/sim_bench/build_mechanism_sweep_memo.R
```

Seeds are deterministic functions of `(scenario, n, mechanism, sim_id)`
inside each driver; reruns reproduce the RDS bit-for-bit on the same
torch/BLAS build. Wall time on the development laptop (Apple MPS,
serial): pilot ≈ 9.5 h, sweep ≈ 11 h (`pigauto_default` ≈ 64 s/rep at
n = 500; `bm_kriging` and `column_mean` are sub-millisecond).
