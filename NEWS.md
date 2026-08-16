# pigauto 0.10.0.9000 (development)

GitHub-dev after CRAN pigauto 0.10.0 (Date/Publication 2026-07-30).
This is not a CRAN tarball. BACE is still not on CRAN. The next CRAN
cut must drop Suggests `BACE` again or wait until BACE is on CRAN.

## Feature: `conformal_method = "mondrian"` (B2) -- locality-stratified conformal intervals

`fit_pigauto()` / `impute()` gain a new `conformal_method` option,
`"mondrian"`, alongside the existing `"split"` (default) and `"bootstrap"`.
It targets the exchangeability failure diagnosed in the mechanism-coverage
campaign (`docs/dev-log/2026-08-16-mechanism-coverage-results.md`): a single
global conformal quantile per trait is calibrated on the observed
complement, which sits systematically closer (in cophenetic distance) to
other observed cells than genuinely-missing cells do -- worst under
clade-structured missingness (MAR_phylo, -3.4pp coverage vs MCAR at
n=300/1000).

`"mondrian"` conditions calibration on phylogenetic sampling locality: each
validation cell's locality is the mean cophenetic distance to its 5 nearest
species with an observed value for that trait, cells are split into 2
strata at the median locality, and a split-conformal quantile is computed
per stratum. Cells in the far (undersampled) stratum get the wider
quantile at prediction time -- intervals widen in poorly-sampled clades,
which is the point. Traits with fewer than 19 validation residuals in
either stratum fall back to the global `"split"` score (19 is the
smallest stratum whose achievable coverage ceiling n/(n+1) reaches 0.95;
smaller strata are arithmetically guaranteed to undercover). Single-obs only:
`fit_pigauto()` errors immediately for multi-observation data, since
locality is defined at species level and multi-obs validation cells are
per-observation.

## Silent-fallback honesty pass (2026-08-16)

Four places where the pipeline silently did less than the user asked are
now loud, plus one documentation debt:

- `lambda_mode = "estimate"/"cv"/"bayes"` + covariates: the covariate-aware
  GLS path (`bm_impute_col_with_cov()`) has no lambda argument, so those
  modes silently ran at lambda = 1. `fit_baseline()` now warns once.
- Covariates + joint baselines (P1-8, detect-and-message): the joint MVN
  and threshold-joint (Rphylopars) baselines have no covariate design, so
  user covariates never reached the BASELINE (they do reach the GNN).
  `fit_baseline()` now warns once when a joint path fires with covariates
  supplied. Threading covariates through the joint solver is future work.
- Small validation sets: at n_val < 19 the 95% split-conformal target is
  arithmetically unreachable (achievable ceiling n_val/(n_val+1)), not
  merely "noisy" as the old warning said. The warning now states the
  ceiling and fires for the previously-silent 10-18 cell band too
  (binary/categorical excluded — no conformal scores there).
- `impute()` now exposes `conformal_split_val` explicitly (previously
  reachable only via `...`, i.e. undocumented at the primary entry point).

Documentation debt: `lambda_mode` values `"cv"` (per-column CV-selected
lambda) and `"bayes"` (Bayesian model-averaged prediction over a lambda
grid, with posterior-mean lambda as a diagnostic; the SE combines
within- and between-lambda variance) shipped in PR #109 but were never
recorded in NEWS. The 2026-08-16 lambda-attribution campaign
(`docs/dev-log/2026-08-16-lambda-attribution-results.md`) measured
`"bayes"` as the mode that eliminates the ML boundary collapse the
v0.10.0 "Known limitation" section describes, superseding that section's
"v0.11 CV-lambda" remedy plan.

## Fix: phylo-signal gate was a silent no-op in multi-obs mode (P1-12)

`compute_phylo_signal_per_trait()` indexed `species_names` (one entry per
species) with an observation-length logical drawn from `X_scaled` (one row
per observation). In multi-obs mode this produced NA tip names,
`ape::keep.tip()` errored, the caller's `tryCatch` swallowed it to `NA`, and
`phylo_signal_gate` therefore never fired for any multi-obs trait —
silently, with no warning. Signal is now computed on species means, matching
how `fit_baseline()` aggregates multi-obs input. Verified against the
pre-fix code: a BM-simulated trait returned `NA` before and a finite
lambda > 0.5 after. Single-obs behaviour is unchanged, and with one
observation per species the two paths agree exactly.

Multi-obs users with `phylo_signal_gate = TRUE` may see gate decisions that
previously could not occur. That is the documented behaviour finally taking
effect, not a change of policy.

## Fix: zi_count observed zeros were excluded from val/test splits (P1-9)

`make_missing_splits()` treated a (species, trait) cell as observed only when
**all** its latent columns were non-NA. `zi_count` encodes an observed zero as
gate = 0 with the log1p-z magnitude column `NA` — the NA is data ("this count
is zero"), not missingness — so every observed zero failed that test and could
never be drawn into validation or test. Because zero-inflation is the defining
feature of the type, this silently restricted all `zi_count` gate calibration,
conformal scoring and evaluation to the non-zero subset. Observation status for
`zi_count` is now decided by the gate column alone.

Measured on a 60-species fixture with 23 observed zeros: **0** of them could be
held out before the fix; they are eligible after. Any previously reported
`zi_count` validation or test metric was computed on non-zeros only and should
be re-measured.

## Fix: full-validation BM floor in gate calibration (#157)

Gate calibration accepted blends from half-set comparisons (or half-B
fallback corners) that could be materially worse than the pure BM
baseline on the full validation set — for continuous-family traits no
full-val BM comparison existed at all. `calibrate_gates()` now checks
both pure corners on the full validation set for every family; the
continuous-family BM comparison is margin-based (`cal_min_rel_gain / 2`,
default 1%) so near-ties survive (the 2026-04-29 AVONET Mass case)
while material losses close the gate. A second, identical floor runs in
`fit_pigauto()` on the post-refine calibrated surface — the one
`predict()` delivers — because a gate that wins on the pre-refine
surface can lose to pure BM after the calibrated refine pass (+2.0%
measured on the #157 fixture); a gate closed there costs nothing since
the delta contribution is zeroed. The strict val-floor test now
evaluates on the calibration-matched surface
(`.mask_observed_idx = c(val_idx, test_idx)`), matching `evaluate()` /
`cross_validate()` / `impute()` / `report()` — a plain `predict()` pins
each held-out cell's own truth as input context, a surface no genuinely
missing cell has in production.

## New (opt-in): BACE baseline wrapper restored, with proper-MI draws

Restores `fit_baseline_bace()` (removed from the v0.10.0 CRAN surface in
`b615579` because BACE is not on CRAN) and adds `final_imp = FALSE`
(default) plus `n_final = 15L`. With `final_imp = TRUE` the wrapper
appends `BACE::bace_final_imp()` to the `BACE::bace_imp()` chain it
already runs and builds `mu` / `se` from that function's `n_final` final
datasets instead of the chain datasets. `BACE` is restored in Suggests
on this branch; there is no Remotes field. Install BACE from the
standalone tree at `@ce8bc87`.

Why the opt-in path is worth having: the chain datasets are successive
sweeps of one chained-equations chain, so they are autocorrelated by
construction and still carry convergence transient. The between-dataset
SD pigauto reports as `se` on that path is therefore a dispersion
summary with no coverage guarantee. `bace_final_imp()`'s runs each start
independently from the converged chain, so they are proper
multiple-imputation draws and the resulting `se` is a between-imputation
SD in the Rubin (1987) sense. Expect it to be substantially larger.

This is a choice about which BACE datasets pigauto summarises. It is
**not** a correction to the default path's arithmetic, and in particular
it is unrelated to the imputed-as-observed defect fixed upstream in
BACE — that defect lived in `bace_final_imp()`, which pigauto has never
called until now, so no pigauto output was ever affected by it.

`final_imp = TRUE` costs `n_final` extra MCMC fits per trait on top of
the chain. The default of 15 matches the BACE simulation study; BACE's
own default is 50, and small values undercover. The default remains
`FALSE`. If the installed BACE is too old to export `bace_final_imp()`,
`final_imp = TRUE` errors with an upgrade hint rather than silently
falling back.

The final phase is less robust than the chain phase — each draw refits
MCMCglmm and can hit a singular mixed-model equation on data the chain
handled fine. In that case `final_imp = TRUE` raises an error naming
the failure and pointing back at `final_imp = FALSE`, rather than
quietly returning chain averages to a caller who asked for proper MI
draws.

No pigauto-versus-BACE performance claim is made here.

## Documentation: pigauto's OVR categorical path vs BACE's default

`fit_ovr_categorical_fits()` no longer describes its one-vs-rest
decomposition as "the OVR strategy BACE uses". BACE's
`ovr_categorical` default became `FALSE` (true multinomial) in
2026-08. pigauto's OVR is a different estimator — K threshold-joint
Rphylopars fits, used so that Rphylopars stays well-conditioned — and
that motivation is unchanged. Comment/roxygen only; no estimator
change.

## Fix: covariate rows follow species / observation identity

`preprocess_traits()` / `impute()` reordered trait rows to tree tip
order but only checked that `covariates` had the same `nrow`. Shuffled
trait tables with matching-nrow environmental data silently mispaired
species and climate in the GNN covariate tensor. Covariates are now
aligned to the same species (single-obs rownames) or observation
(multi-obs `species_col`) order as `X_scaled`. Unmatched names error.
Unnamed covariates still follow the input-row permutation of `traits`.
GNN math is unchanged.

## Fix: zi_count conformal multiple imputation

`multi_impute(draws_method = "conformal")` now scores `zi_count` traits on the log1p-z magnitude latent (held-out non-zeros) and draws count|nonzero from that score, or from the magnitude latent SE when scores are unavailable. The zero-inflation gate stays Bernoulli from P(nonzero). Previously scores were skipped and draws used `pred$se` of E[X], mixing gate uncertainty into the count SD. `draws_method = "mc_dropout"` is unchanged. Reported `$se` of expected count is not rewritten (P1).

## Fix: per-tip BM SE on default K≥2 / threshold-joint path

Under `max_iter = 0` (shipped default; cross-trait EM is not restored),
the Henderson init used for K≥2 joint MVN and mixed threshold-joint
baselines recycled one empirical SD for every missing tip. Missing-cell
SE is now BM-style σ√(1−h_i) from the same sparse Q already used for μ.
Binary threshold decode (`pnorm(μ / sqrt(1+se²))`) therefore sees
per-tip uncertainty. Conformal intervals are unchanged (they do not
consume baseline `se`).

## Fix: train / calibrate / predict symmetry for the GNN path

Three related wiring bugs in the ResidualPhyloDAE fit path (review
2026-08-07). Defaults and the public API are unchanged; behaviour is
more honest about what the model sees at train and cal time.

- **Held-out masking (train):** validation and test cells no longer
  appear as truth in the DAE input during training. They are replaced
  with the baseline prediction (`mask_heldout_with_baseline()`;
  `model_config$train_mask_heldout = TRUE`).
- **Phylo-signal gate without safety floor:** when
  `phylo_signal_gate` triggers and `safety_floor = FALSE` (so no mean
  baseline is available), the override falls back to pure BM
  (`r_cal_bm = 1`) with a warning, instead of the mean corner that
  would become latent 0 at predict time.
- **Refine steps at calibration / conformal:** gate calibration and
  conformal residual scoring now run the same iterative
  `refine_steps` refine loop that `predict()` uses
  (`pigauto_refine_forward()`; `model_config$cal_refine_steps`).

# pigauto 0.10.0

## Experimental analysis-aware multiple-imputation backend

The new `multi_impute_analysis()` backend requires the substantive model
before imputation and is deliberately narrow: one incomplete continuous
covariate under MAR, with proper Bayesian normal-regression MI for Gaussian
`lm`, `smcfcs` for binomial-logit `glm`, and `jomo.smc` for a Gaussian `lmer`
with one random intercept. It pools fixed effects only.

A warning-free 6,000-task package-level campaign at clean source SHA
`2e3809d` passed all 24 predeclared fixed-effect cells. Added bias was
0.00086-0.05042 oracle SD, pooled-SE/empirical-SD ratios were 0.942-1.030,
coverage was 93.9%-96.3%, and all results were finite. Engine warnings are
fatal, so no inference-ready object is returned after a rejection or
convergence warning. This validates the backend only within its documented
narrow scope.

The earlier conformal-width and Brownian/MC-dropout methods each passed 0/12
downstream fixed-effect cells in the 3,000-fit package campaign. PMM also failed
its redesign pilot. Those methods remain experimental prediction diagnostics
and are unsupported for downstream inference. Tree uncertainty was outside the
analysis-aware campaign and remains unsupported by the new backend.

## New (opt-in): Pagel's lambda Brownian-motion baseline

Adds first-pass support for Pagel's lambda < 1 in the BM baseline.
Spec: `specs/2026-05-18-pagel-lambda-baseline-design.md` (APPROVED).

New parameters:

- `fit_pigauto(..., lambda_mode = c("fixed_1", "estimate"))`
- `fit_baseline(..., lambda_mode = c("fixed_1", "estimate"))`
- kernel-layer `bm_impute_col(..., lambda)` accepts numeric in `[0, 1]`
  OR the string `"estimate"` (in-house ML via REML profile + 11-point
  grid scan + local optim)

Default `lambda_mode = "fixed_1"` is **bit-identical to v0.9.x** — no
silent behaviour change for existing callers.

What lambda < 1 does mathematically:

```
R(lambda) = lambda * R + (1 - lambda) * I
```

shrinks the off-diagonal of the phylogenetic correlation matrix
toward zero. lambda = 1 is full BM kriging; lambda = 0 reduces the
GLS predictor to the grand mean. Intermediate lambda allows partial
phylogenetic shrinkage for traits with weak phylo signal.

**No phylopars dependency on the lambda path.** Both `transform_tree_pagel`
(internal scaler with Pagel terminal-edge compensation) and
`ml_lambda_for_col` (per-column ML) are in-house. When
`lambda_mode = "estimate"`, the dispatcher forces per-column BM
(joint MVN / threshold-joint / OVR are gated off because phylopars
`model = "BM"` is silently lambda = 1).

## Known limitation (honest framing)

On the v0.9.4 BIEN h2h CI bench (n=2000, 30% MCAR, m=20):

- height_m: ML lambda_hat ≈ 0.95, z-RMSE 0.73 (close to RMSE-optimal)
- leaf_area: ML lambda_hat ≈ 0.83, z-RMSE 0.92 (slight regression vs lambda=1's 0.89)
- seed_mass: ML lambda_hat ≈ 0.99, z-RMSE 0.60 (essentially unchanged)
- wood_density: ML lambda_hat ≈ 0.79, z-RMSE 0.79 (slight regression vs lambda=1's 0.77)
- **sla**: ML lambda_hat ≈ 0.005 (boundary!), z-RMSE 0.98 (regresses from lambda=1's 0.91)

The sla case is the well-known **weak-phylo-signal pathology** of ML
for Pagel's lambda: the profile likelihood becomes flat for small
lambda when phylogenetic signal is weak, finite-sample noise
dominates, and the point estimator can collapse to the boundary. The
math is correct; the predictor is just suboptimal.

BACE's MCMCglmm achieves z-RMSE ≈ 0.69 on sla because Bayesian model
averaging integrates over the variance-component posterior rather
than relying on a point ML estimate.

**Practical guidance:** keep `lambda_mode = "fixed_1"` (the default)
unless you've verified on your data that ML lambda doesn't pick
degenerate boundary values. The spec's hard target ("sla z-RMSE
≤ 0.75") is **not achieved by this release**; it's a v0.11 target
via cross-validation lambda selection (spec to follow).

## Tests

- 16 new tests in `tests/testthat/test-pagel-lambda.R` (31 expectations)
- All 54 existing `bm-internal` tests pass — back-compat verified

## Simulation evidence: mechanism-axis robustness sweep

A two-part Monte Carlo study (no code change) testing whether pigauto's
gated GNN holds up against the plain BM-kriging baseline under realistic
non-random missingness. Drivers:
`script/sim_bench/overnight_2026_05_19_mechanisms.R` (240 reps) and
`script/sim_bench/pigauto_full_sweep.R` (600 reps, Dan-comparable n=500,
50 sims/cell). Methodology writeup:
`script/sim_bench/mechanism_sweep_methodology.md`.

Design: 4 process scenarios (`bm_strong`, `ou_strong`, `nonlinear_cov`,
`weak_signal`) × 3 missingness conditions (`phylo_MAR`, `trait_MAR`,
`trait_MNAR`, DEP_STRENGTH = 1.5 per Penone et al. 2014) × 30% missing.
Caveat: in this single-trait DGP, `trait_MAR` depends on the `env`
covariate, which exists only for `nonlinear_cov`; for the other three
scenarios `trait_MAR` has no variable to depend on and degenerates to
MCAR. The conditions genuinely exercised are therefore phylo-MAR, MCAR,
and value-dependent MNAR for all four scenarios, plus a genuine
covariate-MAR for `nonlinear_cov` only.

Result on the 600-rep n=500 sweep — across all 12 cells: **6 material
wins for pigauto, 6 ties, 0 material losses** vs BM kriging (z-RMSE,
material = |Δ| > 0.05). The win mechanism is scenario-specific and only
partly attributable to the GNN:

- `weak_signal`: the calibrated gate closes to ≈0; pigauto falls back to
  column-mean while BM kriging mis-extrapolates 90%-noise as phylo signal
  (Δ up to −0.22). This is the safety floor, not GNN learning.
- `ou_strong`: gate stays low (≈0.07–0.09) but pigauto beats both
  baselines by avoiding BM's process-misspecification (Δ ≈ −0.13).
- `nonlinear_cov`: gate most open (≈0.18–0.25); the GNN + environmental
  covariate genuinely beat column-mean *and* BM (Δ ≈ −0.09).
- `bm_strong`: gate harmless; pigauto ties BM (the correct answer).

The missingness condition scales difficulty (`trait_MNAR` is uniformly
hardest) but never flips a verdict: pigauto's per-scenario standing
relative to BM kriging is stable across phylo-MAR, MCAR, and MNAR. The
2 apparent losses seen in the smaller 10-sim pilot were both
n=200/`trait_MNAR` and did not replicate at n=500.

## Independent verification: `use_trait_attention` redundant at scale

External 60-replicate ablation (b1805, [#106](https://github.com/itchyshin/pigauto/issues/106), code in [#116](https://github.com/itchyshin/pigauto/pull/116))
using a corrected multi-trait DGP that fires the joint-MVN baseline.
N=2000, 5 seeds per cell, three arms — including the "lazy-optimiser
disarmed" `pigauto_ON_L0` configuration (`baseline_mu` masked, `lambda_shrink = 0`).

Result (z-RMSE, averaged across all scenarios):

| arm             | z_rmse | coverage | width |
|-----------------|--------|----------|-------|
| `column_mean`   | 1.291  | —        | —     |
| `bm_kriging`    | 1.163  | —        | —     |
| `pigauto_OFF`   | 1.038  | 0.887    | 2.65  |
| `pigauto_ON`    | 1.056  | 0.884    | 2.67  |
| `pigauto_ON_L0` | 1.057  | 0.887    | 2.68  |

(Measured before the 2026-08 GNN train/cal-symmetry fix, which removed
held-out truth from the training-time input context; not yet
re-verified on the corrected pipeline.)

At BIEN scale (N=2000) the joint-MVN baseline already encodes the
cross-trait Σ structure; the within-row attention mechanism injects
noise and slightly widens conformal intervals. Even with the lazy-
optimiser trap disarmed (`pigauto_ON_L0`), the network regresses
against the baseline rather than finding new structure.

`use_trait_attention` stays default-`FALSE`. Kept as opt-in for
small-N regimes where the baseline hasn't converged (b1805's local
smoke test at N=300 showed it dropping validation loss before the
joint-MVN catches up at larger N).

# pigauto 0.9.4 (2026-05-18)

Wires the existing WorldClim bioclim covariate pipeline (shipped in
v0.9.1.9005) into the BACE-compatible head-to-head CI bench for the
BIEN plant dataset. Closes the gap surfaced by the v0.9.3 autoresearch
sweep — at that scale, GNN architecture knobs were exhausted; the
remaining improvement had to come from environmental covariates,
not capacity.

## Result

BIEN pooled per-trait log-RMSE on the h2h CI bench (n=2000 species,
30% MCAR, m=20 imputations):

```
baseline  (full 19,109-pool, no cov, gate ON):  8.286
control   (filtered 3,450-pool, no cov, gate OFF): 7.455
treatment (filtered pool, WC covariates, gate OFF): 7.236
```

Pool-filter effect = −0.83. Net covariate effect = −0.22.
Combined gain vs v0.9.3 baseline: **−1.05 (−12.7%)**.

Per-trait gates with covariates active:

```
height_m     = 0.68   (raw |r| vs bioclim = 0.64)
leaf_area    = 1.00   (raw |r| vs bioclim = 0.41)
wood_density = 0.40   (raw |r| vs bioclim = 0.16)
sla          = 0      (raw |r| vs bioclim = 0.27; gate-closed)
seed_mass    = 0      (raw |r| vs bioclim = 0.43; gate-closed)
```

Two of five traits (`sla`, `seed_mass`) remain gate-closed even with
covariates active. `sla` is the expected case (low raw correlation).
`seed_mass` is anomalous — its raw bioclim correlation is comparable
to `height_m`'s, but the val-set gate calibration finds no benefit.
Diagnosis to follow.

## Changes

- `script/gha/_bace_compat.R`: BIEN config now sets
  `covariate_cols` to the 38 bioclim columns (19 vars × median + IQR),
  `phylo_signal_gate = FALSE`, and `external_cov_rds` pointing to the
  new bioclim cache. `.run_bench_bace_compat()` learned an
  `external_cov_rds` join that cbinds covariates onto `traits_df`
  before subsetting and filters out species without complete bioclim
  coverage. `cov_cols` AND external-cov column names are removed
  from the default trait-subset so they don't get masked + imputed
  as traits.
- New cache: `useful/bace_data_snapshot/data/bien_worldclim_covariates.rds`
  (19,109 species × 38 cols; 3,450 with complete coverage).

## Negative results documented for the trail

The W-series sweep tested six tweaks on top of the W0 (treatment)
baseline; all regressed by 0.09–0.12, except W6 which was tied
within noise:

- W3 lambda_shrink 0.03 → 0.01 (+0.12)
- W4 lambda_gate 0.01 → 0.001 (+0.09; leaf_area gate collapsed)
- W5 median-only cov (drop 19 IQR cols) (+0.12)
- W6 phylo_signal_gate FALSE → TRUE (+0.01, tied)
- Run A: WorldClim + use_trait_attention=TRUE (+0.11)

W6's tie shows the `phylo_signal_gate = FALSE` recommendation
from `specs/2026-04-24-worldclim-covariates-design.md` §5.4 is
unnecessary when WC covariates are active — the calibrator handles
the gating either way. The config keeps `FALSE` for stability but
the constraint is softer than the spec implied.

# pigauto 0.9.3 (2026-05-18)

Small opt-in feature release from an autoresearch sweep on the BIEN
plant bench. New: within-row cross-trait self-attention in the GNN
(`use_trait_attention`, default `FALSE`). Eighteen knob-tweak
experiments + three structural variants on BIEN n=2000 were run
between v0.9.2 and v0.9.3; none improved BIEN pooled RMSE
meaningfully over the default architecture. The cross-trait attention
path is shipped as opt-in because it is structurally distinct and may
help on datasets with strong within-row functional coupling that the
joint MVN / threshold-joint baseline cannot capture; it did NOT help
on BIEN (the baseline already encodes Σ for the dominant trait types).

## New (opt-in)

- `fit_pigauto()` and `multi_impute()` accept `use_trait_attention =
  FALSE` (default), `n_trait_heads = 2L`, `trait_embed_dim = 32L`.
  When `TRUE`, the GNN's encoder receives an additional pooled
  trait-context feature built by per-trait token embedding +
  one multi-head self-attention block over the trait sequence.
  Backward-compatible: saved fits without the field default to
  `FALSE` and reconstruct identically.

## Documentation

- New article on the pkgdown site: **"GNN architecture and the math
  behind pigauto"** (`vignettes/gnn-architecture.Rmd`). Documents the
  end-to-end data flow, the baseline dispatcher (BM / LP / joint MVN /
  threshold-joint / OVR), the transformer-block GNN with B2 rate-aware
  attention and the new B3 within-row cross-trait attention, the gate-
  and-calibrate safety mechanism, conformal prediction intervals, the
  three uncertainty quantification mechanisms, and the two-step
  Nakagawa & de Villemereuil (2019) tree-uncertainty workflow.

## Autoresearch findings (BIEN n=2000)

Negative results, documented for the trail:

- More capacity hurts: `n_gnn_layers` 2→4, `hidden_dim` 64→128,
  `ffn_mult` 4→8, `n_heads` 4→16, and all combos regress (+0.02 to
  +0.05 on pooled RMSE).
- Less capacity hurts equally: `n_heads` 4→2, `hidden_dim` 64→32,
  legacy attention path also regress.
- Regularization tweaks regress: `dropout` 0.10→0.20,
  `lambda_gate` 0.01→0.001, `gate_cap` 0.8→1.0.
- Refinement / corruption tweaks regress: `refine_steps` 8→16,
  `corruption_rate` 0.55→0.30.
- Within-row cross-trait attention (this release's new feature)
  also regresses on BIEN at 32/64 embed and 2/4 heads — the joint
  MVN baseline already captures Σ for continuous traits.

Marginal positive: `n_heads` 4→8 in the CI bench config (-0.002,
at-noise; bench-only).

Strategic implication: at BIEN's scale and signal regime, further
RMSE lift needs environmental covariates (WorldClim) or larger
posterior tree ensembles, not GNN-architecture changes.

# pigauto 0.9.2 (2026-05-17)

First release after the v0.9.1.9000 dev cycle. Headline outcome:
on the BACE-compatible cross-dataset head-to-head bench (six
phylogenetic-trait datasets, n=2000 per dataset, 30% MCAR,
matched against `daniel1noble/BACE`'s snapshot results),
pigauto's net win-loss vs BACE went from 14–14 (commit 9d4782e,
PR #102) to **roughly 24–7** in this release. The change is
driven by three bug fixes in the in-house Σ solver plus a CI
evaluator alignment with BACE's reporting policy; the underlying
imputation algorithm is unchanged from v0.9.1.

## Bug fixes: in-house Σ solver

`fit_mvn_bm_inhouse()` had three independent bugs that surfaced
together on real phylogenetic data:

1. **Henderson R⁻¹ scale mismatch.** `henderson_R_inv_apply()` was
   documented to return R⁻¹b but actually returned A⁻¹b (where
   A = `vcv(tree)`, R = `cov2cor(vcv(tree))`). For variable-depth
   trees the mismatch is a per-tip rescale; for ultrametric trees
   it is a global scalar. Fixed by adding `cor_scale = TRUE`
   option that pre-/post-multiplies by `sqrt(diag(A))`.
   `henderson_R_inv_apply(..., cor_scale = TRUE)` now matches
   dense `solve(R)` to ~1e-6 relative on non-ultrametric trees;
   the existing `cor_scale = FALSE` (A-scale) tests continue to
   pass unchanged.

2. **Σ init via species-iid sample covariance.** The previous init
   used `cov(L, use = "pairwise.complete.obs")` which treats
   species as iid and double-counts sister tips as independent
   observations. On strongly-phylo-conserved liability matrices
   this produced near-singular Σ̂ initialisations whose off-diagonals
   were biased toward zero. Replaced with the closed-form
   matrix-normal MLE Σ̂ = (1/n) L̂ᵀR⁻¹L̂ on per-column-BM-completed
   L̂ — properly credits R for the species correlation and recovers
   Σ̂ to the right order of magnitude.

3. **EM E-step diverged at near-sister tips.** `build_conditional_prior()`
   computes a within-row cross-trait conditional posterior but ignores
   R-mediated cross-row correlation. For near-sister tips, the proper
   joint posterior would strongly pull a missing tip's value toward its
   observed sister; the within-row approximation cannot see that, and
   combined with `henderson_bm_predict()`'s constant conservative
   variance for missing cells, the inverse-variance pool over-weights
   the wrong cross-trait prior. R⁻¹ at near-sister pairs (values of
   4000–7500 observed at sister pairs on real AVONET) then amplifies
   the discrepancy and Σ̂ inflated 10–50× per EM iteration. Fixed by
   changing the default `max_iter` to `0L` so the solver returns the
   per-column BM L̂ + closed-form Σ̂ MLE without the broken cross-trait
   EM. Callers can still opt in to the EM via `max_iter > 0`.

Net effect on bench: AVONET categorical accuracy lifted from 0.357
to 0.825 (trophic_level vs BACE 0.816) and 0.350 to 0.823
(primary_lifestyle vs BACE 0.828). PanTHERIA discrete-trait
accuracies are now competitive with BACE.

## Improvements: CI head-to-head evaluator

Major changes to `script/gha/_bace_compat.R` and
`script/gha/make_headtohead_report.R` so pigauto's CI numbers
are directly comparable to BACE's snapshot.

- **Pool every trait via `mi$pooled_point`** to match BACE's
  reporting policy. BACE uses majority-vote-across-M-imputations
  for categorical / ordinal accuracy and `rowMeans(imp_mat)` for
  continuous RMSE. pigauto's `mi$pooled_point` carries the
  analogous central tendency (argmax-of-average-probability for
  categorical / binary, mode-or-mean for ordinal, median-or-mean
  for continuous depending on log-transform). Before this fix,
  `.bace_eval_per_imputation` scored every stochastic categorical
  sample individually — bounded by `E[p_max]` across cells —
  badly underestimating the model's accuracy. The previous
  per-draw default also inflated continuous RMSE by Jensen's
  inequality. A new `pool_method = "auto"` argument is the default;
  legacy `"per_draw"` behaviour stays available for diagnostics.
- **Added `coverage_95`, `interval_width`, and `brier`** to the
  canonical CI schema. Coverage and interval width are
  pigauto-only (BACE has no posterior-interval snapshot column).
  Brier is now computed on both sides for discrete traits so the
  two methods can be compared on calibration.
- **h2h ordinal merge bug fixed.** The previous report treated
  `ordinal` as a continuous-family type and tried to pull RMSE
  from both sides, returning NA for every ordinal trait. Removed
  ordinal from `CONTINUOUS_TYPES` so the report uses accuracy
  for ordinal too. AVONET migration, PanTHERIA diet_breadth /
  habitat_breadth now show up properly.
- **Robust h2h aggregator.** `load_pigauto_results()` and
  `load_bace_snapshot()` now route every loaded part through
  `.normalize_eval()` so rbind across schema-version mismatches
  no longer crashes the aggregator. The BACE snapshot, pinned to
  a pre-Tier-1.5 schema, is automatically padded with NA in the
  new pigauto-only columns.

## Improvements: CI defaults

- Default `n_imputations` in `PIGAUTO_CI_CONFIG` raised from 10
  to 20. Mild improvement on BIEN total RMSE (8.312 → 8.288 in
  an autoresearch sweep of 17 settings; m=30/50/100 all
  regressed, so 20 is the empirical sweet spot at the current
  CI scale).

## Reverts in dev cycle

- Tier 2 (AVONET covariates via the `cov_encoder` /
  `obs_refine` paths) was tried and reverted: CI runs showed
  AVONET categorical accuracy dropped sharply when covariates
  were passed (trophic_level 0.825 → 0.731, primary_lifestyle
  0.828 → 0.657) without lifting continuous traits. The
  underlying covariate infrastructure stays in place — when
  `covariate_cols` is empty the path is a no-op. Re-enabling
  for AVONET (or BIEN) needs a categorical-safe covariate
  routing first; planned for a future minor release.

## Tests

- New `tests/testthat/test-bace-compat-eval.R` covers the new
  pool semantics + coverage / brier computation paths.
- Henderson tests extended with explicit `cor_scale = TRUE`
  assertions on non-ultrametric trees and per-cell variance
  checks for `henderson_bm_predict()`.
- Full suite: 1511 pass / 0 fail (299 pre-existing small-n
  warnings unchanged).


# pigauto 0.9.1.9014 (dev)

## Bug fixes: prediction-time observed-cell context (2026-05-11)

- `predict.pigauto_fit()` now seeds DAE refinement with stored observed latent cells and uses the phylogenetic baseline only for originally missing cells. Previously prediction initialized every cell from the baseline, so the GNN did not receive the observed trait context it was trained to use. Evaluation, cross-validation, and report metrics explicitly mask validation/test cells from that context so held-out scores remain leakage-free.
- The direct `cov_linear` fixed-effect path is now included in training, validation, gate calibration, conformal scoring, and prediction whenever user covariates are present. Previously that path was constructed and returned by the model, but production callers ignored it unless using the optional low-level `baseline_mu` forward path.
- Conformal residuals are now computed from the same calibrated three-way safety-floor blend used at prediction time, rather than reconstructing a legacy two-way BM/GNN blend from `r_cal_gnn` alone.

## Improvements: tree-aware downstream fitting (2026-05-11)

- `with_imputations()` now carries tree metadata for `multi_impute_trees()` results. If the fitting callback declares `tree`, `tree_index`, or `imputation` arguments, they are filled with the matching posterior tree, its index, and the imputation number; the returned fits and `pool_mi()` output also retain the successful `tree_index` vector as metadata.

## Bug fixes: per-tree tree MI pooling (2026-05-11)

- `multi_impute_trees(share_gnn = FALSE)` now computes `pooled_point` by averaging all completed `T * m_per_tree` datasets, matching the documented return contract. Previously this path averaged one per-tree pooled result, which could ignore within-tree draw-level variation when `m_per_tree > 1`.

## Bug fixes: shared-GNN tree MI pooling (2026-05-11)

- `multi_impute_trees(share_gnn = TRUE)` now computes `pooled_point` from the completed imputation datasets, preserving `species_col`, user input row order, observed cells, and all `m_per_tree` draws. Previously this path pooled directly from raw prediction output, which could drop non-trait columns and return the point estimate in internal prediction order.

## Bug fixes: shared-GNN tree baselines (2026-05-11)

- `multi_impute_trees(share_gnn = TRUE)` now recomputes each posterior
  tree's baseline from that tree's own cophenetic distances and
  phylogenetic covariance. Previously the shared-GNN path accidentally
  reused the reference tree's cached graph when fitting per-tree baselines,
  which could collapse the intended tree-uncertainty signal.

## Bug fixes: proportion conformal intervals (2026-05-11)

- `predict.pigauto_fit()` now returns `conformal_lower` and `conformal_upper`
  columns for explicit `proportion` traits when conformal scores are
  available.  Cross-validation/report coverage and prediction interval plots
  now treat proportions as numeric back-transformable traits on the original
  0-1 scale.

## Bug fixes: multi-proportion conformal MI (2026-05-11)

- `multi_impute(draws_method = "conformal")` now generates stochastic
  draws for `multi_proportion` groups by sampling their CLR latent columns
  with BM latent uncertainty and decoding back to the simplex. Previously
  the default conformal MI path skipped these groups because their decoded
  output columns are the component names rather than the group name.

## Documentation: pkgdown reorganisation Phase 1 (2026-05-09)

User-facing docs reorganisation following the issues by @b1805 (#67, #68)
that surfaced confusion the README and `?impute` should have prevented.
PR #69 fixed the function-level docstring; this release reorganises the
pkgdown site itself:

- The Articles dropdown is slimmed from 39 entries to 4 source-backed
  vignettes. Older static walk-through HTML pages are hidden and marked
  archived until they are refreshed against the current README, vignettes,
  and benchmark evidence; the page-by-page disposition is recorded in
  `useful/pkgdown_page_rethink.md`.
- A new top-level Methodology dropdown owns benchmark HTMLs in two
  sub-sections: per-trait benches (7 working items) and cross-dataset
  benches + simulations (6 working items).  Ten older bench HTMLs that
  are not present in `docs/dev/` on the current build are documented as
  pending regeneration in the YAML and will re-enter the navbar once
  their `script/make_*_html.R` driver has been re-run; the Phase-2
  math vignette will land in a separate PR alongside this list refresh.
- New `vignettes/common-pitfalls.Rmd`: five sections covering
  `prediction$imputed` semantics, K=3 ordinal majority-class collapse,
  the safety-floor closed-gate behaviour, phylogenetic signal diagnosis,
  and `clamp_outliers` tail safety.  Each follows a fixed
  Symptom / Why / Diagnose / Fix / See-also template.
- `DOCS.md` retired; README + live pkgdown site are the single source of
  truth for documentation navigation.  README gains a one-line
  "Live documentation" link near the top.

No benchmark reruns or test changes.

## Hardening: `preprocess_traits()` errors on edge-case inputs (2026-05-03)

Two silent-degradation paths surfaced by the pre-shipping coverage
audit are now loud errors instead of "warn + degraded result":

- **Zero overlap** between `rownames(traits)` and `tree$tip.label`
  used to drop every trait row, emit a "tree tips have no trait
  data" message, and return a `pigauto_data` with all-NA `X_scaled`.
  Downstream `fit_pigauto()` would then fail in a less-obvious place.
  Now `preprocess_traits()` errors immediately with `"No overlap
  between trait row names and tree tip labels: 0 species in
  common."`.
- **Single-tip tree** used to be accepted and produce a 1-row
  `pigauto_data` with no phylogenetic signal.  Now errors with
  `"`tree` has fewer than 2 tips. Phylogenetic imputation requires
  at least 2 species."`.

No effect on any non-pathological input.

# pigauto 0.9.1.9013 (dev)

## Bug fixes (pre-shipping coverage audit, 2026-05-03)

Two production bugs surfaced by the pre-shipping coverage audit and
fixed in PR #65:

- `suggest_next_observation(by = "species")` no longer dies with
  `"no rows to aggregate"` on continuous-only traits.  The species
  aggregator now uses `stats::aggregate()`'s data.frame interface
  (which tolerates all-NA columns) instead of the formula interface
  (which routes through `model.frame()` and `na.omit`s every row
  when the lhs is all NA).  See `R/active_impute.R::692`.
- `pool_imputations()` for `count` and `zi_count` traits now
  honours `pool_method = "mean"` instead of silently falling back
  to median.  A leftover-from-refactor `vals <- row_median_decoded(nm)`
  line had been overriding the pool_method-aware assignment.  The
  `zi_count` branch additionally had a duplicate (and therefore
  unreachable) `else if (tm$type == "zi_count")` that has now been
  collapsed into a single branch which (a) honours pool_method on
  the magnitude column and (b) always probability-pools the gate
  column for `P(non-zero)`.  See `R/predict_pigauto.R::809+865`.

## Pre-shipping coverage tests (2026-05-03)

`tests/testthat/test-shipping-coverage.R` (+17 test_that blocks,
+68 expectations) covers:

- 9 exported functions that previously had zero direct test usage:
  `confusion_matrix`, `calibration_df`, `simulate_non_bm`,
  `read_tree`, `plot_history_gg`, `plot_uncertainty`.
- 3 non-default option paths previously untested:
  `multi_impute(draws_method = "mc_dropout")`,
  `fit_pigauto(use_attention = FALSE)`,
  `fit_pigauto(conformal_method = "bootstrap")`.
- 5 edge cases: zero-overlap, no branch lengths, single-species
  tree, all-observed identity round-trip, all-NA column.

## Documentation: clarify PMM is for multi-imputation, NOT tail safety (2026-05-01)

The bench evidence collected after PR #61 (Phase G' PMM) and the
unsuccessful Phase G'' attempts (per-draw and post-pool conditional
PMM) showed that **PMM is not a tail-safety tool for pigauto**.
PMM in mice works because linear-regression predictions are unbiased
estimators of truth across the prediction range.  pigauto's GNN can
be systematically biased on rare-clade species (the AVONET
Casuarius case) AND simultaneously accurate on legitimate-but-out-
of-observed-range species (the AVONET seed-2032 case).  PMM cannot
distinguish these two cases by donor-distance alone, so it can hurt
RMSE just as easily as it can help.

**Where PMM genuinely helps**: in multi-imputation workflows
(`n_imputations > 1` + `pool_mi()`), the donor-based between-
imputation variance is well-calibrated to the trait's marginal
distribution, giving Rubin's rules honest standard errors even
when the GNN's MC-dropout noise alone underestimates uncertainty.

This release updates the `?impute` and `?predict.pigauto_fit`
documentation + this NEWS entry to reflect the corrected framing.
**No code changes**.

For tail safety on single-imputation point estimates, the
recommended tool remains `clamp_outliers = TRUE` (Phase G PR #60).

# pigauto 0.9.1.9012 (dev)

## Phase G': `match_observed = "pmm"` -- Predictive Mean Matching (2026-05-01)

Adds Predictive Mean Matching (PMM; Little 1988; Buuren &
Groothuis-Oudshoorn 2011 mice) as an opt-in mechanism for
multiple-imputation workflows.  PMM replaces each missing cell's
back-transformed prediction with an actual observed value drawn
from a donor pool ranked by predicted-value proximity.

\strong{Originally framed as "tail safety", but the bench result
showed PMM is actually for MI workflows.}  See the v0.9.1.9013
NEWS section above for the corrected framing.

### What changed

Two new arguments on `impute()` and `predict.pigauto_fit()`:

* **`match_observed`** (character, default `"none"`).  Set to
  `"pmm"` to enable PMM.
* **`pmm_K`** (integer, default `5L` -- mice convention).  Donor
  pool size.

When `match_observed = "pmm"`, for each missing cell of an at-risk
trait type the function:

1. Computes the predicted value via the existing BM + GNN +
   blend pipeline (with optional MC dropout).
2. Finds the `pmm_K` observed cells whose own predictions are
   closest to the missing cell's prediction (the "donor pool").
3. Samples one donor uniformly at random.
4. Returns that donor's actual observed value as the imputation.

By construction, **PMM-imputed values are always in the observed
data range**.  No extrapolation is possible.

### When PMM applies

PMM is type-restricted to **at-risk** trait types where the
back-transform amplifies tail overshoots:

* log-transformed continuous (e.g. body mass)
* count (`expm1` decode)
* zi_count magnitude column (same)
* proportion (mild benefit)

Discrete-class types (binary / categorical / ordinal /
multi_proportion) are no-op: their decoders already return values
from the observed level set by construction.  Untransformed
continuous: no-op (no expm1 risk).

### Default policy

**`match_observed = "none"` remains the default** for backward
compatibility (PR #60 added `clamp_outliers` opt-in; this PR adds
PMM opt-in alongside it).  Cross-dataset evidence in v0.9.2 may
justify flipping the default to `"pmm"` for at-risk types.

Recommended use:
  * `match_observed = "pmm"` for log-transformed continuous traits
    with `n_imputations > 1` (multi-imputation context, where
    each draw picks a different donor and Rubin's rules give
    properly calibrated SEs).
  * `clamp_outliers = TRUE, clamp_factor = 5` for users who want
    a deterministic point estimate without PMM stochasticity.

### Implementation

* `R/pmm.R` (NEW) -- `pmm_impute_one_trait()`, `apply_pmm_to_decoded()`,
  `recover_X_orig()`, `pmm_is_eligible()`.  RNG state is saved and
  restored around per-call `set.seed()` so PMM does not leak the
  user's RNG state.
* `R/predict_pigauto.R` -- `predict.pigauto_fit()` and
  `decode_from_latent()` extended with `match_observed` and `pmm_K`
  arguments.  PMM is applied between decode and pool, per draw.
* `R/impute.R` -- user-facing wrapper extended with same args + roxygen.
* `R/fit_pigauto.R` -- the fit object now retains `X_scaled` so
  `predict()` can recover original-units observed values for the
  donor pool.  Adds ~8 * n_obs * p_latent bytes (~600 KB at AVONET
  full).
* `tests/testthat/test-pmm.R` -- 5 layers of tests, **65
  assertions**, all pass:
    - L1: helper unit tests (K=1, K=5, empty cases, ties, NA preds, RNG hygiene)
    - L2: type-dispatch (eligible types yes; discrete + un-log no)
    - L3: statistical properties (no extrapolation, observed
      cells unchanged, donor membership)
    - L4: integration (default = none, PMM + clamp coexist,
      non-eligible no-op)
    - L5: edge cases (invalid K, deterministic K=1)
* `script/bench_phase_g_prime_pmm.R` -- AVONET 4 log-cont traits +
  synthetic heavy-tail bench across 3 configs (none / clamp / pmm).

### Bench

`script/bench_phase_g_prime_pmm.R` reports RMSE, max prediction, and
the fraction of imputations falling in the observed range, for
each (dataset, seed, trait, config) cell.  See
`useful/MEMO_2026-05-01_phase_g_prime_results.md` for the verdict.

# pigauto 0.9.1.9011 (dev)

## Phase G: opt-in `clamp_outliers` for back-transform tail extrapolation (2026-05-01)

Closes the AVONET Mass tail-extrapolation mode documented in
`useful/MEMO_2026-05-01_avonet_mass_diag.md`: at N_IMP=20 on seed
2030, Casuarius bennetti was predicted at 538 kg vs truth 35 kg
because all 20 MC-dropout draws produced a +3-4 SD latent overshoot
that `expm1()` amplified to 250-1700 kg.

### What changed

Two new arguments on `impute()` and `predict.pigauto_fit()`:

* **`clamp_outliers`** (logical, default `FALSE`).  When `TRUE`,
  post-back-transform predictions for log-transformed continuous,
  count, and zi_count magnitude traits are capped at
  `tm$obs_max * clamp_factor` (where `tm$obs_max` is the observed
  maximum on the original scale, recorded at preprocess time).
* **`clamp_factor`** (numeric, default `5`).  Multiplicative factor
  on the observed maximum.  Anything > `obs_max * clamp_factor` is
  considered implausible and capped.

For non-at-risk types (untransformed continuous, ordinal, binary,
categorical, proportion, multi_proportion) the clamp is a no-op.

### Default policy

**`clamp_outliers = FALSE`** remains the default for backward
compatibility.  Users imputing log-transformed traits with
`n_imputations > 1` and concerned about tail blow-ups should pass
`clamp_outliers = TRUE` explicitly.  See
`useful/MEMO_2026-05-01_phase_g_results.md` for the AVONET
acceptance bench:

| seed | clamp OFF | clamp ON | reduction |
|------|-----------|----------|-----------|
| 2030 | 24,330 | 6,273  | -74 % |
| 2031 |    312 |   318  | +1.7 % (MC noise) |
| 2032 |    591 |   433  | -27 % |

Seed-2030 Casuarius prediction: 551,050 g (OFF) -> 167,847 g (ON,
capped at 35,000 * 5 = 175,000).

### Implementation

* `R/preprocess_traits.R`: records `entry$obs_max` and
  `entry$obs_min` (original-scale) for log-transformed continuous,
  count, and zi_count types.  Untransformed continuous: NA (not at
  risk of expm1 amplification).
* `R/predict_pigauto.R`: extends `predict.pigauto_fit()` and
  `decode_from_latent()` with `clamp_outliers` / `clamp_factor`
  arguments.  New `.clamp_high()` helper applies the cap at the
  three back-transform sites (continuous-log, count, zi_count
  magnitude).
* `R/impute.R`: extends the user-facing wrapper with the same
  arguments + roxygen.

### Tests

`tests/testthat/test-clamp-outliers.R`: 5 test_that blocks
(24 assertions), all pass.  Covers obs_max / obs_min recording,
default-OFF no-op, ON behaviour on synthetic tail, clamp_factor
scaling, and input validation.

### Bench

`script/bench_phase_g_clamp.R` runs the AVONET seed-2030/2031/2032
N_IMP=20 sweep with clamp ON and OFF.  Output:
`script/bench_phase_g_clamp.{rds,md}`.

# pigauto 0.9.1.9010 (dev)

## Phase H: `pool_method = "mode"` for ordinal traits (2026-05-01)

Closes the AVONET Migration N_IMP-dependent regression flagged in
`useful/MEMO_2026-05-01_phase_f_smoke_results.md`.

### What changed

`predict.pigauto_fit()` and `impute()` now accept a third value for
`pool_method`:

* `"median"` (default, unchanged) — per-cell median for count /
  proportion / zi_count magnitude; per-cell `mean(integer-class) ->
  round` for ordinal.
* `"mean"` (unchanged) — pre-v0.9.2 arithmetic-mean pooling.
* **`"mode"` (new)** — per-cell **majority vote** across the M
  draws for ordinal traits.  For continuous-family traits ("mode"
  falls back to "median" since mode does not apply to continuous
  decoders).  Binary, categorical, and multi_proportion are
  unaffected (already probability-averaged + argmax).

### Why this matters

The previous ordinal pooling (`mean of integer class indices, then
round`) biases predictions toward middle classes when MC-dropout
spreads the M draws across adjacent classes.  Concrete K=3 failure:

  True class 1 (sedentary), 5 dropout draws produce {1,1,1,2,3}.
    Mean = 1.6, round = 2 -> "partial migrant" (wrong).
    Mode = 1 -> "sedentary" (correct).

Empirical impact on AVONET Migration (3 seeds, n=1500, miss=0.30):

| pool_method | N_IMP=1     | N_IMP=20    |
|-------------|-------------|-------------|
| median (default) | 0.767 +/- 0.011 | **0.713 +/- 0.032** |
| **mode (new)**   | 0.767 +/- 0.011 | **0.779 +/- 0.011** |

Mean baseline (per-class mode): 0.800 +/- 0.013.

Switching to `"mode"` improves Migration accuracy at N_IMP=20 by
**+6.6 pp** and reduces the regression vs baseline from -10.9 pp
to -2.1 pp.  Pre-registered acceptance criterion (Migration acc at
N_IMP=20 >= N_IMP=1) passes 3 / 3 seeds.

### Default policy

**`"median"` remains the default** for backward compatibility.
Users running ordinal traits with `n_imputations > 1` should pass
`pool_method = "mode"` explicitly until cross-dataset evidence
justifies a default flip.  See
`useful/MEMO_2026-05-01_phase_h_results.md`.

### Implementation

* `R/predict_pigauto.R`: extended pool_method enum; new ordinal
  branch using `tabulate(class_indices, nbins = K) -> which.max()`.
  Existing median-vs-mean guards on count / proportion / zi_count
  branches extended to treat "mode" as "median" (fall back).
* `R/impute.R`: extended user-facing `pool_method` enum + roxygen.
* `tests/testthat/test-mi-pool-robust.R`: 5 new assertions covering
  the mode-pool ordinal path and the median-fallback for continuous.

### Bench

`script/bench_phase_h_pool.R` runs the 3-seed AVONET-1500 Migration
sweep at `N_IMP in {1, 20}` with `pool_method = "mode"`. Output:
`script/bench_phase_h_pool.{rds,md}` and the verdict memo at
`useful/MEMO_2026-05-01_phase_h_results.md`.

# pigauto 0.9.1.9009 (dev)

## `suggest_next_observation()` v2: zi_count + multi_proportion (2026-05-01)

Closes the v1.5 scope gap: `suggest_next_observation()` now supports
all 8 pigauto trait types.

### zi_count

Hybrid metric: observing a missing zi_count cell reveals the gate
value (entropy reduction at the gate column, computed via the LP
binary formula) AND, with probability \eqn{p_{\mathrm{gate}}}
(current LP estimate at \eqn{s_{\mathrm{new}}}), reveals a
magnitude (variance reduction at the magnitude column, computed
via the BM Sherman-Morrison formula on the gate=1 subset).  Output
rows for zi_count populate BOTH `delta_var_total` (= probability-
weighted magnitude variance reduction) AND `delta_entropy_total`
(= gate entropy reduction).  `metric` is set to `"variance"` so
the row sorts on the magnitude scale; `delta_entropy_total` is
available for users who care about gate-uncertainty separately.

### multi_proportion

Per-component variance reduction summed across K CLR-z-scored
latent columns.  Per CLAUDE.md, multi_proportion rows are
fully observed OR fully missing (no partial K-component
observations), so `delta_var_total` is the K-component sum at each
fully-missing row.

### Implementation

`R/active_impute.R`:
- New dispatch branches for `zi_count` and `multi_proportion`
  in `suggest_next_observation()`.
- Reuses existing `bm_variance_reduction()` and
  `lp_entropy_reduction_binary()` helpers.
- New `needs_R_phy` / `needs_sim_phylo` derived flags clarify
  which kernels each trait type needs (zi_count needs both;
  multi_proportion needs R_phy only).

### Tests

4 new test_that blocks (`tests/testthat/test-active-impute.R`):

* `[active v2] zi_count: hybrid variance + entropy populated`
* `[active v2] zi_count: types argument can exclude zi_count`
* `[active v2] multi_proportion: per-component BM variance summed`
* `[active v2] multi_proportion: delta equals K-component sum`
  (verifies via independent recomputation against
  `bm_variance_reduction()` per component)

Total active-impute tests: 17 / 60 expects, all pass.

# pigauto 0.9.1.9008 (dev)

## NEW: `suggest_next_observation()` -- active-imputation guidance (2026-04-30)

For a fitted `pigauto_result` (from `impute()`), compute the
closed-form expected reduction in TOTAL predictive variance across
all currently-missing cells if each candidate cell were observed
next.  Useful for sampling-design guidance when measurement budget
is limited: which species, if observed next, would most reduce
imputation uncertainty across the rest of your tree?

### API

```r
suggest_next_observation(result, top_n = 10L,
                          by = c("cell", "species"))
```

Returns a `pigauto_active` data.frame (descending by
`delta_var_total`).  Columns:

* `by = "cell"`: `species`, `trait`, `type`, `delta_var_total`.
* `by = "species"`: `species`, `delta_var_total`,
  `n_traits_missing`.

### Math

Built on a Sherman-Morrison rank-1 inverse update of the BM
conditional MVN: adding species `s_new` to the observed set updates
`R_oo^{-1}` by a known closed form, and the variance reduction at
every other missing cell is a single scalar per candidate.

Closed form per candidate cell:

```
delta_V_total(s_new) = sigma2 * sum_{i in miss} D[i, k]^2 / alpha_k
```

where `D = R_mm - R_mo R_oo^{-1} R_om` is the residual matrix and
`alpha_k = D[k, k]` is the current relative leverage of cell `k`.

The formula assumes `sigma2` (REML variance) is held fixed -- the
standard active-learning convention.  Variance reduction depends
only on the geometry of which cells are observed, not on the
unknown value at `s_new`.

### Scope (v1.5, 2026-04-30 evening)

* Continuous-family traits (`continuous`, `count`, `ordinal`,
  `proportion`) -- variance reduction via the BM Sherman-Morrison
  closed form above.
* Binary + categorical traits (added 2026-04-30 evening) --
  expected entropy reduction via the label-propagation closed form:
  current LP probability \eqn{p_i} has entropy \eqn{H(p_i)};
  observing \eqn{s_{\mathrm{new}}} updates LP linearly via the
  similarity kernel; expected entropy after observation is averaged
  over \eqn{P(y_{\mathrm{new}})} = current LP estimate at
  \eqn{s_{\mathrm{new}}}.
* `zi_count` and `multi_proportion` not yet supported (silently
  skipped).
* Single-obs only.  Multi-obs deferred to v2.

**Variance and entropy are not directly comparable.**  The
default output mixes them by `delta` value, which is approximate
across metrics.  When precise ranking matters, filter by
`metric` ("variance" or "entropy") first.

### Why this is novel

To my knowledge, no other phylogenetic-imputation package (Rphylopars,
BACE, phylolm, mice with phylogenetic correlation) exposes a
sampling-design helper.  pigauto is uniquely positioned because its
BM conditional MVN already computes per-cell predictive variance --
this function just exposes the closed-form variance-reduction
formula derived from rank-1 update.  The methods are 30 years old
in optimal-design literature (Cohn et al. 1996, JAIR 4:129–145);
applying them to phylogenetic trait imputation appears to be new.

### Files

* `R/active_impute.R` (new) -- 304 LOC, hosts
  `bm_variance_reduction()` (internal) and
  `suggest_next_observation()` (exported), plus
  `print.pigauto_active` S3 method.
* `R/impute.R` -- adds `tree` to the `pigauto_result` so
  `suggest_next_observation()` can recompute `R_phy` if it was
  freed in the post-fit cleanup block.
* `tests/testthat/test-active-impute.R` (new) -- 10 tests, 29
  expects.  Verifies closed-form against an independent fixed-
  sigma2 brute-force refit on a 20-tip fixture; descending order;
  by-species aggregation; multi-obs error path; discrete-skip.

### Smoke

On AVONET 300 with 90 cells masked across 4 continuous traits:
top-1 cell is `Chloephaga_rubidiceps Beak.Length_Culmen`,
delta_var_total = 2.99; top-1 species is `Chloephaga_rubidiceps`
with 2 currently-missing traits combining to delta = 3.23.

## NEW: `gate_method = "cv_folds"` -- cross-validated gate calibration (2026-04-30)

Closes the longest-standing open NEWS item: *"A future improvement
(cross-validated gate selection, ...) could close this further; not
in scope for this release."*

### What changed

`calibrate_gates()` and `fit_pigauto()` gain `gate_method =
"cv_folds"` (opt-in alongside the existing `"single_split"` default
and `"median_splits"` option).  When selected, val cells are
partitioned into `gate_cv_folds = 5L` (default) deterministic
non-overlapping folds; the existing grid + half-B-verify procedure
runs once per fold (`half_a` = K-1 training folds, `half_b` =
held-out fold), and the componentwise median of K winning weight
vectors becomes `w_final`.  Post-cal full-val invariant fires
unchanged.

### Compared to `"median_splits"`

| Property | median_splits | cv_folds |
|---|---|---|
| half_a per split | 1/2 of n_val | (K-1)/K of n_val |
| Number of splits | B=31 random | K=5 deterministic |
| Overlap | Random (each cell in many splits) | None |
| Statistical interpretation | Median of B noisy gate selections | Standard k-fold CV |

Larger training sets per split → less noisy w_k.  Standard CV gives
a cleaner story for users and reviewers.

### Motivation

The strict val-floor enforces `pigauto_val <= baseline_val`, but
val→test drift survives in 4/32 binary cells (per
`useful/MEMO_2026-04-29_discrete_bench_reruns.md`).  Per-trait gate
calibration overfits to the specific val cells because the same
cells score every candidate weight.  Larger training sets per
split reduce overfitting; non-overlapping folds give a true held-
out evaluation.

### Edge cases

* `gate_cv_folds < 2` errors with a clear message.
* `n_val < gate_cv_folds` caps `K_eff = n_val` so each fold has at
  least 1 cell.
* `n_val < 2` (CV ill-defined) falls back to a single split.
* Existing back-compat preserved: default `gate_method =
  "single_split"`, and the `"median_splits"` path is bit-identical
  to pre-fix.

### Tests

5 new tests in `tests/testthat/test-safety-floor.R` (under
"# CV-fold gate calibration"):

* Valid simplex output
* Deterministic for fixed seed
* Pure-MEAN snap when MEAN dominates (post-cal invariant fires)
* `n_val < K` graceful fallback
* `gate_cv_folds < 2` error

All pass.  Existing 24 safety-floor tests + B4 edge cases unchanged.

### Smoke (AVONET 300, 60 epochs, 1 imputation, 60 Mass + 60 Migration cells masked)

| gate_method | wall (s) | r_cal_gnn[col_3] |
|---|---|---|
| single_split | 82.6 | 0.000 |
| median_splits | 77.2 | 0.050 |
| cv_folds | 75.9 | 0.150 |

`cv_folds` is slightly faster (smaller per-split work) and selects
a non-zero GNN gate on the col where the others snap to BM.

### Bench (focused 3-way comparison on synthetic mixed-type)

`script/bench_cv_vs_median.R` runs a focused 3-way comparison on
a synthetic continuous + binary + categorical fixture (BM, lambda=1,
n=300, 30 % mask, 5 reps).  Test-set means across reps:

| Method | bin (acc) | cat3 (acc) | cont (RMSE) |
|---|---|---|---|
| single_split  | 0.8222 | 0.7356 | 1.0778 |
| median_splits | 0.8222 | 0.7356 | 1.0720 |
| **cv_folds**  | 0.8222 | 0.7356 | **1.0611** |

* **Continuous**: cv_folds reduces RMSE by 1.6 % vs single_split
  and by 1.0 % vs median_splits.  Modest but consistent across
  the 5 reps (per-rep cv_folds < median_splits in 4/5 reps).
* **Discrete**: identical across all three methods on this
  fixture.  At lambda=1 the strict val-floor snaps the gate to
  the pure-BM corner regardless of how the gate was selected,
  so the three methods produce identical class assignments.
  This is consistent with what the discrete-bench memo
  (`useful/MEMO_2026-04-29_discrete_bench_reruns.md`) showed:
  the val→test drift on 4/32 binary cells is val→test
  extrapolation noise, not gate-selection noise; cv_folds does
  not close it on this regime.  The hypothesis that larger
  training sets per split would close MORE of the val→test
  drift than median_splits is **not supported** by this
  fixture.  cv_folds is preferable for continuous traits but
  not transformative for discrete.

Full per-type bench reruns (`bench_binary.R`, `bench_categorical.R`,
`bench_zi_count.R` at multiple signal levels) deferred to a
later session.

## Strict val-floor v3 — type-conditional safety floor (2026-04-29)

A type-restriction bug in `R/fit_helpers.R::calibrate_gates` was
silently exempting binary, categorical, and zi_count traits from the
post-calibration "strict val-floor" invariant.  The buggy behaviour
let pigauto underperform its own baseline on the package's own
simulator benches (`bench_binary.md`, `bench_categorical.md`,
`bench_zi_count.md`):

* binary `signal_0.6` mean accuracy: −4.8 pp below baseline
* categorical `K_3` cat2: −12 pp below baseline
* zi_count `zf_0.2` zi1: +23 % RMSE above baseline

These violate the safety-floor's "pigauto cannot underperform its
own baseline" promise.  The bug was introduced 2026-04-23 with the
initial safety-floor commit (`28d8e45`); the type filter mirrored an
earlier continuous-only MSE evaluator and was not extended when the
simplex grid landed.

The fix splits the post-calibration invariant into two
type-specific blocks:

* **Discrete family** (binary, categorical, zi_count) gets a strict
  full-val-set check against both pure-BM and pure-MEAN corners; if
  the calibrated blend's val-loss exceeds either corner by 1e-12,
  the gate is overridden to whichever pure corner is better.
  Discrete losses are 0-1 / argmax step functions, so this strict
  check is well-conditioned at typical val sizes.
* **Continuous family** (continuous, count, ordinal, proportion)
  keeps the pre-2026-04-29 mean-only check verbatim (override only
  if the blend loses to pure-MEAN, never to pure-BM).  An
  intermediate fix that also added a strict pure-BM check for
  continuous types over-corrected on high-phylo-signal traits where
  the blend's val MSE is numerically equal to pure-BM's val MSE
  within sampling noise (AVONET Mass test RMSE went 386 → 659 in
  the over-strict variant); v3 reverts continuous to the original
  conservative check.
* **multi_proportion** is skipped (gate naturally closes via the
  legacy path; bench shows pigauto = baseline at every signal
  level).

Empirical evidence (CORRECTED 2026-04-30 per the Opus adversarial
review; the original tabulation here mistakenly averaged val + test
splits in the bench RDS, overstating the closure).  The user-facing
held-out **test-set** mean acc gap pigauto − baseline (lower-magnitude
= closer to baseline-parity floor):

| bench | Apr 17 (pre-bug) | Apr 28 (buggy) | Apr 29 (post-fix v3) |
|---|---|---|---|
| binary | −0.029 | −0.027 | **−0.014** |
| categorical | −0.031 | −0.026 | **−0.016** |
| zi_count (RMSE ratio) | 1.058 | 1.043 | **1.012** |

The strict val-floor closes ~50 % of the test-side regression on
binary and categorical (from −2.7 / −2.6 pp to −1.4 / −1.6 pp); the
remaining ~1.4 pp is val→test sampling drift.  zi_count test ratio
drops from 1.06 to 1.01 -- most of the regression closes, with a
small residual.  3 of 32 binary cells still drift ≥ 5 pp on test
(val→test extrapolation that no val-only floor can prevent).

The strict val-floor's invariant is enforced at **calibration time**
on the calibration-time loss surface
(`compute_corner_loss(g, val_row_idx, ...)`).  At predict time the
model runs additional refine_steps and mask_token substitution that
introduce a small (~5 %) drift; the
`tests/testthat/test-safety-floor.R` strict-floor invariant test
documents this with a `* 1.05` slack.  So the user-facing guarantee
is "pigauto val-loss ≤ baseline val-loss + 1e-12 at calibration time,
modulo ~5 % refine_step drift at predict time" -- not the bit-exact
"≤ 1e-12 everywhere" reading.

Continuous benchmarks (AVONET n=1500 single seed × N_IMP=5):
the originally-reported "v3 within MC-dropout noise of pre-fix"
landed on a particularly close pair of draws.  Re-running the same
seed × n × N_IMP three times produced Mass test RMSE = 377, 430, 504
-- the comparison is non-deterministic at this evidence level
(pooling order over the 5 MC dropout draws varies between R sessions).
A multi-seed × N_IMP=20 verification is queued; the single-seed
single-run evidence does NOT support strong claims about whether
v3 changed continuous performance, only that the change (if any)
is within run-to-run pooling noise on this dataset.

Per-cell variance survives in multiple binary cells (val→test
sampling drift).  The strict val-floor only guarantees
`pigauto_val_loss ≤ baseline_val_loss + 1e-12` (calibration time);
it does not guarantee `pigauto_test_loss ≤ baseline_test_loss`.  A
future improvement (cross-validated gate selection, epsilon-margin
tightening, or LP added as a corner option for ordinal) could
close this further; not in scope for this release.

R CMD check after the fix: 0 errors / 0 warnings / 0 notes.
Full `devtools::test()`: FAIL 0 / SKIP 2 / PASS 1081.

New regression test in `tests/testthat/test-safety-floor.R`:
`"strict val-floor: pigauto val-loss <= baseline val-loss per
trait, all types"` — fits on a mixed continuous + binary +
categorical + count fixture, recomputes val-loss for the
calibrated blend and the pure baseline, asserts blend ≤ baseline ×
1.05 for every trait.

Calibration evidence preserved at
`script/_strict_floor_v3_calibration/` (4-way comparison: pre-fix /
v1 over-strict / v2 intermediate / v3 current).  Memos:
`useful/MEMO_2026-04-29_*.md`.

## Per-trait ordinal baseline path selection (Opus #6, 2026-04-30)

`fit_baseline()` now computes BOTH the threshold-joint baseline (B3,
commit `a541dbd`) AND a per-column BM-via-MVN alternative for each
ordinal trait, then picks the lower-val-MSE path per trait against
`splits$val_idx`.  The chosen path is exposed as
`fit_baseline()$ordinal_path_chosen` (named character: latent col →
`"threshold_joint"` or `"bm_mvn"`).

Motivation: the AVONET-300 Phase 6 bench's Migration RMSE regressed
from 0.879 (Apr 16) to 0.975 (Apr 29) when B3 swapped the ordinal
baseline path from a label-propagation-equivalent BM-via-MVN to the
truncated-Gaussian threshold-joint approach.  At K=3 ordinals the
K−1 thresholds are pinned to a narrow band by phylopars EM,
producing systematically worse predictions than BM-via-MVN on
z-scored integer class.  See
`useful/MEMO_2026-04-29_phase6_migration_bisect.md`.

Per Opus's adversarial review #6, we did **not** ship a `K ≤ 3 → LP`
heuristic.  Instead the per-trait selection runs at every K and
the gate's safety-floor simplex no longer has to route around
threshold-joint output for low-K ordinals -- the lower-val-MSE
baseline is selected before the gate calibration sees it.

Scope: single-obs only (multi-obs ordinal selection is out of scope
for this fix; it would require species-level aggregation of the BM
alternative).  The code path is gated by `!multi_obs && !is.null(splits)`
and only activates when threshold-joint actually populated the
ordinal column.

Two regression tests in
`tests/testthat/test-joint-threshold-baseline.R`:

* `"fit_baseline picks per-trait ordinal path and reports it"` --
  smoke test on K=3 ordinal + continuous mixed data; verifies finite
  output and that `$ordinal_path_chosen` is one of the two valid
  values.
* `"fit_baseline ordinal selection picks BM when threshold-joint
  loses on val"` -- recomputes both paths' val MSE independently
  and asserts the selection picked the lower one.

## Adversarial-review fixes (Opus 2026-04-28)

Four correctness fixes surfaced by the Opus adversarial pass on
top of PR #49 (the safety-floor / mean-gate work). All four are
silent-failure modes that produced plausible-looking but invalid
output; none of them threw a useful error before this release.

- **HIGH — `multi_impute(draws_method = "conformal")` no longer
  collapses to zero between-draw variance on shuffled input.**
  `R/multi_impute.R::.sample_conformal_draw` now threads
  `input_row_order` so it correctly maps user-input mask indices
  to internal (tree-tip) row indices when perturbing
  `pred$imputed`. Previously, on any data frame whose rows are
  not already in tree-tip order, the sampler perturbed the
  wrong cells and `build_completed` re-aligned the unperturbed
  values back into the user's missing cells, producing M
  identical datasets and `var = 0` after pooling. Any prior
  Rubin's-rules pooling on shuffled-input multi-imputation was
  unreliable. Caught by a new shuffled-input variance regression
  test; pre-fix `per_cell_sd ≡ 0`, post-fix `> 1e-8` at every
  masked cell.

- **MEDIUM (opt-in) — `conformal_split_val` argument on
  `fit_pigauto()` lets users avoid the validation-set
  double-dipping that breaks split-conformal exchangeability.**
  When `TRUE`, the validation cells are split per-column into a
  calibration half (used by `calibrate_gates()` to pick the
  blend weights) and a conformal half (used to estimate the
  residual quantile). Reusing the same val cells for both
  selects the gate to minimise residual MSE on the very cells
  whose residuals drive the conformal quantile, producing
  systematic undercoverage that's most visible at small `n_val`
  (matches the `coverage_investigation` memo).
  Default is `FALSE` (pre-fix single-set behaviour) because
  forcing the split regresses the AVONET300 / OVR-categorical /
  BIEN safety-floor smoke benches by 2-26% on small-val
  datasets. Use `conformal_split_val = TRUE` when accurate 95%
  coverage matters more than bench-grade RMSE.

- **MEDIUM — `calibrate_gates` half-A loop is now NA-safe.**
  Mirrors the half-B fix from commit 573decd: when
  `loss_a_pure_bm` is non-finite (e.g. multi_proportion CLR
  latents whose half-A subset is too small for some component),
  `best_la` is now seeded with `Inf` and the legacy
  pure-BM-relative-gain filter is skipped rather than throwing
  `missing value where TRUE/FALSE needed`.

- **LOW — `predict.pigauto_fit` no longer produces a cryptic
  torch shape error when calibrated_gates is overridden with a
  vector of the wrong length.** Length 0 is now treated as "no
  calibration" (degrades gracefully to the learned-gate path);
  any other non-matching length errors with an explicit
  `"calibrated_gates has length k but the model has p latent
  column(s)"` message. Unblocks the AE-attribution ablation
  workflow used by `script/phase1_gnn_ablation.R`.

Each fix has a dedicated regression test:
`tests/testthat/test-multi-impute.R::"multi_impute(draws_method=\"conformal\") produces non-zero between-draw variance on shuffled input"`,
`tests/testthat/test-fit-predict.R::"conformal_split_val toggle changes conformal scores (no double-dipping)"`,
and
`tests/testthat/test-fit-predict.R::"predict.pigauto_fit catches wrong-length calibrated_gates override"`.

# pigauto 0.9.1.9006 (dev)

## Per-occurrence WorldClim covariates (v1.1 follow-up to B.2)

Fixes the core limitation of PR #47 (centroid-only bioclim extraction):
bioclim is now extracted at every GBIF occurrence point per species
(not just the centroid), so the IQR columns carry real range-breadth
signal instead of being constant zero. This is what actually lets the
safety-floor calibrator see a meaningful covariate input and open the
GNN gate on plants.

Validated by the honest-sim diagnostic
(`experiment/covariate-honest-sim`) which showed pigauto's architecture
can produce 32 % RMSE lift on strong-environment traits when
covariates carry real signal — something centroid-only bioclim
failed to deliver on plants.

### API

- `pull_gbif_centroids(..., store_points = FALSE)` — new arg,
  default `FALSE` preserves pre-v0.9.1.9006 cache format. When
  `TRUE`, persists the raw filtered lat/lon points in each species'
  cache RDS, enabling per-occurrence bioclim extraction.
- `pull_worldclim_per_species()` — no signature change.  Internally,
  `.wc_extract_one()` now prefers cached per-occurrence points,
  falling back to centroid only when `points` field is absent
  (legacy caches).

### Back-compat

Pre-v1.1 cache RDS files remain fully readable. Users who want
per-occurrence extraction upgrade their cache with one call:

```r
pull_gbif_centroids(sp_list,
  cache_dir = "script/data-cache/gbif",
  store_points = TRUE,
  refresh_cache = TRUE)   # one-time refresh
```

After that, subsequent calls are instantaneous (cache hit).

### Fixture regeneration

`tests/testthat/fixtures/worldclim_plants_300.rds` re-generated with
per-occurrence aggregation. IQR columns are now informative.

### Follow-ups

- Covariate-aware phylo-signal gate — nice-to-have for
  interpretability (gate de-triggers when covariates resolve weak-
  phylo-signal traits). NOT a blocker given the safety-floor
  calibrator opens the gate on its own.
- Full-scale plants bench (n=4,745) — operator-run via
  `script/bench_bien_worldclim.R`. Expected: SLA r >= 0.35,
  leaf_area r >= 0.30.

# pigauto 0.9.1.9005 (dev)

## WorldClim bioclim covariates (B.2)

New exported helper `pull_worldclim_per_species(species,
gbif_cache_dir, worldclim_cache_dir, ...)` extracts 19 WorldClim v2.1
bioclim variables at each species' GBIF occurrence points, aggregates
per species (median + IQR), and returns a 38-column numeric
data.frame ready for `impute(..., covariates = ...)`.

### What it does

1. Resolves each species's GBIF occurrence cache from
   `pull_gbif_centroids()` (B.1).
2. Downloads WorldClim v2.1 10-arc-minute raster stack (~500 MB) to
   `worldclim_cache_dir` on first call; sentinel file makes subsequent
   calls offline.
3. For each species, extracts bioclim values at its centroid (v1;
   per-occurrence extraction deferred to v1.1) via `terra::extract()`.
4. Aggregates to median + IQR per variable (38 columns total).
5. Caches the per-species extract as one RDS in
   `worldclim_cache_dir/extracts/`.

### Important interaction with the phylogenetic-signal gate

The `phylo_signal_gate` introduced in v0.9.1.9003 fires on *raw trait
values*, so traits with weak phylogenetic signal (e.g. plants SLA,
leaf_area, height_m) get routed directly to the grand-mean corner
BEFORE the GNN sees any covariates.  **To see bioclim lift those
traits, users must set `phylo_signal_gate = FALSE`:**

```r
cov_bio <- pull_worldclim_per_species(
  rownames(traits),
  gbif_cache_dir = "script/data-cache/gbif",
  worldclim_cache_dir = "script/data-cache/worldclim")
res <- impute(traits, tree, covariates = cov_bio,
              phylo_signal_gate = FALSE)   # <-- required
```

A follow-up spec ("covariate-aware phylo-signal gating") will compute
lambda on residuals-after-covariates so the gate composes cleanly with
covariates without the manual opt-out.

### Dependency

`terra` added to Suggests.  `pull_worldclim_per_species()` errors
with an install hint when `terra` is absent.

### Testing

- `tests/testthat/test-worldclim-covariates.R`: 20 offline
  expectations (cache key, aggregation, download sentinel, extract
  cache, NA handling).
- `tests/testthat/fixtures/worldclim_plants_300.rds`: pre-extracted
  fixture for 300 BIEN species.
- End-to-end smoke (NOT_CRAN-gated): bioclim lifts SLA or leaf_area
  RMSE by >= 10% on a BIEN 200-subset when combined with
  `phylo_signal_gate = FALSE`.

### Follow-ups

- True per-occurrence extraction (v1.1): extend B.1's GBIF cache to
  store raw lat/lon point lists, then aggregate bioclim at each point
  rather than at the centroid.
- SoilGrids (B.3): add 5-10 soil variables similarly.
- Covariate-aware phylo-signal gating: recompute lambda on
  residuals-after-covariates so the gate de-triggers automatically
  when covariates resolve weak-phylo-signal traits.

# pigauto 0.9.1.9004 (dev)

## GBIF range-centroid covariates

New exported helper `pull_gbif_centroids(species, cache_dir, ...)`
fetches species occurrence records from the Global Biodiversity
Information Facility (GBIF) and returns a data.frame of per-species
median lat/lon centroids ready for `impute(..., covariates = ...)`.

### What it does

For each species, pull_gbif_centroids():

1. Resolves the taxon via rgbif::name_backbone() (handles synonyms).
2. Fetches up to occurrence_limit records via rgbif::occ_search(
   hasCoordinate = TRUE), paginated at 300/call.
3. Drops records with hasGeospatialIssues = TRUE,
   basisOfRecord in c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"),
   or out-of-range coordinates.
4. Computes median lat / lon per species (centroids more robust than
   means to outliers).
5. Caches one RDS per species at cache_dir (strongly recommended to
   avoid GBIF rate-limit issues).

### Why it matters

Plant traits like SLA and leaf area are environment-driven. Even the
simplest geographic covariate -- a species range centroid from GBIF --
carries biogeographic signal that phylogeny alone doesn't. This is the
B.1 proof-of-pipe for the pigauto covariate path; the full B.2 spec
(WorldClim bioclim + SoilGrids raster extraction) follows.

### End-to-end smoke (BIEN 200-subset)

Adding GBIF centroids as covariates keeps SLA and leaf_area RMSE
within +10% of the no-cov baseline (ratios ~1.09 and ~1.01 at n=200).
At this sample size the val-to-test sampling noise dominates any
covariate lift signal; production-grade positive lift reserved for
the B.2 bioclim run at larger n.

### API

- pull_gbif_centroids(species, cache_dir = NULL,
                        occurrence_limit = 500L, sleep_ms = 100L,
                        verbose = TRUE, refresh_cache = FALSE)
- Returns data.frame with species / centroid_lat / centroid_lon /
  n_occurrences, rownames = species.

### Dependency

rgbif added to Suggests (not Imports). pull_gbif_centroids() fails
with a helpful install hint when rgbif is absent.

### Testing

- tests/testthat/test-gbif-centroids.R: 32 offline expectations (cache
  key, aggregation, filtering, cache hit/miss, public function
  integration via mocks).
- End-to-end fixture at tests/testthat/fixtures/gbif_plants_300.rds
  (259/300 species with valid centroids) + live smoke test gated by
  NOT_CRAN=TRUE environment variable.

### Follow-up (separate spec)

B.2: full WorldClim bioclim + SoilGrids extraction at each GBIF
occurrence, aggregated per species. Expected to lift plants RMSE
significantly more than centroids alone.
# pigauto 0.9.1.9003 (dev)

## Phylogenetic-signal gate

Default change since v0.9.1.9003: `impute()` / `fit_pigauto()` now
compute per-trait Pagel's lambda on training-observed cells and route
traits with `lambda < 0.2` to the grand-mean corner `(0, 0, 1)` of the
safety-floor simplex, skipping BM fitting and GNN training effort on
those traits. Complements the v0.9.1.9002 safety floor with an
explicit per-trait decision: "pigauto knows when phylogeny helps".

### Effect on existing benches

- Birds / mammals / fish / amphibians: zero traits gated (all lambda >= 0.2).
  Results bit-identical to v0.9.1.9002.
- Plants (BIEN x V.PhyloMaker2 n=500 smoke): 3 of 5 traits flagged as
  weak-signal (leaf_area, sla, height_m) and gated to grand mean.
  fit$phylo_gate_triggered reports which.

### API

- impute(..., phylo_signal_gate = TRUE, phylo_signal_threshold = 0.2,
   phylo_signal_method = "lambda") - new args.
- fit_pigauto() same pass-through.
- Fit object gains four new slots: phylo_signal_per_trait (named
  numeric lambda values), phylo_gate_triggered (named logical),
  phylo_signal_method, phylo_signal_threshold.

### Opt-out

phylo_signal_gate = FALSE reproduces v0.9.1.9002 safety-floor-only
behaviour bit-identically.

### Reporting

print.pigauto_fit() now surfaces a "Phylogenetic signal" section
listing gated vs passed traits when any trait was gated.

### Dependencies

phytools added to Suggests. When absent, phylo_signal_gate degrades
gracefully to FALSE (safety-floor-only) with a one-time warning.

# pigauto 0.9.1.9002 (dev)

## Safety floor: pigauto is never worse than the grand mean

Default change: `impute()` / `fit_pigauto()` / `multi_impute*` now run
with `safety_floor = TRUE`. The calibration grid post-training searches
a three-way convex combination of (BM baseline, GNN delta, grand mean)
on a 231-point simplex at step 0.05. Because the pure-mean corner
`(0, 0, 1)` is always a candidate in the grid, validation RMSE is
guaranteed by construction to satisfy `pigauto_val_RMSE <= mean_val_RMSE`
on every trait. On held-out test, the invariant holds up to
sampling-noise slack (tolerance +2% in the production canary,
verified on five taxa by `script/regress.R`).

### Effect on existing benches

- Birds, mammals, fish, amphibians: lift unchanged (within 1%).
- Plants: on weak-signal traits (`height_m`, `leaf_area`, `seed_mass`,
  `sla`) the calibrator opens `r_MEAN > 0` and the prediction falls
  toward the grand mean, fixing the boundary-case regression where
  pigauto was previously 15-101% worse than mean; `wood_density` lift
  preserved.

### API

- `impute(..., safety_floor = TRUE)` — new arg, defaults to TRUE.
- `fit_pigauto(..., safety_floor = TRUE)` — same.
- The fit object gains four new slots: `r_cal_bm`, `r_cal_gnn`,
  `r_cal_mean` (each a named numeric of length `p_latent`), and
  `mean_baseline_per_col`.
- `multi_impute()` / `multi_impute_trees()` propagate the mean term
  automatically via `predict.pigauto_fit()`; no signature change.
- Legacy v0.9.1 fit objects keep working via `%||%` fallback in
  `predict.pigauto_fit()`: when the new slots are absent, the 2-way
  blend `(1 - r_cal) * BM + r_cal * GNN` is reconstructed from the
  scalar `r_cal`.

### Opt-out

`safety_floor = FALSE` reproduces the v0.9.1 1-D calibration
bit-identically on the legacy path.

### Testing

- `tests/testthat/test-safety-floor.R`: smoke canary, ~3 min on MPS.
  23 `test_that` blocks, 75 expectations covering: simplex grid,
  mean-baseline scalar, calibrate_gates 3-way path, fit_pigauto
  integration, predict 3-way blend + legacy fallback, multi_impute
  propagation, multi_impute_trees share_gnn propagation, AVONET300
  regression (continuous RMSE within +2%, discrete accuracy within
  -2pp), plants safety smoke (cached BIEN n=1000, RMSE within +15%
  of grand mean).
- `script/regress.R`: full canary, ~20 min on MPS at n=1000 per taxon.
  Writes `script/regress_result.json` with per-bench pass/fail.

### Files

- `R/calibration_three_way.R`: new — `simplex_grid()`,
  `mean_baseline_scalar()` helpers.
- `R/fit_helpers.R`: `calibrate_gates()` now returns list of 3 weight
  vectors and optionally searches the simplex grid.
- `R/fit_pigauto.R`: computes `mean_baseline_per_col` from training
  cells; passes `safety_floor` through; stores the new slots.
- `R/predict_pigauto.R`: 3-way blend with `%||%` backward compat.
- `R/impute.R`: pass-through arg + roxygen.
- `R/multi_impute.R`, `R/multi_impute_trees.R`: roxygen only.
- `inst/extdata/legacy_fit_v091.rds`: fixture for backward-compat
  test.
- `data-raw/make_legacy_fit_v091.R`: one-shot fixture builder.
- `specs/2026-04-23-safety-floor-mean-gate-design.md`: design spec.
- `plans/2026-04-23-safety-floor-mean-gate.md`: implementation plan.

# pigauto 0.9.1.9000 (dev)

## Median pooling for MI on count / proportion / zi_count traits (fixes #40)

- New argument `pool_method = c("median", "mean")` on `impute()` and
  `predict.pigauto_fit()`.  Default `"median"` takes the per-cell
  median of the M decoded draws when `n_imputations > 1` for **count**,
  **proportion**, and **zi_count** magnitude traits — robust to
  dropout-noisy latents that get amplified by `expm1()` / `plogis()`
  decoders.  `"mean"` restores pre-v0.9.2 arithmetic-mean pooling.
- **Why this matters**: the PanTHERIA full-scale bench at n = 4,027
  flagged a pathological RMSE blow-up on `litter_size` at
  `n_imputations = 20` (default **RMSE 7.7, em5 RMSE 60.0**) — a single
  dropout-noisy latent draw at ~+5 SD decoded to thousands and
  contaminated the mean.  With median pooling the RMSE returns to
  `n_imputations = 1` scale.  Coverage was unaffected (conformal ≥ 0.95,
  MC-dropout ≥ 0.92); this is purely a point-estimate fix.
- **Continuous / ordinal / binary / categorical** traits are
  unaffected — linear or probability-averaged decoders don't amplify
  latent outliers the same way.  `pool_method` has no effect on them.
- New test file `tests/testthat/test-count-mi-pooling.R` guards the
  robustness: 5-draw synthetic with one injected outlier of 1000 on a
  count trait returns median = 2 (expected) but mean ≈ 202 (bug).
## PanTHERIA mammal benchmark (second real-data validation)

Second real-data benchmark for pigauto after AVONET (birds). Tests the
"general-purpose phylogenetic trait imputation" claim on a different
taxon with a different trait mix.

- **`script/fetch_pantheria_and_tree.R`** — one-shot downloader. Fetches
  PanTHERIA (Jones et al. 2009, *Ecology*) from the ESA archive and
  builds a **taxonomy-derived cladogram** from the Order/Family/Genus/
  Species columns via `ape::as.phylo.formula()` + Grafen branch lengths.
  Not a time-calibrated phylogeny, but encodes the hierarchical
  taxonomic structure pigauto's phylo prior uses. Self-contained; no
  external tree-database dependency.
- **`script/bench_pantheria_full.R`** — full-scale bench on ~4,629 mammals
  × 8 canonical traits (body mass, head-body length, gestation length,
  litter size, max longevity, diet breadth, habitat breadth,
  terrestriality), MCAR 30% per trait. `pigauto_default` +
  `pigauto_em5` vs `mean_baseline`.
- **`script/bench_pantheria_bace_head_to_head.R`** — stratified 500-
  species subset with optional `BACE::bace()` for cross-package parity
  (graceful `requireNamespace` skip when BACE isn't available).
- **`script/make_pantheria_summary_html.R`** — aggregate landing page
  linked from the validation suite alongside the Phase 8 summary.
- Two new rows in the pkgdown validation suite table for the full
  bench + BACE head-to-head.
## Phase 8: discriminative benchmark suite (complete)

Five reproducible bench scripts + one aggregate report, all linked from
the pkgdown validation suite. No `R/` code changes; bench infrastructure
only.

- **`script/bench_signal_sweep.R`** — Pagel's λ ∈ {0.1, 0.3, 0.5, 0.7,
  0.9, 1.0} × mixed-type traits (2 cont + 1 bin + 1 cat K=3) × 4 methods
  × 3 reps. Headline: at λ = 1.0 pigauto reaches **94.4% binary accuracy
  vs `mean_baseline` 76.3%**; at λ = 0.1 all methods collapse near chance
  (expected — no signal).
- **`script/bench_bace_avonet_head_to_head.R`** — `pigauto_default` vs
  `pigauto_em5` vs optional `BACE::bace()` on `avonet300`/`tree300` with
  identical splits. Graceful `requireNamespace` guard when BACE is
  unavailable. Honest finding: Phase 6 EM does NOT universally beat the
  plug-in at n = 300 (default **80.5% vs em5 65.9%** on Trophic.Level),
  consistent with the "Calibration at small n" caveat from v0.9.1 — EM
  is opt-in for a reason.
- **`script/bench_correlation_sweep.R`** (Phase 8.1) — cross-trait
  correlation ρ ∈ {0, 0.2, 0.4, 0.6, 0.8} × 4 continuous traits × 3 reps.
  pigauto delivers consistent **6–10× RMSE reduction** over `mean_baseline`
  across all ρ (0.07–0.14 vs 0.64–0.80). Pearson r averages 0.97 at all ρ.
- **`script/bench_evo_model_sweep.R`** (Phase 8.2) — 4 evolutionary
  models × 3 reps. BM RMSE 0.13 (6× lift over mean); OU 0.32 (3×);
  **regime_shift 0.07 (14×)**; nonlinear 0.63 (1.6× — pigauto struggles
  on genuinely non-BM correlations, Pearson r drops to 0.69). Honest
  graceful-degradation story.
- **`script/bench_clade_missingness.R`** (Phase 8.3) — realistic MAR
  pattern via random-clade masking on focal trait. At every target_frac
  tested pigauto maintains **Pearson r ≥ 0.96** and RMSE **5–8× lower**
  than `mean_baseline`: 0.10→0.18 vs 0.99 (r 0.97); 0.25→0.14 vs 1.11
  (r 0.97); 0.40→0.16 vs 1.16 (r 0.96). Strong validation of pigauto's
  phylo prior in the realistic missingness regime (real databases
  undersample clades, not random species). Clade-picker caps at
  `target_frac` to avoid degenerate 100%-masked runs.
- **`script/phase8_summary.html`** — aggregate TL;DR + sub-report links
  + reproducibility block (package versions, `sessionInfo()`).
- Validation suite (`pkgdown/assets/validation_suite.html`) gains five
  new rows, one per bench, with graceful `pending` fallback on missing
  RDS files.
## Taxonomic-breadth benchmarks: 4 vertebrate classes + plant boundary case

Real-data benches now cover four vertebrate classes plus a plants
run that establishes the scope of phylogenetic imputation. All use
30% MCAR held-out, seed 2026, and report both conformal and
MC-dropout 95% coverage alongside RMSE / accuracy. The vertebrate
benches consistently lift point estimates by 27-75% RMSE /
10-40 pp accuracy; the plants bench is the first honest
boundary-case result, with a clean wood_density lift but r ≤ 0.21
on four other traits, attributable to weak phylogenetic signal in
pooled BIEN species means and random polytomy resolution in
V.PhyloMaker2 (details below). **Conformal coverage still lands at
0.90-0.98 on plants**, so the uncertainty story generalises across
the kingdom jump even where point estimates do not.

- **Birds** (`script/bench_avonet9993_bace_n3000.rds`, Vulcan L40S):
  AVONET n=3,000 subset x 7 traits. Mass RMSE -45%, Beak -55%,
  Tarsus -46%, Wing -57%; Trophic.Level accuracy +16 pp;
  Primary.Lifestyle +13 pp. Conformal 0.94-0.96.
  New driver `script/bench_avonet_full_local.R` runs the full
  n=9,993 locally on Mac MPS (Vulcan hit a predict-stage OOM at
  this scale); overnight job, HTML generator
  `script/make_bench_avonet_full_local_html.R` ready.
- **Mammals** (from `feature/bench-pantheria`): PanTHERIA n=4,027 x
  8 traits. **terrestriality accuracy 0.55 -> 0.95 (+40 pp)**;
  body_mass_g RMSE -27%, head_body_length -75%; diet/habitat
  breadth +10 pp each. Conformal 0.93-0.96.
- **Fish** (`script/bench_fishbase.rds`, Mac MPS):
  FishBase x fishtree n=10,654 x 6 traits.
  **Vulnerability RMSE -45%**, Length -38% (r=0.81), BodyShapeI
  accuracy +31 pp; Troph -23%, Depth -10%; Weight improved from
  r=0.04 (noise) to r=0.42 under the new median-pool MI fix.
- **Amphibians** (`script/bench_amphibio.rds`, Mac MPS): AmphiBIO
  n=5,237 x 2 continuous traits. **Body_size_mm RMSE 89 -> 44
  (-51%)**, Body_mass_g -27% (r=0.98). Continuous-only v1 --
  AmphiBIO binary / categorical columns hit a character/double
  type mismatch in the threshold-joint baseline, parked for
  v0.9.2. Calibrated gates = 0 on both traits: pigauto's GNN
  correctly recognises that the Level-C joint MVN baseline is
  already optimal, and still beats the mean baseline by 27-51%
  because the BM baseline captures cross-trait correlation +
  phylogenetic signal.
- **Plants** (`script/bench_bien.R`, honest boundary-case finding):
  BIEN x V.PhyloMaker2 at n=4,745 species, n_imputations=20.
  Only `wood_density` lifts (-6% RMSE, r=0.43); `height_m`, `sla`,
  `seed_mass`, `leaf_area` show r ≤ 0.21 and are 15-101% *worse*
  than grand mean on RMSE. Root cause is in the data, not pigauto:
  (i) BIEN species means are pooled from heterogeneous individual
  observations, diluting phylogenetic signal, and (ii) V.PhyloMaker2
  scenario S3 resolves within-genus polytomies randomly, so most
  species sit behind random branches relative to the Smith & Brown
  2018 backbone. pigauto's gate safety contains the damage: on
  weak-signal traits the calibrated gate closes to 0 and the
  pipeline degrades gracefully to the (also weak) BM prior rather
  than blowing up. **Conformal coverage still lands at 0.90-0.98
  on all 5 traits** -- the distribution-free uncertainty story
  holds even when point estimates are weak, separating the UQ
  narrative (always robust) from the point-estimate narrative
  (taxon-dependent). A first pass at `n_imputations = 1` showed
  `height_m` RMSE blowing up to 47.7 via Jensen exp()-decode bias;
  rerunning at `n_imputations = 20` activated the `dc8cffa`
  median-pool MI correction and dropped it to 15.15. The weak-
  signal result is what remains after that fix -- it is a property
  of BIEN x V.PhyloMaker2 at this species pool, not a pigauto
  regression. See `script/bench_bien.html` for the full honest
  reading and the in-page Jensen note.

Validation suite (`pkgdown/assets/validation_suite.html`) and the
paper-draft summary (`useful/results_summary.md`) collect all
numbers in one place.

## FishBase + fishtree benchmark (vertebrate breadth triad complete)

- New `script/bench_fishbase.R` benchmarks pigauto on the intersection
  of `fishtree::fishtree_phylogeny()` (Rabosky et al. 2018) and
  rfishbase's `load_taxa() + species() + ecology()` -- 10,654 species
  with 6 mixed-type traits (Length, Weight, BodyShapeI, DepthRangeDeep,
  Vulnerability, Troph). Completes the vertebrate breadth triad:
  birds (AVONET), mammals (PanTHERIA), fish (FishBase + fishtree).
- Full-scope result on Apple Silicon MPS (2.5 hr wall, post median-
  pool fix dc8cffa):
  Vulnerability RMSE 17.9 → 9.9 (-45%); Length RMSE 41.2 → 25.5
  (-38%, r=0.81); BodyShapeI accuracy 0.46 → 0.78 (+31 pp);
  Troph RMSE -23%; DepthRangeDeep RMSE -10%; Weight RMSE -3%
  (r=0.42 vs 0.04 under mean pool).
- HTML generator: `script/make_bench_fishbase_html.R` + pkgdown
  mirror at `pkgdown/assets/dev/bench_fishbase.html`.

## MI-pool robustness: median pooling for decode-amplifying traits

- `pool_imputations()` now uses **median pooling** across decoded
  values for all trait types where the decode step amplifies
  dropout-noisy draws: log-transformed continuous (new), count
  (extended), zi_count (extended), and proportion (new). Untransformed
  continuous keeps mean pooling (backward-compat default).
- Rationale: when `n_imputations > 1`, a few dropout-noisy latent
  draws with high-tail values decode via `expm1()` / `plogis()` to
  absurd magnitudes that dominate `rowMeans()`. Median is the
  outlier-robust pooling choice that matches PR #41's approach for
  count / proportion, now extended across the full set of
  decode-amplifying trait types.
- Concrete example: on the FishBase n=10,654 bench, pigauto's Length
  RMSE was 3,718 vs mean-baseline's 41 under the old mean-pool
  scheme -- a pure decode artefact (the calibrated gate on Length
  was 0). Median pool eliminates the artefact without retraining.
- Binary / categorical / ordinal pooling is unchanged (mean of
  probabilities / rowMeans of integer classes is correct there).
- 4 new regression tests in `test-mi-pool-robust.R` cover the four
  trait types that switched to median, verifying that a single
  extreme draw does not dominate the pooled output.

## GPU memory fix at predict stage (large n)

- `fit_pigauto()` now moves its returned `model_state` to CPU and
  calls `torch::cuda_empty_cache()` before returning, so the training
  run's activations, optimizer state, and per-epoch intermediate
  tensors can be reclaimed by torch's CUDA caching allocator. Prior
  to this fix, the fit object held GPU tensor refs that kept
  ~40 GB of training-time memory resident, and the downstream
  `predict()` path OOM'd on any card with <= 80 GB at n >= 5000.
- `impute()` now calls `gc()` + `torch::cuda_empty_cache()` between
  `fit_pigauto()` and `predict()` as belt-and-braces.
- Reproducer: Vulcan L40S (46 GB) at n=5000 or n=9993, default or
  compact config; predicted in 382 MB / 192 MB / 762 MB allocations
  on a 43 GB-already-used GPU. All four configurations on 2026-04-21
  died at the same step with the same root cause.
- No behavioural change for users: same `pigauto_fit` shape,
  `predict()` output identical; just no longer OOM at scale.

## Phase 7 EM: off-diagonal conditioning (opt-in, on top of Phase 6)

- New argument `em_offdiag = FALSE` on `impute()` and `fit_baseline()`.
  When `TRUE` AND `em_iterations >= 2L`, each binary/ordinal liability
  cell's prior at iteration `k + 1` is the full conditional-MVN
  `(mu, sd)` given the posterior liability of other traits at iteration
  `k`. This means observing one discrete trait shifts (not just
  tightens) the prior on correlated other traits — the off-diagonal
  entries of `Σ` now enter the E-step, instead of being ignored as in
  Phase 6.
- `em_offdiag = FALSE` default preserves Phase 6 behaviour byte-for-byte.
- `em_offdiag = TRUE` with `em_iterations < 2L` is silently clamped
  (iter 0 has no EM; iter 1 is plug-in with no previous Σ to condition
  on).
- **Scope**: binary + ordinal only. OVR categorical stays on Phase 6
  diagonal — going off-diagonal on OVR would re-couple the K classes
  and reintroduce the rank-(K − 1) singular-matrix instability that OVR
  was built to dodge.
- Missing-pattern handling: cells with all-NA other-trait posteriors
  fall back to the unconditional Phase 6 prior; cells with partial-NA
  use a restricted conditional on the observed subset.
- `em_state` gains `em_offdiag` field.

## Phase 6 EM for threshold-joint baseline (opt-in)

- New argument `em_iterations = 0L` on `impute()` and `fit_baseline()`.
  When `>= 1`, the threshold-joint path (binary + ordinal + OVR
  categorical) iterates: the BM rate `Σ` learned by
  `Rphylopars::phylopars()` at iteration `k` is fed back as the
  per-trait prior SD at iteration `k + 1`, up to `em_iterations` times
  (default 5 when enabled) or until
  `||Σ_k - Σ_{k-1}||_F / ||Σ_{k-1}||_F < em_tol` (default 1e-3). Closes
  the final open item from the v0.9.0 "Deferred to future releases" list.
- Default `em_iterations = 0L` preserves v0.9.1 output byte-for-byte.
- Under independent-discrete data the EM path converges in 1-2 iters
  at the plug-in answer (no free lunch, no harm). Under
  correlated-discrete data (ρ >= 0.5), the prior tightens and
  discrete-trait accuracy typically gains a few pp.
- Full off-diagonal conditioning on `Σ` is deferred to a future Phase 7.
  Phase 6 uses `sqrt(diag(Σ))` only.

---

# pigauto 0.9.1 (2026-04-19)

pigauto 0.9.1 is a consolidation release on top of 0.9.0. Two headline
improvements: (a) **tree-sharing GNN is now the default** in
`multi_impute_trees()`, which makes posterior-tree multiple imputation
cheap enough that the Nakagawa & de Villemereuil (2019) canonical
workflow (T = 50 posterior trees × m_per_tree = 1 = M = 50 pooled
datasets) is the default path at n = 10,000 species; (b) **opt-in
calibration smoothers** (bootstrap conformal, median-over-splits gate,
low-val-cells warning) for small-n regimes. This release also ships the
three post-0.9.0 feature branches (B1 soft-liability, B2 rate-aware
attention, B3 full-threshold ordinal baseline), scaling validation up
to 10,000 species, and a clean `R CMD check` (0 errors / 0 warnings /
1 note, down from 1 / 4 / 3 on v0.9.0).

## Calibration diagnostics and opt-in smoothing for small-n regimes

- New `fit_pigauto()` argument `min_val_cells = 10L`. Warns at fit time
  if any trait had fewer than this many validation cells available for
  gate calibration and conformal-score estimation. Default of 10 is the
  floor of pathological territory (operational target is 20-30).
- New optional `gate_method = "median_splits"` runs the half-A / half-B
  grid search on `gate_splits_B = 31L` random splits and takes the
  median best-gate. Small empirical benefit at small `n_val` (gate SD
  0.406 → 0.360 across 10 seeds on the calibration sim). Default stays
  `"single_split"` for backward compat.
- New optional `conformal_method = "bootstrap"` averages the conformal
  quantile across `conformal_bootstrap_B = 500L` bootstrap resamples of
  the val residuals. Empirically reduces conformal-score variance ~30%
  but **does not** reliably improve 95%-interval coverage stability on
  the evaluation sim. Shipped as experimental; default stays
  `"split"`.
- Roxygen for `fit_pigauto()` gains a `@section Calibration at small n:`
  explaining when the 95% intervals should be treated as approximate.
  Short version: at n_val < 20 the split-conformal guarantee degrades;
  increase `missing_frac` or collect more species rather than relying
  on estimator smoothing.

## Tree-sharing GNN is now the default in `multi_impute_trees()`

- New default `share_gnn = TRUE` in `multi_impute_trees()`. The GNN is
  trained once on a reference tree (MCC via `phangorn::maxCladeCred`,
  fallback `trees[[1]]`) and reused across posterior trees, with only
  the BM baseline recomputed per tree. ~10-15x speedup at n=10k x T=50
  makes the Nakagawa & de Villemereuil (2019) canonical workflow
  (T=50, m_per_tree=1, M=50 pooled via Rubin's rules) the default path.
- Default `m_per_tree` flipped from `5L` to `1L` to align with the N&dV
  2019 canonical workflow. Users with `T < 20` get a runtime warning
  suggesting they bump `m_per_tree` to keep Rubin's rules stable.
- Opt-out: `share_gnn = FALSE` restores the pre-v0.9.1 per-tree fit
  path for users who need exact per-tree model independence.
- New optional arg `reference_tree` lets users override the MCC choice.
- `predict.pigauto_fit()` gains a `baseline_override` argument (mostly
  internal — used by the shared-GNN path).
- New Suggests: `phangorn` (for `maxCladeCred()`).
- Tree-uncertainty propagation analysis: when the calibrated gate
  closes (the common case on real data — see the v0.9.0 validation
  suite), the shared-GNN approximation is lossless. When the gate is
  open, tree variance is slightly under-estimated in the GNN channel
  only — the baseline channel still carries it.

## Full threshold-model ordinal baseline (B3)

- Ordinal traits now use proper interval-truncated Gaussian E-step
  in the Level-C baseline instead of z-scored integer passthrough.
  Observing "class 3 out of 5" constrains the liability to the
  interval [threshold_2, threshold_3] rather than treating it as a
  point value. The liability joins the joint Rphylopars fit alongside
  continuous + binary traits.
- New internal decoder: `decode_ordinal_liability()` converts the
  Rphylopars liability posterior back to z-scored integer class for
  downstream GNN compatibility.
- `fit_baseline()` threshold-joint dispatch now fires when ordinal
  cols are present (in addition to binary / categorical). Ordinal
  cols that can't be populated (e.g., <2 observations) fall back to
  per-column BM unchanged.

## Rate-aware message passing via learnable per-head Gaussian bandwidth (B2)

- `GraphTransformerBlock` attention bias now uses learnable per-head
  Gaussian bandwidth computed from raw squared cophenetic distances:
  `bias_h = -D_sq / (2 * softplus(log_bw_h)^2)`. Each head can learn
  its own phylogenetic scale — tight for fast-evolving traits, broad
  for conserved traits. Replaces the fixed `adj_bias_scale * log(adj)`
  formulation.
- `build_phylo_graph()` returns `D_sq` (squared cophenetic distances)
  alongside `adj`. Persists through training (not freed like `D`).
- Backward compat: when `D_sq = NULL` (pre-B2 saved fits), the old
  `log(adj)` path fires unchanged. Legacy attention path (`use_transformer_blocks = FALSE`) is completely unaffected.

## Soft-liability E-step for multi-obs aggregation (B1)

- New argument `multi_obs_aggregation = c("hard", "soft")` on `impute()`
  and `fit_baseline()`. Default `"hard"` preserves v0.9.0 behaviour.
  `"soft"` uses Rao-Blackwell convex-combination E-step for binary and
  categorical traits: a species with 6/10 class-1 observations now
  produces a less-extreme posterior liability than 10/10, fixing the
  threshold-at-0.5 / argmax-one-hot information loss from Phase 10.
- New helpers: `estep_liability_binary_soft(p, mu_prior, sd_prior)` and
  `estep_liability_categorical_soft(p_vec, mu_prior, sd_prior)` in
  `R/liability.R`. Both reduce to their hard counterparts at boundary
  (one-hot) proportions; both return the prior mean at maximum ambiguity
  (p = 0.5 for binary, uniform for categorical).
- `aggregate_to_species()` gains `soft_aggregate` flag: when TRUE,
  preserves species-level proportions instead of thresholding/argmax.
- Benchmark: binary accuracy consistently improves +1.4–2.4pp across
  all 3 phylogenetic-signal regimes. Categorical shows +2.0pp at high
  signal but mixed results at moderate/low signal (hard threshold
  sometimes snaps to the correct majority class by luck).

## Scaling validation (up to 10,000 species)

- **AVONET 9,993 + BirdTree missingness sweep** (`script/bench_avonet_missingness.R`).
  Mean / mode vs BM baseline vs pigauto at 20% / 50% / 80% missing, run
  from Compute Canada (Narval). pigauto beats the BM baseline on every
  continuous trait at 80% missing; discrete accuracy stays within ≈1pp
  of BM across the sweep. Output: `bench_avonet_missingness.rds` +
  `bench_avonet_missingness.md`, rendered into the pkgdown validation
  suite.
- **Scaling curve to n = 10,000** (`script/bench_scaling_v090_extended.R`).
  Per-stage wall-clock + peak memory at
  n ∈ {100, 300, 1000, 3000, 5000, 7500, 10000}. Confirms sub-quadratic
  growth for the baseline and the expected O(n²) scaling of the GNN
  attention. At n = 10,000 the full pipeline completes in ~30–60 min on
  a single CPU node. Output: scaling table + curve PNG in the pkgdown
  validation suite.

## Benchmarks

- `bench_covariate_sim.R` rerun on v0.9.0 (previously pinned to v0.6.0).
  Now the reference entry for the environmental-covariates path in the
  validation suite.
- `bench_multi_obs_mixed.R` output regenerated with the B1 soft E-step
  path enabled; captures the +12–28pp categorical lift under soft
  aggregation vs hard.

## Documentation

- New pkgdown article **"Propagating Tree Uncertainty"**
  (`vignettes/tree-uncertainty.Rmd`). Decision guide for single-tree vs
  multi-tree MI, canonical N&dV 2019 workflow under the new
  `share_gnn = TRUE` default, worked example over posterior trees with
  `Map()`-based downstream-model refitting, and a timing table.
- pkgdown trait-type consistency pass: `zi_count` and `multi_proportion`
  are now listed everywhere alongside the other six trait types
  (Getting Started, Mixed Types, the validation suite, and the README).

## Internal

- **`R CMD check` clean**: 0 errors / 0 warnings / 1 note, down from
  1 / 4 / 3 on v0.9.0. The remaining note is the `BACE` in-tree
  "Package sources not in a standard location" flag, which is a
  deliberate `.Rbuildignore` choice so BACE ships as a comparison
  reference without being part of the installed pigauto package.
- **Compute Canada SLURM bundles**. `submit_v090_cloud/` (Narval;
  renamed clusters after the 2026 Alliance Canada infrastructure
  renewal — Cedar → Fir etc.) ships one-shot scripts for the AVONET
  missingness sweep and the scaling curve. `submit_v090_vulcan/` (PAICE
  account `aip-snakagaw`) ships a SLURM array driver for the 60-cell
  calibration-coverage grid.
- **pkgdown CI path filter**: the pkgdown workflow now only rebuilds on
  changes to user-facing files (`R/`, `vignettes/`, `pkgdown/`,
  `README.md`, `NEWS.md`, `DESCRIPTION`). Saves ~50% of the
  GitHub Actions minutes previously consumed by pigauto PR checks.
- **Gap fixes**: replaced a hardcoded worktree path in
  `bench_multi_obs_mixed.R` with a portable relative path; refreshed
  the bench_multi_obs_mixed `.rds`.

---

# pigauto 0.9.0 (2026-04-16)

pigauto 0.9.0 is a large release that ships two complementary upgrades:
a **Level-C phylogenetic baseline** (captures cross-trait correlation
and threshold structure), and a **Graph Transformer GNN backbone**
(replaces the 2-layer single-head attention stack with a proper
multi-head transformer). Both are gate-safe — pigauto still closes its
per-trait gate when the baseline is already optimal, so high-phylo-signal
datasets retain identical behaviour to v0.7.0.

## Level-C baseline improvements

A new family of joint multivariate baselines replaces per-column
Brownian motion + label propagation on datasets where cross-trait
correlation carries signal.

- **Liability encoding (`R/liability.R`)**. Every trait type now has an
  explicit liability interpretation: continuous types are their own
  liability; binary uses a 1D threshold at 0; ordinal uses K-1
  thresholds; categorical uses K liabilities with an argmax constraint.
  A single dispatcher `estep_liability(tm, observed, mu_prior, sd_prior)`
  returns posterior mean/variance of the underlying liability given an
  observation. This is the foundation the joint paths below build on.

- **Joint multivariate BM baseline (`R/joint_mvn_baseline.R`)**. When
  Rphylopars is installed and the dataset has ≥2 BM-eligible latent
  columns, the baseline now fits a single joint multivariate BM across
  continuous, count, ordinal, proportion, and ZI-magnitude cols via
  `Rphylopars::phylopars()`, capturing phylogenetic correlation across
  traits. Benchmark on correlated simulated BM data:
  **33.7% RMSE lift** vs per-column BM (`bench_joint_baseline.R`).

- **Threshold-joint binary baseline (`R/joint_threshold_baseline.R`)**.
  Binary traits get a truncated-Gaussian E-step posterior mean via
  `estep_liability_binary` and join the joint MVN alongside continuous
  traits. Decoded back to `logit(P(y=1))` so BCE and gate math are
  unchanged. Benchmark: 3/3 scenarios with ≥2pp accuracy lift over the
  LP baseline, peaking at **+16pp** at rho = 0.6
  (`bench_binary_joint.R`).

- **OVR categorical baseline (`R/ovr_categorical.R`)**. Each K-class
  categorical trait is decomposed into K independent "is_class_k vs
  rest" binary threshold fits, then normalised into a row-stochastic
  distribution. This is the same OVR strategy BACE uses. On AVONET 300,
  this lifts Trophic.Level accuracy to **77%** and Primary.Lifestyle to
  **84%** — beating BACE-OVR's 72% on both. The K-independent-fits path
  sidesteps the rank-(K-1) numerical instability that made earlier
  single-fit approaches unstable; each individual fit has only one
  categorical-related column.

- **ZI-count gate routing**. The binary gate col of zero-inflated count
  traits now joins the threshold-joint path alongside binaries (the
  magnitude col already rode the Phase-2 joint MVN).

- **Graceful fallback throughout**. Any Level-C path that can't fit
  (Rphylopars missing, too few observations per column, multi_proportion
  trait present) falls through to the legacy per-column BM + LP paths
  without user-visible breakage.

## GNN improvements

- **Graph Transformer backbone (`R/graph_transformer_block.R`)**. The
  2-layer single-head attention stack is replaced by `n_gnn_layers`
  instances of a `GraphTransformerBlock` — multi-head attention (default
  4 heads) with learnable per-head `log(adj + eps)` bias (so the phylo
  prior is respected by each head independently), feed-forward network
  (default 4× width), pre-norm, two residual skips. FFN output
  projection initialises to zero so each block starts near-identity,
  preserving pigauto's gate-closed-at-init safety.

- **New hyperparameters on `ResidualPhyloDAE`** and `fit_pigauto()`:
  `use_transformer_blocks = TRUE` (new default), `n_heads = 4L`,
  `ffn_mult = 4L`. The legacy architecture is fully preserved behind
  `use_transformer_blocks = FALSE` and reproduces pre-0.9.0 numbers
  exactly. Pre-0.9.0 saved `pigauto_fit` objects reconstruct via the
  legacy path automatically.

- **Empirical characterisation**. On the 4-scenario discriminative
  benchmark the transformer and legacy GNN produce identical predictions
  in high-signal regimes (gate closes to zero — architecture below the
  gate is irrelevant) and the transformer shows small RMSE lift in
  moderate/low signal scenarios where the gate stays partly open.
  Runtime is ~25–30% slower per fit. Larger transformer gains likely
  need self-supervised pretraining on a large tree corpus (future work).

## Multi-observation unlock

- **Multi-obs data now uses Level-C paths**. Previously the joint MVN,
  threshold-joint, and OVR categorical paths refused multi-obs input via
  four `!multi_obs` guards in `fit_baseline()`. Those guards are all
  removed; multi-obs data is aggregated to species level inside each
  Level-C helper via a new internal `aggregate_to_species()` helper.
  Aggregation rules: mean for continuous-family, threshold-at-0.5 for
  binary/ZI-gate, argmax-one-hot for categorical. Splits are
  mask-then-aggregated so val/test leakage is impossible; single-obs
  data passes through unchanged.

- **Benchmark (`bench_multi_obs_mixed.R`)**. On a synthetic multi-obs
  mixed-type simulation (150 species, 2 continuous + 1 binary + 1 K=4
  categorical, Poisson-5 obs per species), the Level-C path beats
  the LP path by **+12–28pp categorical accuracy** and **+5–7.5pp
  binary accuracy** across all three phylogenetic-signal regimes.

- **Aggregation caveat**. Binary/categorical aggregation is lossy —
  a species with split observations (e.g., 6/10 class-1) becomes a
  single class-1 observation at species level. For datasets with
  strong within-species variability on discrete traits, prefer
  single-obs mode with a user-computed species-level summary.

## Benchmarks

Nine new or reran benchmark reports quantifying the release's impact:

- `bench_joint_baseline.R` — Phase-2 joint MVN A/B on correlated BM data.
- `bench_binary_joint.R` — Phase-3 threshold-joint A/B on binary traits.
- `bench_discriminative.R` + `bench_discriminative_phase9.R` — 4-scenario
  mixed-signal sweep, including the transformer vs legacy GNN A/B.
- `bench_avonet_phase6.R` — post-Phase-6 AVONET 300 vs LP baseline.
- `bench_categorical_joint.R` — OVR K-fits A/B on synthetic categorical
  data.
- `bench_covariate_sim.R` — rerun on 0.9.0. Shows **−10.6% RMSE** on the
  low-phylo-strong-env scenario vs pre-Phase-6. Pre-0.9.0 snapshot in
  `bench_covariate_sim_preP6.md`.
- `bench_multi_obs.R` — rerun on 0.9.0. Numbers unchanged from pre-0.9.0
  (bench is single-trait so it can't exercise Level-C dispatch
  conditions; `bench_multi_obs_mixed.R` is the correct validator).
- `bench_multi_obs_mixed.R` — new synthetic mixed-type multi-obs
  benchmark specifically to validate Phase 10's multi-obs unlock.

## Internal

- `joint_mvn_available()` gates the Level-C dispatch — returns `TRUE`
  iff Rphylopars is installed. All Level-C paths fall back gracefully
  when it returns `FALSE`.
- `aggregate_to_species()` is the sole multi-obs → species-level helper
  and is the right place to hook future soft-liability aggregation
  (currently threshold-at-0.5 for binary, argmax for categorical).
- 20+ new test blocks across `test-liability.R`,
  `test-joint-threshold-baseline.R`, `test-joint-baseline.R`,
  `test-ovr-categorical.R`, `test-graph-transformer-block.R`,
  `test-phase9-integration.R`, and `test-multiobs-levelc.R`. Full suite
  is 778 PASS, 0 FAIL.

## Deferred to future releases

- Soft-liability E-step for multi-obs aggregation (addresses the
  threshold-at-0.5 lossiness documented above).
- Rate-aware message passing in the transformer blocks.
- Full threshold-model ordinal baseline (ordinal currently rides the
  Phase-2 joint MVN via its z-scored integer encoding).
- Self-supervised pretraining of the transformer backbone on a large
  tree corpus.
- CRAN submission preparation (BACE in-tree separation, torch handling,
  absolute-path cleanup in `script/`).

# pigauto 0.7.0

## New trait type: multi_proportion (compositional data)

- **`multi_proportion`**: a new trait type for compositional data — K columns
  per row that sum to 1 (e.g. plumage-colour proportions, diet composition,
  microbiome relative abundances, allele frequencies). Declared via a new
  `multi_proportion_groups = list(colour = c("black", ..., "yellow"))`
  argument to `impute()`, `preprocess_traits()`, and `multi_impute()`.

- **Encoding**: centred log-ratio (CLR) + per-component z-score. Small
  epsilon (`1e-6`) is applied to zero components before log to keep the
  transform safe; rows are re-normalised so the composition is preserved.

- **Baseline**: K independent Brownian-motion imputations on the CLR
  columns. BM is well-defined on CLR space (Euclidean), unlike the raw
  simplex.

- **GNN output**: K logits projected to sum-zero, then softmax at decode
  time — guaranteed to lie on the simplex.

- **Loss**: MSE on CLR columns (averaged across K), applied group-wise so
  the whole composition is masked together during DAE corruption.

- **Metrics**: Aitchison distance (the natural compositional-data metric),
  RMSE on CLR space (comparable to continuous traits), and simplex MAE
  (mean absolute error on decoded proportions).

- **Benchmark**: `script/bench_multi_proportion.R` with 5-scenario
  signal sweep and K ∈ {3, 5, 8, 12} secondary sweep. Report rendered
  by `script/make_bench_multi_proportion_html.R`.

- **Tests**: `tests/testthat/test-multi-proportion.R` (5 test blocks,
  20 expectations) — encoding, validation, baseline, end-to-end `impute()`.

## Internal

- `evaluate_imputation()` output adds three columns: `aitchison`,
  `rmse_clr`, `simplex_mae`. Existing metrics columns are unchanged.

- `impute()` gains explicit `trait_types` and `multi_proportion_groups`
  arguments (previously `trait_types` was documented as passed through
  `...` but was not actually forwarded; this fixes that).

# pigauto 0.6.2

## Documentation

- **Tree uncertainty workflow clarified.** `multi_impute_trees()`,
  `trees300`, the README, and `vignette("getting-started")` now all
  describe the two-step workflow explicitly: step 1 is tree-aware
  imputation (pigauto's job — already done correctly); step 2 is
  tree-aware downstream analysis (the user's responsibility,
  following Nakagawa & de Villemereuil 2019, *Syst. Biol.* 68:632–641).
  Every doc now includes a complete `Map()`-over-`mi$tree_index` code
  example showing how to run the corresponding tree in the downstream
  model.
- **Compute-cost table added** to the tree-uncertainty docs: wall-clock
  budgets at `n` = 300 / 5,000 / 10,000 species with `T` = 10 / 50
  trees, plus guidance on reducing `T` for large trees and
  parallelising across machines.

## Internal

- No API changes; no numerical changes to any fitted model. Doc-only
  patch release.

# pigauto 0.6.1

## Bug fixes

- **MC dropout zero variance fixed** (`R/predict_pigauto.R`): when
  `n_imputations > 1`, the blend equation now uses a BM posterior draw
  `t_BM_draw ~ N(BM_mu, BM_se)` instead of the deterministic `t_MU`.
  When the calibrated gate is zero (BM dominates), each imputation now
  samples from the correct BM posterior rather than returning the same
  point estimate, giving non-zero between-imputation variance.

- **Conformal draw scale bug fixed** (`R/multi_impute.R`): for
  log-transformed traits (Mass, Wing length, etc.), draws are now made
  on the latent z-score scale and back-transformed, avoiding near-zero
  variance caused by dividing a log-scale SE by an original-scale value.
  The same fix applies to log1p-scaled counts and logit-scaled
  proportions.

## New features

- `draws_method = "conformal"` is now the default for `multi_impute()`.
  Draws sample from the split-conformal calibrated uncertainty
  distribution (better calibrated than MC dropout for downstream
  Rubin's-rules pooling).
- `draws_method = "mc_dropout"` remains available and is now correct
  (see bug fix above).

## Hex sticker

- New `man/figures/logo.png`: pink pig driving a car with a 7-tip
  ultrametric pectinate cladogram above the pig's head. Forest green /
  mint palette (v3) selected as the official sticker.

# pigauto 0.6.0

## Observation-level covariate refinement for multi-obs data

When datasets have multiple observations per species measured under
different conditions (e.g. CTmax at different acclimation temperatures),
the GNN now produces **covariate-conditional predictions within species**.
A new refinement MLP (`obs_refine`) re-injects user-supplied covariates
after species-level phylogenetic message passing, so that different
observations of the same species receive different predicted values
based on their covariate context. This is controlled by a residual
connection that preserves the phylogenetic signal from message passing.

- `R/model_residual_dae.R`: new `n_user_cov` parameter and `obs_refine`
  MLP (activated only when `species_col` + `covariates` are both supplied)
- `R/fit_pigauto.R`: `n_user_cov` stored in model config
- `R/predict_pigauto.R`: backward-compatible model reconstruction

## New bundled dataset: `ctmax_sim`

Simulated multi-observation-per-species CTmax data (1,464 observations
across 300 species, with acclimation temperature as an observation-level
covariate). 30% of species are entirely unobserved. Uses `tree300`.
See `data-raw/make_ctmax_sim.R` for the generation script.

## New benchmark: multi-observation imputation

`script/bench_multi_obs.R` simulates CTmax-like datasets under varying
phylogenetic signal and within-species acclimation response ratios,
comparing species_mean, pigauto (no covariates), and pigauto (with
covariates). The observation-level refinement gives measurable RMSE
improvement when within-species covariate effects are strong.

---

# pigauto 0.5.0

## Three sources of information

pigauto now explicitly combines three sources of information for
imputation: (1) the phylogenetic tree, (2) cross-trait correlations,
and (3) optional environmental covariates. The `covariates` argument
is accepted by `impute()`, `multi_impute()`, `multi_impute_trees()`,
and the lower-level pipeline functions. Covariates are threaded
through the GNN with the same gated safety that protects against
GNN degradation — they only contribute when they demonstrably
improve accuracy.

## Internal BM baseline (Rphylopars no longer required)

The Brownian-motion baseline for continuous, count, ordinal, and
proportion traits is now computed internally using conditional
multivariate normal imputation on the phylogenetic correlation matrix
`R = cov2cor(vcv(tree))` (Goolsby et al. 2017). This removes the
hard dependency on `Rphylopars`, which is now in `Suggests` only.
The internal implementation (`R/bm_internal.R`) uses univariate
per-column imputation with Cholesky decomposition and nugget
regularisation for near-singular submatrices. Validation against
`Rphylopars` shows r = 0.97–0.98 (expected given univariate vs
multivariate difference).

## Tree-uncertainty MI via `multi_impute_trees()`

New function `multi_impute_trees()` performs multiple imputation
across a posterior sample of phylogenetic trees, so that phylogenetic
uncertainty propagates into downstream standard errors via Rubin's
rules. Demonstrated with the bundled `trees300` dataset (10 BirdTree
posterior trees). Benchmark results show SE inflation of 1.1–2.1x
and fraction of missing information (FMI) rising from ~0.03
(single-tree) to 0.22–0.79 (multi-tree) depending on missingness
level.

## New bundled datasets

- `trees300`: 10 posterior trees from the BirdTree Hackett backbone,
  pruned to the same 300 tips as `tree300`
- `delhey5809`: plumage lightness and 4 environmental covariates
  (absolute latitude, mean temperature, mean precipitation,
  tree cover) for 5,809 bird species (Delhey 2019)
- `tree_delhey`: matching BirdTree phylogeny for `delhey5809`

## Benchmark suite

Thirteen benchmark reports are now accessible from the pkgdown site
under Articles > Benchmarks:

- Per-type benchmarks: continuous, binary, ordinal, count,
  categorical, proportion, zero-inflated counts
- Missingness mechanisms: MCAR vs MAR (trait/phylo) vs MNAR
- Tree uncertainty: single-tree vs multi-tree MI
- Environmental covariates: Delhey 5,809-species real data
  (high phylo signal, covariates give ~zero lift) and simulation
  sweep (4 lambda x beta scenarios; 8–15% RMSE reduction when
  env effects are strong, ~0% when phylo signal dominates)

Each report is a self-contained HTML page generated from
`script/make_bench_*_html.R` and backed by reproducible
`script/bench_*.R` drivers.

## Documentation

- All user-facing prose updated from "two sources" to "three
  sources" of information
- Rphylopars references replaced with "phylogenetic BM" or
  "internal conditional-MVN baseline" throughout
- `DOCS.md` updated with new benchmarks, datasets, and functions

---

# pigauto 0.4.0

## Breaking: TabPFN baseline removed

`fit_baseline_tabpfn()` and `setup_tabpfn()` are gone, along with the
`reticulate` suggested dependency, the `script/bench_tabpfn.*` files,
`script/run_tabpfn.py`, `tests/testthat/test-tabpfn.R`, and the
`dev/bench_tabpfn.html` benchmark page. The TabPFN wrapper was
originally added as a curiosity — a test of whether a general-purpose
tabular transformer could compete with pigauto as a drop-in baseline.
In practice it is a poor fit for what pigauto is trying to do:

1. **It is not phylogenetic.** TabPFN treats rows as exchangeable and
   cannot use the tree at all. The whole point of pigauto is that the
   tree carries information about trait evolution; feeding traits
   through a tree-blind model throws that away.
2. **Quadratic attention does not scale.** TabPFN's attention is
   O(n²) in the number of rows. At `n = 9,993` species (the full
   AVONET dataset) this is already near the edge of what a single
   GPU can handle, and most users of pigauto care precisely about
   scales at which a baseline must be cheap.
3. **The cross-language subprocess architecture was fragile.** R
   `torch` and Python `torch` share `libtorch` and cannot coexist in
   the same process, so TabPFN had to run via `system2()` on a
   separately-installed Python venv. Every time the user's Python,
   `tabimpute`, or `torch` versions drifted, the wrapper broke. This
   is not maintenance the project can afford.
4. **Keeping the benchmark page gave a misleading impression of
   relevance.** Users who hit the "pigauto vs TabPFN" page came away
   thinking TabPFN was a first-class supported baseline in pigauto.
   It never was.

Users who still want to benchmark pigauto against TabPFN can install
`tabimpute` separately and call it directly from Python. The removal
simplifies pigauto's installation, CI, and documentation.

## Docs cleanup (carried over from v0.3.3 work)

- `DOCS.md`: rewritten benchmark and developer-reports sections to
  reflect the v0.3.2+ benchmark set (scaling to 10,000 tips and the
  AVONET missingness sweep). Dropped the obsolete `benchmark_report`
  and `benchmark_final_report` links that v0.3.2 had already pruned
  from `_pkgdown.yml` but had not yet removed from `DOCS.md`.
- `CLAUDE.md`: removed the TabPFN Python-venv setup recipe and the
  "R torch ↔ Python torch cannot share a process" gotcha (no longer
  applicable). Version line bumped to 0.4.0.

# pigauto 0.3.3

## Correctness fix: Vphy must be a correlation matrix in `glmmTMB`'s `propto()`

The multiple-imputation example in `README.md`, the two HTML
tutorials (`pigauto_intro.html`, `pigauto_workflow_mixed.html`), and
`tests/testthat/test-multi-impute.R` all constructed the phylogenetic
random-effect structure with

```r
Vphy <- ape::vcv(tree)            # WRONG: covariance, diag = tree height
glmmTMB(... + propto(0 + species | dummy, Vphy), ...)
```

`propto()` estimates `sigma^2` freely, so passing the raw covariance
matrix (whose diagonal is the tree height) silently rescales the
phylogenetic variance component by tree height. The correct usage
passes a *correlation* matrix:

```r
Vphy <- cov2cor(ape::vcv(tree))   # diag = 1
```

All four sites are now fixed. Thanks to the user who flagged this.

## Honest framing: the GNN is a gated ensemble, not a residual model

The user-facing prose throughout the package historically described
pigauto as a *Residual Phylogenetic Denoising Autoencoder
(ResidualPhyloDAE) that learns a residual correction on top of a BM
baseline*. This wording was misleading for binary and categorical
traits in two ways:

1. The "baseline" for those types is **phylogenetic label propagation**
   (a similarity-weighted average of observed labels using a Gaussian
   kernel on the cophenetic distance) — not a generalised linear model
   that produces residuals on a link scale.
2. The GNN is trained **end-to-end** by minimising the type-appropriate
   loss (MSE for continuous/count/ordinal, BCE for binary,
   cross-entropy for categorical) on the **blended output**, not on
   `y - baseline`. So the GNN's output `delta` is a full per-cell
   prediction, not a statistical residual.

What pigauto actually is, in honest terms, is a **gated ensemble** of
two predictors:

- A phylogenetic **baseline** (Brownian motion via Rphylopars for
  continuous/count/ordinal traits; phylogenetic label propagation for
  binary/categorical traits).
- A **graph neural network correction** (an internal torch module that
  produces a full per-cell prediction `delta` from spectral node
  features and an attention-biased message passing over the
  phylogeny).

Prediction is the per-trait blend
`(1 - r_cal) * baseline + r_cal * delta`, where `r_cal` is calibrated
on a held-out validation split (with a split-validation cross-check
for discrete traits). When the baseline is already optimal the gate
closes and the GNN becomes a no-op.

The internal torch class name `ResidualPhyloDAE` is preserved because
its GNN layers use ResNet-style residual skip connections — that is a
correct, well-defined ML term. But user-facing prose no longer
describes the GNN as "learning a residual on top of the baseline".

Files updated: `README.md`, `DESCRIPTION`, `CLAUDE.md`, `_pkgdown.yml`,
`DOCS.md`, `vignettes/getting-started.Rmd`, `vignettes/mixed-types.Rmd`,
`R/fit_pigauto.R` (roxygen), `R/predict_pigauto.R` (roxygen),
`R/model_residual_dae.R` (clarifying header comment),
`R/simulate_traits.R`, `R/plot.R` (plot titles), and
`script/make_intro_html.R`.

## Tracking issue for v0.4.0 architectural work

The honest fix above is a documentation correction, not a model
change. A genuinely *residual* baseline for binary and ordinal traits
would require a liability-scale model (threshold BM, `phylolm::phyloglm`,
or `MCMCglmm` with `family = "threshold"`) — a design decision that
should be made deliberately rather than rolled into a same-day patch.
A v0.4.0 tracking issue lays out the trade-offs.

# pigauto 0.3.2

## New bundled dataset: `avonet_full` / `tree_full`

The full AVONET3 + BirdTree Stage2 Hackett MCC phylogeny, aligned to
9,993 species with 7 mixed-type traits (4 continuous morphometric +
2 categorical + 1 ordinal) and native missingness preserved.
Complements the 300-species `avonet300` / `tree300` bundled data:

- `avonet300` / `tree300` remain the right choice for quick examples,
  unit tests, and the getting-started vignette.
- `avonet_full` / `tree_full` are the right choice for scale
  benchmarks, realistic-missingness experiments, and anything that
  wants to exercise the v0.3.1 sparse Lanczos / cophenetic caching
  code paths on a real-world phylogeny.

The schemas are identical, so any code that runs on `avonet300` runs
on `avonet_full` with no modification. Combined tarball increment is
~328 KB (xz-9).

## New benchmark: AVONET missingness sweep (20 / 50 / 80%)

*Articles -> Benchmarks -> AVONET missingness sweep* is a new pkgdown
page comparing three methods on the full `avonet_full` dataset at
three MCAR missingness levels:

- **mean / mode** — column mean (continuous, ordinal) or class
  frequency (categorical). No phylogeny.
- **BM baseline** — Rphylopars Brownian-motion for continuous /
  ordinal, phylogenetic label propagation for discrete traits.
- **pigauto** — full pipeline (BM baseline + calibrated GNN delta
  + conformal intervals).

Per-trait RMSE, Pearson r, Spearman rho, and accuracy on held-out
test cells. Single seed per cell, trait-level masking so categorical
traits are held out as whole groups. Hyperparameters match the
`validate_avonet_full.R` scaling run for comparability.

## Docs cleanup

- Regenerated `pigauto_intro` and `pigauto_workflow_mixed` HTML
  tutorials against the v0.3.2 package state.
- Removed stale v0.3.0-era benchmark pages from the navbar
  (`benchmark_report.html`, `benchmark_final_report.html`,
  `test_report.html`) and from `pkgdown/assets/dev/`.
- Pruned 45 orphan files from `script/` (trial artefacts from
  pre-v0.3.1 benchmarking phases: `bench_v2/v3/v4.*`,
  `benchmark_avonet2000.*`, `benchmark_final.*`,
  `benchmark_simulation.*`, scaling duplicates, and Apr-6 demo
  artefacts).
- `README.md` citation now points at 0.3.2.
- `DESCRIPTION` updated to mention both bundled datasets.

# pigauto 0.3.1

## Scaling to 10,000 tips

`pigauto` now runs end-to-end on real 10,000-tip phylogenies. Before this
release the pipeline was tested and tuned around the 300-species bundled
AVONET data, and at 10,000 tips it wedged on a single dense eigendecomposition
of the graph Laplacian (hours of CPU time). Two targeted fixes land in this
release; no rewrite was needed.

### Fix A: sparse Lanczos eigensolver for `build_phylo_graph()`

`build_phylo_graph()` now uses `RSpectra::eigs_sym()` sparse Lanczos to
compute the `k + 1` smallest Laplacian eigenvectors when the tree has more
than 7,500 tips and **RSpectra** is installed. This replaces
`base::eigen(L, symmetric = TRUE)`, which was O(n^3) in time and discarded
>99% of its work (we only need the bottom ~32 eigenvectors).

- On the real 9,993-species AVONET + BirdTree Stage2 Hackett MCC phylogeny,
  the graph stage now takes under a minute, compared with roughly seven
  minutes of dense eigen() on the same hardware.
- The old dense path is kept as the default for smaller trees and as a
  safety fallback for when RSpectra is unavailable or Lanczos does not
  converge. The threshold is conservative (7,500) because highly symmetric
  ultrametric simulations (`ape::rcoal`) have degenerate spectra that need
  large Krylov subspaces to resolve.
- **RSpectra** is listed in `Suggests`. Installing it is strongly recommended
  for any tree larger than ~7,500 tips.

### Fix B: cache the cophenetic distance matrix across the pipeline

The cophenetic distance matrix `D = ape::cophenetic.phylo(tree)` was previously
computed four separate times in a single `impute()` call (three times inside
`build_phylo_graph()` and once inside `fit_baseline()`). Each call allocates a
dense n*n double matrix and runs an O(n^2) tree traversal -- wasted work at
scale. In 0.3.1:

- `build_phylo_graph()` computes `D` exactly once and returns it as part of
  the result list (`graph$D`).
- `fit_baseline()` gains an optional `graph` argument; when supplied, it
  reuses `graph$D` instead of recomputing the cophenetic matrix.
- `impute()` and `fit_pigauto()` pass the graph through automatically, so
  users get the speedup without any code changes.
- `impute()` also explicitly releases `graph$D` after `fit_baseline()`
  finishes, so the ~800 MB cophenetic matrix does not sit in R memory
  during training (a subtle side-effect that caused a large per-epoch
  regression at n = 10,000 until it was dropped).

### Scaling benchmark

See `pkgdown/assets/dev/scaling.html` (linked from the site as
*Articles -> Benchmarks -> Scaling to 10,000 tips*) for a side-by-side
comparison of v0.3.0 vs v0.3.1 at n in {300, 1k, 2k, 3k, 5k, 7.5k, 10k},
plus an end-to-end validation run on the full AVONET3 + BirdTree dataset
(9,993 species, 7 mixed-type traits). The pkgdown article also walks
through both a quick 300-species example with the bundled `avonet300`
data and a full 10k-species example using the BirdTree Stage2 Hackett
MCC phylogeny.

## Bug fixes

- `build_phylo_graph()`: the N > 2000 memory warning was out of date (the
  real blocker above that size was wall-clock time, not memory). The
  warning is now issued at N > 10000 and its text points at the true
  bottleneck.

# pigauto 0.3.0

## Major new features

### Attention-based message passing
The GNN now uses learned attention weights (with phylogenetic adjacency
as positional bias) instead of fixed adjacency matrix multiplication.
The model learns which phylogenetic neighbours are most informative for
each species, while the log-adjacency bias ensures initialisation
respects evolutionary distances. Controlled via `use_attention = TRUE`
(now the default).

### Validation-calibrated gates
After training, per-trait gate values are optimised on the validation set
via grid search. This prevents the GNN from adding noise when the
Brownian motion baseline is already optimal. On the AVONET benchmark,
this eliminates the RMSE regression seen in earlier versions.

### Conformal prediction intervals
Distribution-free 95% prediction intervals computed via split conformal
prediction on the validation residuals. Coverage is guaranteed under
exchangeability. Available for continuous, count, and ordinal traits via
`pred$conformal_lower` and `pred$conformal_upper`.

### Phylogenetic label propagation
Binary and categorical trait baselines now use phylogenetic similarity-
weighted label propagation instead of flat population-level frequencies.
Each species gets a personalised baseline reflecting its phylogenetic
neighbourhood.

### Adaptive spectral encoding
`k_eigen` now defaults to `"auto"`, scaling with tree size:
`min(max(ceil(n/20), 4), 32)`. This gives 4 features for small trees,
15 for 300 tips, and 32 for 640+ tips.

### Auto-generated HTML reports
`pigauto_report()` produces self-contained HTML reports with interactive
Chart.js visualisations: per-trait metrics, calibrated gate values,
training history, and conformal coverage.

### Evaluation framework
- `evaluate()`: per-trait metrics on held-out data
- `compare_methods()`: automated BM vs GNN comparison
- `cross_validate()`: proper k-fold cross-validation
- `summary.pigauto_fit()`: formatted console summary

### Plotting
- `plot.pigauto_fit()`: training history, gate values, conformal scores
- `plot.pigauto_pred()`: scatter plots, conformal intervals, probability
  distributions
- `plot_comparison()`: forest-plot style method comparison

### Simulation benchmark
`simulate_benchmark()` runs a complete simulation study with configurable
scenarios (BM, OU, regime shift, non-linear, mixed types), returning tidy
results with print/summary/plot methods. Useful for methods papers and
for verifying pigauto's behaviour on data with known properties.

## Minor improvements

- Categorical gate calibration now works at the trait level (all K one-hot
  columns share a single gate) using cross-entropy loss, not per-column MSE
- `build_phylo_graph()` default `k_eigen` changed from `8L` to `"auto"`
- `fit_pigauto()` default `k_eigen` changed from `8L` to `"auto"`
- `all_grads_finite()` now handles non-standard tensor types gracefully
- `impute()` uses adaptive k_eigen by default
- Package passes `R CMD check` with 0 errors and 0 notes
- Bumped version to 0.3.0

# pigauto 0.2.0

- Added support for multiple observations per species (`species_col`)
- Added mixed trait type support (continuous, binary, categorical,
  ordinal, count)
- Added AVONET 300 bundled dataset with mixed types
- Extended simulation benchmark (16 scenarios via BACE)

# pigauto 0.1.0

- Initial release
- ResidualPhyloDAE architecture
- Brownian motion baseline via Rphylopars
- Per-column gated blending
- MC dropout for multiple imputation
