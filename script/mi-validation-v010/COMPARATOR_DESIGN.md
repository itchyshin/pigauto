# Positive-control design for the MI correctness gate

This is a harness-only diagnostic response to the failed v0.10.0 MI campaign.
It does not add a public estimator. The purpose is to prove that the frozen
downstream gate is attainable before spending on another pigauto campaign.

## A — Aims

Determine whether (1) the exact known-DGP conditional distribution and (2)
established substantive-model-compatible MI recover the downstream `x` and `z`
fixed effects under the existing Gaussian, binomial, and random-intercept DGPs.

## D — Data-generating mechanism

The DGP, MAR mask, sample sizes, regimes, and release gates remain those in
`README.md` and `0_prepare.R`. Exactly 30% of `x` is selected as missing using
fully observed `y` and `z`.

Before standardisation, the covariate model is

`x_raw,i = a(regime, phylogeny_i, z_i) + epsilon_x,i`,
`epsilon_x,i ~ Normal(0, sigma_x^2)`.

The realised standardisation constants are retained, giving the conditional
prior `x_i ~ Normal(mu_x,i, s_x^2)`. For Gaussian outcomes,

`y_i = 1 + beta_x x_i + beta_z z_i + b_species[i] + epsilon_y,i`,
`epsilon_y,i ~ Normal(0, 1)`.

The oracle conditional sampler uses the resulting conjugate Normal posterior.
For the logit DGP it samples from the Normal prior and accepts with the exact
Bernoulli likelihood. The mixed-DGP oracle conditions on the simulated random
intercept; this is intentionally a positive control, not a deployable method.

The standard analysis-compatible comparator dispatches by frozen analysis:

- Gaussian `lm`: proper Bayesian normal-regression MI from
  `x ~ y + z + I(z^2)`, using the Jeffreys posterior
  `sigma^2 = RSS / chi-square(n_obs - p)`,
  `beta | sigma^2 ~ Normal(beta_hat, sigma^2 (X'X)^-1)`, followed by an
  independent residual draw for every missing cell;
- logit `glm`: `smcfcs` with the substantive model `y ~ x + z`;
- Gaussian random-intercept `lmer`: `jomo::jomo.smc(model = "lmer")`.

The Gaussian branch is deliberately narrow: it supports this complete-rank
Gaussian main-effects validation cell and is not a general SMC implementation.
All branches see `z^2` as a fully observed auxiliary. `jomo` is not used for
the logit comparator because its documented binomial substantive model uses a
probit link.

## E — Estimands

The estimands remain `beta_x` and `beta_z` from `lm(y ~ x + z)`,
`glm(y ~ x + z, binomial())`, and
`lmer(y ~ x + z + (1 | species))`. Variance components remain descriptive and
are never Rubin-pooled.

## M — Methods

The first run is an economical controls-only pilot: 50 replicates per each of
the six DGP-by-regime cells, `M = 50`, paired seeds, and no pigauto training.
Only if both controls pass the directional gate across all 12 fixed-effect cells
does implementation proceed to an analysis-aware package backend.

## P — Performance measures

The metrics remain bias, empirical SD, mean pooled SE, SE/empirical-SD ratio,
95% coverage, finite-valid Rubin rate, and downstream fit rate. The bias gate is
the paired mean imputation-minus-complete-data-oracle difference divided by the
oracle empirical SD. Absolute truth-based bias remains reported, but cannot be
the imputation gate because the frozen 500-replicate complete-data logit oracle
itself has finite-sample standardized bias above 0.10 in one cell. Coverage and
SE calibration remain truth-based. The frozen campaign cannot satisfy the
prospectively amended 1,000-replicate design requirement, so it cannot authorize
release.

## Symbolic alignment table

| Symbol | Harness object | Imputation role | Recovery extractor | Truth |
|---|---|---|---|---|
| `x_i` | completed column `x` | missing covariate draw | downstream coefficient `x` | 0.60 (`lm`/`lmer`), 0.80 (`glm`) |
| `y_i` | observed column `y` | substantive likelihood | downstream outcome | not an imputation estimand |
| `z_i` | observed column `z` | substantive covariate | downstream coefficient `z` | -0.40 (`lm`/`lmer`), -0.50 (`glm`) |
| `z_i^2` | auxiliary `z_sq` | imputation model only | not extracted | not an estimand |
| `mu_x,i`, `s_x` | DGP prior fields | oracle conditional sampler | internal | known DGP values |
| `b_species` | stored simulated intercept | oracle mixed conditional sampler | internal | known DGP value |
| `tau^2`, `sigma^2` | `VarCorr()` diagnostics | descriptive only | mean across fits | 0.64, 1.00 |

## Williams 11-item self-audit

| Item | Coverage |
|---|---|
| 1. Aims | Gate-attainability question stated above |
| 2. DGP | Frozen DGP and MAR selection retained |
| 3. Estimands | `beta_x`, `beta_z`; variance components diagnostic only |
| 4. Methods | Exact oracle, Bayesian normal-regression MI, SMCFCS, and jomo |
| 5. Performance | Frozen bias, SE, coverage, validity thresholds |
| 6. Environment | Git SHA and dependency versions stored per task |
| 7. Code | Three-stage harness in this repository |
| 8. Reproducibility | Deterministic manifest and method seeds |
| 9. Applied example | Out of scope for this diagnostic |
| 10. Reporting | One row per DGP x regime x method x term |
| 11. MC uncertainty | 50-rep pilot is directional; 500 reps required for release |

## Focused calibration amendment after the frozen control campaign

The frozen 500-replicate campaign at `3c0e413` used `smcfcs_numit = 20`,
`jomo_nburn = 1000`, and `jomo_nbetween = 100`. It left three standard-SMC
coverage misses while all added-bias and SE-ratio cells passed. Before changing
the public package, run paired, controls-only calibration on exactly the affected
DGP/regime cells:

1. binomial phylogeny: compare SMCFCS 20 versus 100 iterations;
2. mixed phylogeny and mixed auxiliary: compare jomo burn/spacing
   `(1000, 100)`, `(1000, 1000)`, and `(5000, 1000)`;
3. use 200 paired replicates per cell for setting selection, then 1,000 new,
   independent replicates per cell for the selected setting;
4. retain `M = 50`, the same DGP seeds, downstream models, and all bias/SE/
   coverage criteria;
5. select a stronger setting only when it does not worsen added bias or SE
   calibration and moves coverage toward 0.95 consistently across the two mixed
   regimes. Do not select on a single term crossing 0.925 by chance.

This amendment is prospective relative to the stronger-setting runs. It does
not reclassify the frozen campaign or authorize a release.

The final gate uses the unchanged 92.5--97.5% performance band and replaces
observed-MCSE filtering with the planned worst-band precision requirement
`sqrt(0.925 * 0.075 / n_planned) <= 0.01`. Thus final decisions require 1,000
replicates per cell. Observed MCSE and exact binomial confidence intervals remain
diagnostics. At `n = 500`, the former observed-MCSE rule implicitly required
coverage of about 94.7% or higher and contradicted the declared band.
