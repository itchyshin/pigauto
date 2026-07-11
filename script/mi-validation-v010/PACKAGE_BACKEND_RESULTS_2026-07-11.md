# Package-level analysis-aware MI validation

## Decision

**PASS within the documented fixed-effect scope.** The public
`multi_impute_analysis()` implementation passed every predeclared core cell.
This validates the narrow backend; it does not validate legacy conformal,
MC-dropout, PMM, or posterior-tree draws, and it makes no random-effect pooling
claim.

## Frozen evidence

- Source SHA: `2e3809da0a0cbc42032688406224fb0492a21514`
- Platform: `x86_64-pc-linux-gnu` on Totoro
- Campaign: 3 DGPs x 2 regimes x 1,000 replicates = 6,000 tasks
- Imputations: `M = 50`
- Completion: 6,000/6,000 successful; all pooled results finite
- Captured engine warnings: zero; engine warnings are fatal in the public API
- Gaussian `lm`: public proper Bayesian normal-regression route
- Logit `glm`: public `smcfcs` 2.0.2 route, `numit = 20`,
  `rjlimit = 100000`
- Random-intercept `lmer`: public `jomo` 2.7-6 route, `nburn = 1000`,
  `nbetween = 100`
- Raw output:
  `script/mi-validation-v010/results/package-2e3809d-reps1000/` on Totoro
  (ignored by git)

## Results

All 24 method-by-term cells passed the added-bias, SE-ratio, coverage,
finite-result, fit-success, and planned-MCSE gates.

| Method | Passing cells | Added-bias range | SE-ratio range | Coverage range | Finite-valid rate |
|---|---:|---:|---:|---:|---:|
| Exact known-DGP conditional | 12/12 | 0.000002-0.005397 | 0.974-1.038 | 93.6%-96.4% | 100% |
| Public analysis-aware backend | 12/12 | 0.000860-0.050415 | 0.942-1.030 | 93.9%-96.3% | 100% |

The planned worst-band coverage Monte Carlo SE was 0.00833 in every cell,
below the 0.01 design threshold.

## Validated boundary

The evidence covers fixed effects with one partially observed continuous
covariate under MAR, a fully observed outcome and auxiliaries, and one of:

- an additive Gaussian `lm`;
- an additive binomial-logit `glm`;
- an additive Gaussian `lmer` with one random intercept.

It does not cover multiple incomplete variables, missing discrete covariates,
MNAR, transformations or interactions involving the incomplete covariate,
random slopes, variance components, correlations, BLUPs, latent loadings,
posterior-tree propagation, or arbitrary downstream model classes.

The statistical gate is complete. CRAN submission still requires the final
clean package, pkgdown, cross-platform, and win-builder engineering checks.
