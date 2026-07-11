# Final MI positive-control result

## Decision

**CONTROL PASS ONLY.** The frozen downstream gate is attainable with the
analysis-compatible positive controls. This does not validate pigauto's
conformal, MC-dropout, or PMM draws, does not validate random-effect quantities,
and does not by itself authorize CRAN submission.

## Frozen evidence

- Source SHA: `430b2c9415f9fd097f26835687a61b77b6e33ce1`
- Platform: `x86_64-pc-linux-gnu` on Totoro
- Campaign: 3 DGPs x 2 regimes x 1,000 replicates = 6,000 tasks
- Imputations: `M = 50`
- Completion: 6,000/6,000 successful; all pooled results finite
- Captured engine warnings: zero
- Gaussian `lm`: proper Bayesian normal-regression MI from
  `x ~ y + z + I(z^2)`
- Logit `glm`: `smcfcs` 2.0.2, `numit = 20`, `rjlimit = 100000`
- Random-intercept `lmer`: `jomo` 2.7-6, `nburn = 1000`,
  `nbetween = 100`
- Raw output:
  `script/mi-validation-v010/results/controls-430b2c9-rj100k-reps1000/`
  on Totoro (ignored by git)

The `rjlimit = 100000` setting was selected prospectively after the otherwise
passing 10,000-limit campaign emitted one logit rejection-limit warning. The
exact warning seed completed without warnings at 100,000, and the entire new
campaign was rerun from a clean checkout at the same source SHA.

## Fixed-effect results

All 24 method-by-term cells passed the preregistered added-bias, SE-ratio,
coverage, finite-result, fit-success, and planned-MCSE gates.

| Method | Passing cells | Added-bias range | SE-ratio range | Coverage range | Finite-valid rate |
|---|---:|---:|---:|---:|---:|
| Exact known-DGP conditional | 12/12 | 0.000002-0.005397 | 0.974-1.038 | 93.6%-96.4% | 100% |
| Analysis-compatible standard control | 12/12 | 0.000860-0.044771 | 0.942-1.030 | 93.9%-96.1% | 100% |

The planned worst-band coverage Monte Carlo SE was 0.00833 in every cell,
below the 0.01 design threshold.

## Claim boundary

This campaign supports implementation and package-level testing of a narrow
analysis-aware backend for fixed effects under the validated designs: one
partially observed continuous covariate under MAR, fully observed outcome and
auxiliaries, Gaussian `lm`, logit `glm`, or Gaussian random-intercept `lmer`.
It does not support claims for multiple incomplete variables, missing discrete
covariates, MNAR, nonlinear or interaction analyses involving the incomplete
covariate, random slopes, variance components, correlations, BLUPs, latent
loadings, posterior-tree propagation, or arbitrary downstream model classes.

The public backend must still pass package unit/integration tests and the CRAN
engineering gate. Until then the release decision remains blocked.
