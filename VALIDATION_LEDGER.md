# Validation ledger

Evidence tiers follow the package validation-harness convention:
`validated` means known-DGP recovery has passed; `experimental` means the path is
implemented but recovery is pending; `assumed` marks an unverified
parameterization.

| Claim | Evidence | Tier | Date |
|---|---|---|---|
| Rubin pooling for downstream `lm`/`glm` fixed effects has calibrated bias, SE, and 95% coverage after pigauto MI | Full 500-replicate-per-cell campaign at `800e0cb` failed every conformal and MC-dropout core cell; see `script/mi-validation-v010/FULL_RESULTS_2026-07-11.md` | experimental | 2026-07-11 |
| Rubin pooling for downstream random-intercept `lmer` fixed effects has calibrated bias, SE, and 95% coverage after pigauto MI | Full 500-replicate-per-cell campaign at `800e0cb` failed the joint fixed-effect gate; variance components were diagnostic only and showed no boundary/singularity flags | experimental | 2026-07-11 |
| Conformal-width Normal draws form an inferentially valid imputation distribution | Full paired campaign: 0/12 core cells passed; coverage 35.4% to 93.6% | experimental | 2026-07-11 |
| Brownian-posterior/MC-dropout draws form an inferentially valid imputation distribution | Full paired campaign: 0/12 core cells passed; coverage 52.2% to 93.6% | experimental | 2026-07-11 |
| Stochastic predictive mean matching repairs the failed downstream fixed-effect behavior | Exploratory 60-task paired pilot on `568e0e2`: all six `x` cells retained standardized bias above 0.10; not promoted to the public API | experimental | 2026-07-11 |
| The frozen downstream MI gate is attainable with proper outcome-compatible draws | Controls-only 3,000-task campaign at `3c0e413`: exact DGP conditional draws had negligible added bias and 11/12 coverage cells in range; one cell was 92.4%, so this remains a near-pass rather than validated | experimental | 2026-07-11 |
| Off-the-shelf substantive-compatible MI is calibrated for all frozen `lm`/logit `glm`/`lmer` cells | `smcfcs`/`jomo.smc` campaign at `3c0e413`: all added-bias and SE-ratio cells passed, but 3/12 coverage cells were 92.0–92.2%; see `script/mi-validation-v010/CONTROL_RESULTS_2026-07-11.md` | experimental | 2026-07-11 |
| Random-effect variances, correlations, and BLUPs can be pooled by `pool_mi()` | Explicitly outside v0.10.0 scope; variance components are diagnostic only | experimental | 2026-07-10 |

The v0.10.0 CRAN gate is blocked. No inferential claim above should be promoted
to `validated`, and no CRAN submission should be made, until a redesigned method
passes every pre-specified core cell on a new frozen package SHA.
