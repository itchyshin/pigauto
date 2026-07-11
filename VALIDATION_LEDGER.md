# Validation ledger

Evidence tiers follow the package validation-harness convention:
`validated` means known-DGP recovery has passed; `experimental` means the path is
implemented but recovery is pending; `assumed` marks an unverified
parameterization.

| Claim | Evidence | Tier | Date |
|---|---|---|---|
| Rubin pooling for downstream `lm`/`glm` fixed effects has calibrated bias, SE, and 95% coverage after pigauto MI | `script/mi-validation-v010/`; full 500-replicate campaign pending | experimental | 2026-07-10 |
| Rubin pooling for downstream random-intercept `lmer` fixed effects has calibrated bias, SE, and 95% coverage after pigauto MI | `script/mi-validation-v010/`; full 500-replicate campaign pending | experimental | 2026-07-10 |
| Conformal-width Normal draws form an inferentially valid imputation distribution | Paired conformal versus Brownian/MC-dropout campaign specified but not yet run | experimental | 2026-07-10 |
| Brownian-posterior/MC-dropout draws form an inferentially valid imputation distribution | Paired conformal versus Brownian/MC-dropout campaign specified but not yet run | experimental | 2026-07-10 |
| Random-effect variances, correlations, and BLUPs can be pooled by `pool_mi()` | Explicitly outside v0.10.0 scope; variance components are diagnostic only | experimental | 2026-07-10 |

No inferential claim above should be promoted to `validated` until the frozen
package SHA passes every pre-specified core cell and the evidence artifact records
that SHA, platform, package versions, manifest, failures, and summary.
