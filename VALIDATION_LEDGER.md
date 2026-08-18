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
| Analysis-aware MI is an attainable target for one incomplete continuous covariate under MAR with Gaussian `lm`, binomial-logit `glm`, and one-random-intercept Gaussian `lmer` | Warning-free 6,000-task controls-only campaign at clean SHA `430b2c9` passed all 24 method-by-term cells using proper Bayesian normal MI, SMCFCS, and jomo; this validates the target, not the package implementation | experimental | 2026-07-11 |
| `multi_impute_analysis()` package implementation passes the analysis-aware fixed-effect gate | Warning-free 6,000-task package campaign at clean SHA `2e3809d` passed all 24 method-by-term cells; added bias 0.00086-0.05042 oracle SD, SE ratio 0.942-1.030, coverage 93.9%-96.3%, 100% finite results | validated | 2026-07-11 |
| Random-effect variances, correlations, and BLUPs can be pooled by `pool_mi()` | Explicitly outside v0.10.0 scope; variance components are diagnostic only | experimental | 2026-07-10 |
| One recommendation from `suggest_next_observation()` improves realised follow-up imputation over random acquisition | Audited Stage-A recovery campaign: 8,000 receipts / 32,000 policy rows, zero failures, and correct treatment/provenance receipts. In the continuous BM lambda=1 one-step cell-level regime, active-minus-random normalized MSE is negative at n=100 (-0.01223, paired 95% MC interval -0.02190 to -0.00256) and n=300 (-0.00431, -0.00715 to -0.00147). Both binary lambda=1 intervals include zero, so the broad cross-family claim is not validated. | experimental | 2026-08-18 |

The warning-free package campaign validates the narrow fixed-effect backend at
its documented boundary. The v0.10.0 CRAN gate now depends on package checks
and cross-platform release engineering. Legacy conformal-width,
Brownian/MC-dropout, PMM, and posterior-tree draws remain unsupported for
downstream inference.
