# MI positive-control validation result

## Decision

**NEAR-PASS, NOT A RELEASE PASS — the gate is broadly attainable, but the
standard substantive-compatible comparator still has three undercoverage
cells.**

This campaign validates the new controls-only harness and narrows the redesign.
It does not validate pigauto's conformal, MC-dropout, or PMM draw methods, and it
does not authorize a CRAN submission.

## Frozen evidence

- Harness SHA: `3c0e41317e1734b36c4e7299c16e5a9a6836b6bd`
- Platform: `x86_64-pc-linux-gnu` on Totoro
- Campaign: 3 DGPs x 2 regimes x 500 replicates = 3,000 tasks
- Imputations: `M = 50`
- Completion: 3,000/3,000 successful; zero error logs
- Methods: exact known-DGP conditional draws and standard SMC
- Standard SMC engines: `smcfcs` for `lm`/logit `glm`; `jomo.smc` for `lmer`
- Raw output: `~/pigauto_v010/pigauto/script/mi-validation-v010/results/controls-3c0e413-reps500/`

The controls-only path skipped torch and pigauto training. The full campaign
completed with 96 workers, so no further neural campaign is justified until the
remaining inferential design issue is resolved.

## Gate amendment forced by the controls

The earlier gate used absolute truth-based standardized bias. The frozen
500-replicate complete-data logit oracle itself has standardized bias `0.134`
for auxiliary-regime `x`, despite 95.0% coverage. No imputation method can
remove inherent finite-sample bias in the downstream estimator without changing
the estimand. The harness therefore now gates the paired imputation-added bias:

`abs(mean(beta_MI - beta_complete)) / SD(beta_complete) <= 0.10`.

Absolute truth-based bias remains reported. Coverage and SE calibration remain
truth-based. This is a correction to an unattainable criterion, not a relaxation
of the imputation requirement.

The rule requiring observed coverage MCSE no greater than `0.01` is also
internally inconsistent with allowing coverage as low as `0.925` at 500
replicates: the corresponding MCSE is about `0.0118`. The current summary still
applies the rule and therefore reports an overall fail. A future preregistered
gate must either increase the replicate count or treat simulation precision as
a design property rather than a method-performance criterion.

## Results

All 24 method-by-term cells had 100% successful pooled analyses and finite Rubin
quantities.

| Method | Added-bias range | SE-ratio range | Coverage range | Coverage cells inside 92.5–97.5% |
|---|---:|---:|---:|---:|
| Exact DGP conditional | 0.0001–0.0086 oracle SD | 0.943–1.054 | 92.4%–96.4% | 11/12 |
| Standard SMC | 0.00002–0.0425 oracle SD | 0.914–1.036 | 92.0%–96.4% | 9/12 |

The exact conditional control missed only mixed auxiliary `z`, at 92.4%
coverage—one covered replicate below the 92.5% boundary. Standard SMC missed:

- mixed phylogeny `x`: 92.0%;
- binomial phylogeny `z`: 92.2%;
- mixed auxiliary `z`: 92.2%.

Thus outcome compatibility and parameter draws remove the large bias failure,
but the off-the-shelf comparator is not yet fully calibrated under the frozen
gate. The `jomo.smc` run used a 1,000-iteration burn-in and 100 iterations
between imputations; convergence/mixing and the single-level binomial miss need
focused investigation before a public backend is designed.

## Consequence

Merge the harness and evidence because they are correct and useful. Keep CRAN
blocked. The next economical slice is a focused comparator-calibration study,
not another pigauto model campaign: check longer `jomo.smc` spacing/convergence,
audit the binomial SMCFCS cell, and preregister a mathematically coherent MCSE
rule. Only after the standard control passes should pigauto expose an
analysis-aware inferential backend.
