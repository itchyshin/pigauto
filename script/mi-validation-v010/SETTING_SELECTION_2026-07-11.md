# SMC setting-selection result

## Decision

Retain the economical frozen settings for independent confirmation:

- `smcfcs_numit = 20`;
- `jomo_nburn = 1000`;
- `jomo_nbetween = 100`.

Longer iteration schedules did not materially improve bias, SE calibration, or
coverage on paired DGP replicates. This is a calibration-setting decision, not a
release pass.

## Frozen experiment

- Harness SHA: `70fa938b102e9bfcd00a548c5e05958f92b5906b`
- Replicates: 200 per affected DGP/regime cell, identical seeds across settings
- Imputations: `M = 50`
- Completion: 1,600/1,600 tasks, zero errors, zero captured engine warnings
- No pigauto training was run
- Totoro output: `script/mi-validation-v010/results/calibration-70fa938/`

## SMCFCS

For binomial phylogeny, `numit = 20`, `50`, and `100` gave effectively the same
paired result. Coverage was 93.5%, 93.0%, and 93.0% for `x`, and 96.0%, 96.0%,
and 96.5% for `z`. The `z` cell that was 92.2% in the earlier 500-replicate
campaign was therefore not repaired by longer chains; its first 200 paired
replicates were already at 96.0%. Median runtime rose to about 12 seconds at 50
iterations and 24 seconds at 100 iterations. Retain 20.

## jomo

The paired schedules `(nburn, nbetween) = (1000, 100)`, `(5000, 100)`,
`(1000, 1000)`, and `(5000, 1000)` produced nearly identical results in both
mixed-model regimes. For example, auxiliary `z` coverage was 90.5%, 91.0%,
91.0%, and 91.0%, with SE ratios 0.8884--0.8888. Phylogeny `x` coverage was
93.5%, 94.5%, 93.5%, and 93.5%. Longer spacing raised median runtime from about
3 seconds to 16--17 seconds without improving calibration. Retain `(1000, 100)`.

The low auxiliary-`z` SE ratio in this 200-replicate subset is substantially
lower than in the earlier 500-replicate campaign, consistent with appreciable
outer Monte Carlo variation. Burn-in and spacing are not the explanation.

## Prospective confirmation

Run the selected defaults on independent replicate IDs 501--1500:

- binomial phylogeny: 1,000 replicates;
- mixed phylogeny and mixed auxiliary: 1,000 replicates each;
- exact conditional and standard SMC controls paired within every replicate;
- unchanged bias, SE-ratio, coverage, fit, and finite-result thresholds;
- planned worst-band coverage MCSE below 0.01 by design.

The restricted confirmation is calibration-only and cannot issue a package
release pass. It determines whether the three earlier coverage misses persist
before any public analysis-aware backend is implemented.

## Rejection-limit amendment after the full-grid confirmation

The independent 1,000-replicate full grid passed every statistical control cell,
but SMCFCS emitted rejection-limit warnings in 106/2,000 `lm` tasks at
`rjlimit = 10000`. Some warnings reported hundreds of exhausted proposals.
Warnings are now captured and a warning-free campaign is required; the
statistical pass is not accepted while these remain.

`smcfcs_rjlimit` is therefore added to the manifest and task provenance. Test
larger limits on the two worst seeds (`lm` phylogeny replicate 539 and auxiliary
replicate 793) before rerunning the single-level grid. A larger limit is accepted
only if both tasks are warning-free and preserve finite imputations, observed
values, and row order. This amendment changes a numerical rejection-sampling
safeguard, not the imputation model or performance thresholds.

## Superseding Gaussian-engine decision

The rejection-limit experiment did not repair the Gaussian `lm` branch:
`rjlimit = 100000` reduced but did not remove warnings on the two worst seeds,
and `rjlimit = 1000000` was not economically viable. The next prospective
confirmation therefore replaces SMCFCS only for the frozen Gaussian `lm` cell
with proper Bayesian normal-regression MI using `x ~ y + z + I(z^2)`. SMCFCS
remains the logit `glm` engine and jomo remains the `lmer` engine. Historical
results above retain their original engine labels.
