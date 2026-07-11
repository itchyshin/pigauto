# pigauto v0.10.0 multiple-imputation validation result

## Decision

**FAIL — block CRAN submission and redesign the multiple-imputation method.**

Neither conformal-width Normal draws nor Brownian-posterior/MC-dropout draws
passed the predeclared downstream fixed-effect gate. The package remains
experimental for multiple-imputation inference. Variance components, random
effects, correlations, BLUPs, and latent loadings remain unsupported for Rubin
pooling.

## Frozen evidence

- Release-candidate SHA: `800e0cb7be1006d8fe842e3e9a83624d096a59a4`
- Package version: `0.10.0`
- Platform: `x86_64-pc-linux-gnu`
- R: `4.5.3`
- torch: `0.17.0`
- Campaign: 3 DGPs x 2 regimes x 500 replicates = 3,000 fitted imputers
- Imputations: `M = 50` paired datasets per draw method and replicate
- Missingness: 30% of `x`, MAR through fully observed `y` and `z`
- Epoch cap: 1,000; 2.2% suggested more training, below the predeclared 5% escalation trigger
- Completion: 3,000/3,000 task files with `status = "success"`
- Provenance: 0 dirty records; 0 SHA mismatches

The three DGPs were Gaussian `lm`, binomial `glm`, and random-intercept
`lmer`. Each was evaluated in phylogeny-dominant and auxiliary-nonlinear
regimes. The complete-data oracle and complete-case analyses were retained as
comparators. The raw task and summary RDS files are retained on Totoro under
`~/pigauto_v010/pigauto/script/mi-validation-v010/results/full-800e0cb-epochs1000/`.

## Predeclared fixed-effect gate

Each method had 12 core cells: `x` and `z` for every DGP x regime combination.
Every cell required standardized bias no greater than 0.10, pooled-SE to
empirical-SD ratio from 0.90 to 1.10, coverage from 92.5% to 97.5%, coverage
Monte Carlo SE no greater than 1 percentage point, at least 95% successful
analyses/fits, and at least 99% finite valid Rubin quantities.

| Method | Cells passing all criteria | Standardized-bias range | SE-ratio range | Coverage range |
|---|---:|---:|---:|---:|
| Conformal-width Normal | 0/12 | 0.051–2.054 | 0.670–1.053 | 35.4%–93.6% |
| Brownian/MC-dropout | 0/12 | 0.077–1.320 | 0.796–1.037 | 52.2%–93.6% |

Only 1/12 cells per method met the bias threshold. Conformal met the SE-ratio
criterion in 6/12 cells and the coverage criterion in 3/12. MC-dropout met
those criteria in 7/12 and 4/12 cells, respectively. Success and finite-valid
rates were 100%, so numerical fitting failure does not explain the result.

The bias changed by DGP and regime. For example, conformal `x` bias was
`-0.354` in auxiliary `glm`, while MC-dropout `x` bias was `-0.176` in
phylogeny-dominant `glm`. This pattern, together with undercoverage, indicates
that the problem is not an interval-multiplier adjustment alone: the completed
data do not preserve the downstream fixed-effect relationship sufficiently.

## Variance-component diagnostics

These estimates were not Rubin-pooled and do not create a support claim.
Nevertheless, behavior was stable:

- 0 boundary fits and 0 singular fits in every mixed-model method/regime cell;
- residual-variance relative bias from 0.9% to 5.3% for the two MI methods;
- random-intercept-variance relative bias from -1.0% to -3.3%;
- maximum added absolute relative bias over the oracle: 4.5 percentage points;
- no predeclared variance-component diagnostic flag fired.

Thus the release blocker is downstream fixed-effect calibration, not obvious
variance-component instability.

## PMM redesign pilot

A harness-only stochastic PMM candidate was tested on 10 replicates per DGP x
regime cell (60 paired fitted imputers) at SHA `568e0e2c4ba994f0dd21d4aa17502a800dc4f005`.
It used the same Brownian/MC-dropout predictions and selected among the five
nearest observed donors. This was directional evidence only.

PMM did not repair the problem: standardized bias for `x` remained above 0.10
in all six cells (approximately 0.43–1.70), and coverage remained unstable.
PMM was therefore not promoted into `multi_impute()`.

## Consequence and next decision

Do not tag v0.10.0 and do not submit it to CRAN under the current release gate.
The package documentation should continue to call MI fixed-effect inference
experimental. The next method-development step requires an explicit scope
decision: either build a substantive-model-compatible imputation strategy and
rerun the full gate, or defer inferential MI from the first CRAN release and
revise the predeclared release scope. Repeating the same campaign with more
epochs or PMM is not supported by these results.
