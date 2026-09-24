# Note for the lambda lane: plug-in covariance shrinks cross-trait correlation under missingness

2026-09-24. From the MI-GLS sweep on `arc/mi-gls-attenuation` (`docs/dev-log/mi-gls/results.md`,
finding 3). Written here because `R/joint_mvn_solver.R` belongs to the lambda lane; this branch does
not edit it.

## What was measured

- Setting: two continuous traits under bivariate Brownian motion, true correlation 0.7, n = 300,
  30% of both traits missing, six trees.
- Result: the plug-in cross-trait correlation from the joint-MVN solver was 0.46, against 0.67 from the
  same data before masking (true value 0.7). Raising `max_iter` to 50 gave the same value, so this is not a convergence artefact.
- Consequence: conditional draws that used this covariance biased a downstream PGLS slope by -0.08
  to -0.11 in the both-missing regimes. With an EM covariance estimate
  (`.dcb_sigma_em()` in `R/draws_conditional.R`) the bias was -0.009 to +0.053.

## Why it matters to the lambda lane

The same estimator feeds pigauto's joint baseline, so the baseline's own predictions borrow
cross-trait information with a correlation that is too weak whenever several traits are missing.
The PR #187 lambda work does not change this; it changes the phylogenetic scale, not the
cross-trait covariance estimator.

## Not covered

- More than two traits.
- Missingness other than 30%.
- lambda below 1 in this specific check.

## A check for the lambda lane

On the same six trees, compare the solver's cross-trait correlation with EM's across 0% to 50% missing.
If the shrinkage grows with missingness, replace the plug-in step with EM or with the full-data
likelihood.

The posterior sampler on `arc/mi-posterior` estimates Sigma_P and Sigma_E by Gibbs, so it does not
inherit this problem; its recovery gate (G3) will show that directly.
