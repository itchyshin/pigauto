# Methods: simulation study comparing the frequentist stack with BACE

Draft for the BACE paper. Written to be pasted into the manuscript and edited, so it states every
setting a reader would need to reproduce the study. Structure follows ADEMP (Morris, White and
Crowther 2019, *Stat. Med.* 38: 2074) and the reporting items of Williams, Bates and colleagues
(2024, *MEE* 15: 1926). Numbers are filled in from the completed campaign; this file carries the
design, which was fixed before any result was seen.

Scope here is arms 1 and 2 only. The four-arm comparison that adds pigauto is reported separately.

## Aim

Does joint Bayesian imputation improve accuracy, probability calibration and interval coverage over
the standard frequentist stack, across trait types, phylogenetic signal, cross-trait correlation,
missingness mechanism and sample size?

The primary contrast was named before any result was inspected: **BACE against the frequentist stack
on z-RMSE and 95% interval coverage, with the interval score beside them, on the core slice, pooled
over trait types.** Everything else in this study is secondary and is labelled as such.

## Data-generating mechanism

One tree per replicate, one row per species, eight traits per species.

#### Tree

`ape::rcoal(n)` rescaled to unit height, so it is ultrametric. BACE requires an ultrametric
tree; `ape::rtree` is not used anywhere.

#### Latent traits

For a tree correlation matrix `V = cov2cor(vcv(tree))`, Pagel's lambda enters as

```
V_lambda = lambda * V + (1 - lambda) * I
```

with unit diagonal and no separate residual term, so the marginal variance of every latent is exactly
1 at every lambda. Cross-trait correlation enters through a K x K matrix `Sigma_rho` with unit
diagonal and off-diagonal `rho`. Latents are drawn as

```
L = t(chol(V_lambda)) %*% Z %*% chol(Sigma_rho),    Z an n x K matrix of standard normals
```

which gives `Var(L_ik) = 1`, `Cov(L_ik, L_il) = rho` and `Cov(L_ik, L_jk) = lambda * V_ij` exactly.
Under the Ornstein-Uhlenbeck sensitivity, `V` is replaced by the stationary OU tip correlation
`exp(-alpha * d_ij)` with `alpha = 2` and `d_ij` the patristic distance, standardised to unit
diagonal; on an ultrametric tree that expression is the stationary OU correlation, not an
approximation.

#### Observed traits

Two continuous traits are the latents themselves. A count is
`Poisson(exp(1.5 + 0.8 L))`. A proportion is `plogis(L + N(0, 0.3^2))`, whose extra noise sits
outside the lambda bookkeeping and is stated rather than absorbed. A binary trait, a four-level
ordinal and a three-level categorical are cut from their latents at **fixed population thresholds**
(0; the quartiles of the standard normal; its terciles). Fixed thresholds were chosen over sample
quantiles so that class balance varies between replicates, which is realistic and exercises the
rare-class regime where both castor and MCMCglmm can fail; sample quantiles would force exact balance
every replicate and leak information from the cells that are later masked.

#### The missingness driver

An eighth continuous trait, `d1`, is always observed and never masked. It
carries its own correlation of 0.35 to every other trait, independent of `rho`. That value is the
largest that keeps the correlation matrix positive definite with seven other traits, and the code
checks the smallest eigenvalue at run time. Without this, a replicate at `rho = 0` would have a
driver independent of everything it is supposed to predict, and the missing-at-random condition would
be missing-completely-at-random wearing another label.

`d1` is present in every cell of every mechanism, including those with no MAR component. That is worth stating
plainly because it changes what the rho factor measures. A fully observed trait correlated 0.35 with
all seven scored traits is available to any method that models traits jointly, so even at rho = 0 the
joint arms have cross-trait information to exploit. The rho factor should therefore be read as the
trait-to-trait correlation over and above that shared observed covariate, rather than as the
presence or absence of any cross-trait information at all. Every arm sees `d1` on equal terms, so
the comparison between arms is unaffected; what changes is the interpretation of the rho axis.

#### Missingness

Applied to the complete data, with the same mask given to every arm within a
replicate; masks are seeded from `seed + 1000`.

- MCAR at 10% and at 30%.
- MAR at 30%: `P(missing) = plogis(a + log(3) * z(d1))`, an odds ratio of 3 per standard deviation of
  the driver, with the intercept `a` solved by `uniroot` in each replicate so the expected fraction
  matches the target. The realised fraction is stored as a diagnostic.
- Phylogenetically biased at 30%, following Gendre and colleagues (2024): random subtrees of 5 to 15%
  of the tips are masked whole until the target fraction is reached, subject to at least five
  observed cells remaining in every column.

#### Factors

Sample size n in {100, 300, 1000}; lambda in {0.3, 0.7, 1.0}; rho in {0, 0.5};
missingness in the four mechanisms above; evolutionary model BM or OU. The **core slice** is BM with
MCAR 0.30 and no covariates, which is lambda x rho x n = 18 cells. The factorial adds 56 further
cells; n = 300 and lambda = 0.7 appear in the core slice only, to spend the compute on the contrasts
that separate the arms rather than on filling a grid. A covariate sensitivity runs the 18 core cells
with two environmental covariates.

#### Replicates

200 per cell for the frequentist stack, 100 for BACE, with BACE's seeds nested inside
the others (1 to 100 of 1 to 200) so every contrast between the two arms is paired on a common mask.

## Estimands and targets

For every masked cell the truth is stored with the replicate, so all targets are computed on masked
cells only.

- Point accuracy: z-RMSE on continuous-family traits, standardised by the training standard
  deviation; accuracy and macro-F1 on discrete traits. Macro-F1 averages only over classes present in
  the masked truth, and the number of classes scored is recorded: scoring an absent class as zero and
  averaging it in deflates the metric, and deflates it most in exactly the low-prevalence cells where
  the arms differ.
- Probability calibration: Brier score with its reliability component, and a calibration slope.
  Expected calibration error is reported as a secondary measure, pooled over masked cells across
  replicates within a cell and binned by equal mass; per-replicate ECE on a few dozen cells has a bias
  of roughly 0.29 that grows with how many bins an arm occupies, which would penalise the sharper arm.
- Interval coverage and width: Empirical 95% coverage and mean standardised width, with the
  interval score of Gneiting and Raftery (2007) as the single proper scoring rule that trades the two
  against each other.
- Cost and robustness: Wall time per fit, and the failure rate per cell. A fit that fails is
  scored at the mean or mode floor and its failure is reported beside every metric. No fit is dropped.

The three intervals are different objects, and the paper says so. Rphylopars supplies a
model-based plug-in variance for an unobserved tip; the check that it already includes the phenotypic
component was run empirically (five replicates at n = 300 with phenotypic error, giving 0.933
coverage against a nominal 0.95, where omitting that component would give roughly 0.70). BACE's
interval is a Bayesian posterior predictive interval taken as the 2.5 and 97.5 percentiles over its
final imputed datasets, each of which is one posterior predictive draw. The two are never described
as interchangeable, and coverage is never called better without the standardised width beside it.

## Methods compared

#### Arm 1, the frequentist stack

`Rphylopars::phylopars(model = "BM", phylo_correlated = TRUE,
pheno_correlated = TRUE, REML = TRUE)` fitted jointly to the continuous-family columns, with counts
and proportions transformed and back-transformed. `castor::hsp_mk_model` per discrete trait, with
equal rates for binary and categorical traits and a stepwise model for the ordinal trait. The evolutionary model is the one setting on which this arm is reported twice.
Rphylopars' documented default is `model = "BM"`, which pins the phylogenetic signal at its maximum
even in cells whose data were generated with Pagel's lambda below 1. Because a reader could
reasonably take either the default or the signal-estimating variant as the thing being compared,
both were run on every cell as separate arms: `freq` with `model = "BM"` and `freq_lambda` with
`model = "lambda"`, which estimates the signal from the data. The difference is large enough that it
changes the conclusion, and it is reported in full below rather than left as a caveat. Counts go
through `phylolm::phyloglm(method = "poisson_GEE")` rather than a log1p Gaussian fit, because
Rphylopars on log1p counts fell below the mean floor at n = 100 in the feasibility test; its interval
is a marginal Poisson band, the same for every masked cell of that trait, and the results table says
so. This stack carries no information between trait types.

#### Arm 2, BACE

`BACE::bace()` with `nitt = 50000`, `burnin = 10000`, `thin = 25`, `runs = 5`,
`n_final = 20`, one-versus-rest categorical handling, and `n_cores = 1`.

Two settings deserve a sentence each, because both were chosen from measurement rather than default.
`runs` is the number of initial imputation iterations that BACE's own `assess_convergence()` reads;
that function requires at least three, and a pilot at `runs = 3` failed to converge in both replicates
while `runs = 5` converged in both, at 6% more cost. `n_final` is the number of final imputation runs,
each of which refits the model per response, so it is a cost parameter and not a chain length; 20 was
chosen to give usable percentile intervals at a measured 2 h 46 m per replicate at n = 1000.

Convergence is not asserted. `skip_conv` is left at `TRUE`, so BACE assesses convergence but does not
retry, which keeps the cost of a cell bounded across thousands of cells; its verdict is recorded for
every fit and the convergence rate is reported as a result. Effective sample sizes over the fixed
effects are reported alongside, with the median rather than the minimum used as the summary, because
MCMCglmm's threshold and categorical models mix slowly by construction and a single badly-mixing
parameter is expected.

#### A failure mode both arms share

At lambda = 1 with fixed population thresholds, a discrete trait can come out monomorphic among the
observed cells: maximum phylogenetic signal plus a threshold cut can leave every observed species on
one side. Measured over the complete core slice, this happened in 77 replicates, all at lambda = 1,
spread over n = 100 (32), n = 300 (29) and n = 1000 (16), about 6% of the lambda = 1 replicates and
none anywhere else.

castor cannot fit an Mk model to a trait with one observed state, and reports so. That is a real
property of the frequentist stack and is reported as a failure, scored at the floor. pigauto's
GNN-on arm fails on exactly the same cells for an unrelated reason, a dimension error on the
collapsed one-hot encoding, which is a robustness defect in the package rather than a statement
about the method; it is recorded here and fixed separately, outside this study.

Because both arms fail on the same cells and both are scored at the floor, the paired contrast
between them is not distorted, but the absolute figures for both at lambda = 1 are pulled toward the
floor. Any lambda = 1 number for those two arms should be read with that in mind.

#### Reference

A mean or mode floor, so that an arm doing worse than ignoring the phylogeny entirely
is visible as such.

## Performance measures and Monte Carlo error

Every reported number carries a Monte Carlo standard error. For a mean over replicates it is
`sd / sqrt(n_sim)`; for coverage, `sqrt(p (1 - p) / n_sim)`; for pooled ECE, a bootstrap over
replicates. Because the arms share a mask within a replicate, the headline contrasts are reported as
paired differences with the MCSE of the per-replicate difference, which is materially tighter than
differencing two independently-summarised means: in an earlier campaign on these data the paired MCSE
was 0.02 where the unpaired was 0.04.

At 200 replicates the coverage MCSE is 1.5 percentage points, so a 5-point difference is three MCSE;
at BACE's 100 replicates it is 2.2 points.

## Reproducibility

Every replicate stores its latent matrix, the realised lambda, rho and missing fraction, its seeds,
the random number generator kind, the pigauto commit, the host and the wall time, together with
`sessionInfo()`. The generator is set to `L'Ecuyer-CMRG` so that a replicate is reproducible across
hosts; this is checked directly, by running overlapping cells on two machines and requiring the data,
the mask and the deterministic frequentist arm to agree.

Software: R 4.5 and 4.6; Rphylopars 0.3.10; castor 1.8.7; phylolm 2.6.5; MCMCglmm 2.36; BACE
0.0.0.9000. Runners and drivers are in `script/` of the pigauto repository
(`campaign_sim_cell.R`, `campaign_gnn_off_lib.R`, `campaign_sim_design.R`), and the aggregated results
with their per-cell seed and failure tables are committed beside them.

Compute ran on a 384-core shared server and on three Digital Research Alliance of Canada clusters,
one process per cell, with BLAS threads pinned to one.

## Results

All figures below come from the core slice: types_mixed, Brownian motion, MCAR at 0.30, no
covariates, lambda in {0.3, 0.7, 1.0}, rho in {0, 0.5}, n in {100, 300, 1000}. Rho is pooled because
it moved nothing of substance. The frequentist arms ran 200 replicates per cell and BACE ran 100.
Monte Carlo standard errors are in parentheses. Continuous-family figures pool the two continuous
traits, the count and the proportion; discrete-family figures pool the binary, ordinal and
three-level categorical traits.

### Accuracy on the continuous family, z-RMSE

| n | lambda | BACE | freq (BM) | freq_lambda | floor |
|---|---|---|---|---|---|
| 100 | 0.3 | 1.004 (0.0089) | 1.158 (0.0074) | 0.913 (0.0056) | 1.014 |
| 100 | 0.7 | 0.874 (0.0105) | 0.984 (0.0072) | 0.801 (0.0053) | 1.015 |
| 100 | 1.0 | 0.794 (0.0123) | 0.584 (0.0065) | 0.566 (0.0062) | 1.018 |
| 300 | 0.3 | 0.923 (0.0057) | 1.127 (0.0055) | 0.888 (0.0036) | 1.007 |
| 300 | 0.7 | 0.777 (0.0058) | 0.965 (0.0057) | 0.778 (0.0035) | 1.004 |
| 300 | 1.0 | 0.737 (0.0127) | 0.531 (0.0058) | 0.515 (0.0054) | 1.008 |
| 1000 | 0.3 | 0.892 (0.0036) | 1.106 (0.0043) | 0.878 (0.0023) | 1.006 |
| 1000 | 0.7 | 0.743 (0.0044) | 0.943 (0.0048) | 0.766 (0.0030) | 1.004 |
| 1000 | 1.0 | 0.659 (0.0132) | 0.475 (0.0048) | 0.464 (0.0044) | 1.002 |

The primary contrast, stated in advance, is BACE against the frequentist stack on z-RMSE. It has two
answers, and which one is reported depends entirely on how the frequentist stack is specified.

Against `freq`, the documented default, BACE wins wherever the phylogenetic signal is below its
maximum. Paired differences, negative favouring BACE: 0.141 (0.0094) at n = 100 and lambda = 0.3,
0.207 (0.0074) at n = 300, 0.205 (0.0045) at n = 1000, with smaller but same-signed differences at
lambda = 0.7. At lambda = 1 the sign reverses and BACE loses by 0.168 to 0.213, because Brownian
motion is then the true process and the default specification is correct.

Against `freq_lambda`, which estimates the same signal parameter BACE estimates, the advantage
disappears. At lambda = 0.3 freq_lambda is ahead by 0.10 at n = 100 and by 0.03 at n = 300 and
n = 1000. At lambda = 0.7 the two are within one to three MCSE of each other in both directions. At
lambda = 1 freq_lambda is ahead by 0.20 to 0.24.

The honest summary is that BACE's apparent advantage on continuous traits is an advantage over a
misspecified competitor rather than over the frequentist approach. Once both estimate the signal,
the two are close at low and moderate signal and BACE is behind at high signal.

Note also that `freq` at lambda = 0.3 sits above the mean floor at every n, between 1.106 and 1.158
against a floor near 1.01. Used at its default, the stack does worse than ignoring the phylogeny.
`freq_lambda` sits below the floor in every cell.

### Accuracy on the discrete family

| n | lambda | BACE | freq (BM) | freq_lambda | floor |
|---|---|---|---|---|---|
| 100 | 0.3 | 0.528 (0.0043) | 0.421 (0.0036) | 0.421 (0.0035) | 0.458 |
| 100 | 0.7 | 0.641 (0.0059) | 0.563 (0.0054) | 0.568 (0.0053) | 0.545 |
| 100 | 1.0 | 0.767 (0.0075) | 0.901 (0.0030) | 0.901 (0.0030) | 0.630 |
| 300 | 0.3 | 0.565 (0.0031) | 0.414 (0.0028) | 0.414 (0.0028) | 0.472 |
| 300 | 0.7 | 0.687 (0.0040) | 0.560 (0.0049) | 0.562 (0.0049) | 0.557 |
| 300 | 1.0 | 0.794 (0.0076) | 0.938 (0.0023) | 0.938 (0.0023) | 0.638 |
| 1000 | 0.3 | 0.587 (0.0025) | 0.422 (0.0025) | 0.422 (0.0025) | 0.481 |
| 1000 | 0.7 | 0.705 (0.0039) | 0.575 (0.0047) | 0.574 (0.0047) | 0.562 |
| 1000 | 1.0 | 0.840 (0.0076) | 0.964 (0.0019) | 0.964 (0.0019) | 0.636 |

This is where BACE has a real and substantial advantage, and it is the opposite pattern to the
continuous one. At lambda = 0.3 BACE is ahead by 10.7, 15.1 and 16.5 accuracy points at n = 100, 300
and 1000, and the frequentist stack is below the mode floor in every one of those cells. At
lambda = 0.7 BACE is ahead by 7 to 13 points. At lambda = 1 the frequentist stack is ahead by 12 to
14 points, because `castor`'s Mk model recovers a strongly conserved discrete trait very well.

The two frequentist arms are identical on discrete traits to within 0.5 of a point, as they must be:
the `model` argument changes only the Rphylopars fit on the continuous family, and the discrete path
through `castor` is shared. This serves as an internal check that the two arms differ only where
intended.

### Interval coverage and interval score

| n | lambda | BACE | freq (BM) | freq_lambda |
|---|---|---|---|---|
| 100 | 0.3 | 0.803 | 0.799 | 0.881 |
| 100 | 0.7 | 0.813 | 0.807 | 0.887 |
| 100 | 1.0 | 0.871 | 0.883 | 0.885 |
| 300 | 0.3 | 0.824 | 0.845 | 0.894 |
| 300 | 0.7 | 0.833 | 0.853 | 0.902 |
| 300 | 1.0 | 0.879 | 0.902 | 0.899 |
| 1000 | 0.3 | 0.837 | 0.867 | 0.897 |
| 1000 | 0.7 | 0.842 | 0.875 | 0.906 |
| 1000 | 1.0 | 0.885 | 0.908 | 0.902 |

Neither the Bayesian nor the frequentist arm reaches the nominal 0.95 anywhere in the core slice.
That is the second finding of this study, and it applies to both routes compared here. It is not a
property of phylogenetic imputation in general: in the companion four-arm study, pigauto's split
conformal intervals reach 0.95 to 0.96 at n = 300 and n = 1000, because they are calibrated on
held-out residuals rather than derived from the fitted model. Coverage improves with n and with
lambda for every arm, which is the signature of intervals that are too narrow because they condition
on a fitted model rather than integrating over it. `freq_lambda` is closest to nominal across the
board, between 0.881 and 0.906, and BACE is furthest away at low signal, at 0.803 where nominal is
0.95. At 200 replicates the coverage MCSE is 1.5 points and at BACE's 100 replicates it is 2.2
points, so these gaps are many MCSE wide.

Interval score, which penalises width and non-coverage together and is the one-number comparison,
lower being better:

| n | lambda | BACE | freq (BM) | freq_lambda |
|---|---|---|---|---|
| 100 | 0.3 | 6.77 | 9.86 | 5.46 |
| 100 | 0.7 | 5.66 | 8.20 | 4.77 |
| 100 | 1.0 | 2.69 | 3.63 | 3.44 |
| 300 | 0.3 | 5.87 | 8.95 | 5.07 |
| 300 | 0.7 | 4.72 | 7.60 | 4.39 |
| 300 | 1.0 | 1.98 | 3.27 | 2.99 |
| 1000 | 0.3 | 5.45 | 8.84 | 5.03 |
| 1000 | 0.7 | 4.40 | 7.54 | 4.41 |
| 1000 | 1.0 | 1.74 | 3.22 | 2.91 |

BACE's intervals are better than either frequentist arm's at lambda = 1 and are beaten by
`freq_lambda` at lambda = 0.3 and 0.7. The default `freq` intervals are worst everywhere.

### Failures and convergence

Failed fits are scored at the mean or mode floor and are never dropped, so every figure above
includes them. Counts are errored plus divergent (a finite but absurd value, defined in the factorial
section below), written as errored + divergent.

| arm | n = 100 | n = 300 | n = 1000 |
|---|---|---|---|
| BACE | 112 + 7 of 600 (19.8%) | 113 of 600 (18.8%) | 96 of 600 (16.0%) |
| freq | 32 + 4 of 1200 (3.0%) | 29 + 6 of 1200 (2.9%) | 16 of 1200 (1.3%) |
| freq_lambda | 32 + 1 of 1200 (2.8%) | 29 of 1200 (2.4%) | 16 of 1200 (1.3%) |

BACE fails on roughly one replicate in six, with "mixed model equations singular" the dominant
message, and the rate barely improves with sample size.

Convergence is a separate matter from failure, and it was promised above as a reported rate. Every fit
from the runner of 2026-09-20 onward stores BACE's own `assess_convergence()` verdict over its five
imputation iterations, its drift, and the effective sample sizes of the final fits; 1,774 successful
fits across the core slice and the factorial carry them. By BACE's own verdict, **15% of fits
converged at lambda = 0.3 in the core and 20% in the factorial, 33% at lambda = 0.7, and 81% and 67%
at lambda = 1**; by sample size, 23 to 28% at n = 100, 43% at n = 300 and 48 to 61% at n = 1000. Median
effective sample size over the fixed effects is comfortable throughout, 572 at n = 100 rising to 1,600
at n = 1000, with a tenth percentile of 66 to 315. So the chains mix; what BACE's autocorrelation,
percent-change, trend and Geweke criteria object to is drift of the imputed values across the five
sequential iterations, most often where the phylogenetic signal is weak. Because `skip_conv = TRUE`, the
imputations scored in every table are the ones BACE produced regardless of that verdict. A user who
followed the package's retry path would pay an unbounded multiple of the 3-hour n = 1000 fit for it; a
user who did not would be using imputations BACE itself flags as unconverged in most low-signal cells.
Either way it is a cost of the Bayesian route that belongs beside its accuracy figures. The frequentist failures are the
monomorphic-discrete case described above and are confined to lambda = 1. This difference in
robustness is a practical cost of the Bayesian route that a user will meet, and it is reported here
rather than hidden by dropping the failed cells.

## Factorial extension

The factorial varies what the core slice held fixed: the evolutionary model (Brownian motion or
Ornstein-Uhlenbeck with alpha = 2), the missingness mechanism (MCAR at 10% and 30%, MAR at 30%
driven by an always-observed trait, and clade-biased at 30%), lambda in {0.3, 1.0}, rho in {0, 0.5}
and n in {100, 1000}, minus the eight cells already in the core: 56 cells. The frequentist arms ran
200 replicates per cell. BACE ran 100 at n = 100 and 30 at n = 1000, a reduction Shinichi approved
on 2026-09-21 after a successful n = 1000 fit measured at about 3 hours and over 16 GB, with the
MCSE widening by about 1.8x at n = 1000 and reported as such. Strata below give equal weight to
each cell and pool over rho, which again moved nothing of substance. Every figure comes from the
committed aggregate (`script/campaign_sim_results/summary.csv`), the same file the pkgdown article
renders from.

Two things were added to the scoring rule for this stage, and both are reported rather than
absorbed. First, a fit that returns a finite but absurd value is treated as a failure that did not
throw: an arm is divergent on a replicate if its z-RMSE exceeds three times the mean floor or its
interval score exceeds 1,000, where the sane range across every arm and cell is 2 to 12. A divergent
replicate is scored at the floor, like an errored one, and its interval is dropped. This caught 15
`freq_lambda`, 34 Rphylopars-solver and 2 BACE replicates out of 11,200, values up to 8.6 x 10^18,
any one of which would otherwise dominate a stratum mean. Second, `freq_lambda` acquired a failure
mode of its own: 51 replicates, all in the OU, lambda = 1, MCAR 10%, n = 1000 cell, where Rphylopars'
singular-solve fallback ends in a type error. These are floored and counted.

### z-RMSE on the continuous family

| evolutionary model | lambda | BACE | freq (BM) | freq_lambda | floor |
|---|---|---|---|---|---|
| Brownian motion | 0.3 | 1.027 | 1.125 | 0.924 | 1.043 |
| Brownian motion | 1.0 | 0.872 | 0.579 | 0.572 | 1.070 |
| Ornstein-Uhlenbeck | 0.3 | 0.989 | 1.106 | 0.912 | 1.035 |
| Ornstein-Uhlenbeck | 1.0 | 0.608 | 0.568 | 0.562 | 1.051 |

| mechanism | BACE | freq (BM) | freq_lambda | floor |
|---|---|---|---|---|
| MCAR 10% | 0.725 | 0.783 | 0.691 | 1.000 |
| MCAR 30% | 0.743 | 0.814 | 0.697 | 1.009 |
| MAR 30% | 0.938 | 0.895 | 0.774 | 1.113 |
| clade-biased 30% | 0.986 | 0.867 | 0.783 | 1.054 |

The core-slice conclusion holds in every stratum. Against the BM default BACE wins at lambda = 0.3
and loses at lambda = 1; against the signal-estimating specification it loses everywhere on the
continuous family. Paired on BACE's own replicates, `freq_lambda` is ahead by 0.123 (MCSE 0.0085)
under Brownian motion at lambda = 0.3, by 0.095 (0.0058) under Ornstein-Uhlenbeck at lambda = 0.3,
and by 0.258 (0.0111) and 0.046 (0.0070) at lambda = 1. Ornstein-Uhlenbeck changes the magnitudes
and not the signs, so this is about estimating the signal, not about which process generated the
data.

The mechanism table carries the one new result. Under MAR and clade-biased missingness BACE falls to
0.938 and 0.986, close to the floors of 1.113 and 1.054, while the frequentist arms are unaffected in
their ranking. Partitioned by cause, as it should be before it is read: at lambda = 0.3 BACE fails
on no replicate under any mechanism, and its clade-biased figure of 1.181 at n = 100 and 1.140 at
n = 1000 under Brownian motion is a genuine accuracy result, worse than predicting the mean. At
lambda = 1 the picture inverts and becomes a failure result: BACE errors with "mixed model equations
singular" on 35% of MCAR 10% replicates, 39% of MAR and 39% of clade-biased ones, rising to 53% and
86% in the Brownian clade-biased cells at n = 100 and n = 1000. On the replicates where it does fit
there, it reaches 0.816 and 0.545, behind the frequentist arms but not collapsed. The two halves
belong side by side: at low signal BACE is inaccurate under phylogenetically biased missingness, and
at high signal it is fragile under it.

### Accuracy on the discrete family

| evolutionary model | lambda | BACE | freq (BM) | freq_lambda | floor |
|---|---|---|---|---|---|
| Brownian motion | 0.3 | 0.552 | 0.418 | 0.418 | 0.461 |
| Brownian motion | 1.0 | 0.752 | 0.910 | 0.908 | 0.605 |
| Ornstein-Uhlenbeck | 0.3 | 0.544 | 0.406 | 0.406 | 0.436 |
| Ornstein-Uhlenbeck | 1.0 | 0.824 | 0.875 | 0.872 | 0.541 |

BACE's discrete advantage at low signal survives the factorial intact: 13 to 14 points over the
frequentist stack at lambda = 0.3 under both processes, with the frequentist stack below the mode
floor in every one of those cells. At lambda = 1 the frequentist stack is ahead by 16 points under
Brownian motion and 5 under Ornstein-Uhlenbeck. Across mechanisms BACE's lead is 3 to 5 points under
MCAR and vanishes under clade-biased missingness, where all three arms sit at 0.62.

### Coverage and interval score

| mechanism | BACE | freq (BM) | freq_lambda |
|---|---|---|---|
| MCAR 10% | 0.874 | 0.869 | 0.890 |
| MCAR 30% | 0.842 | 0.856 | 0.884 |
| MAR 30% | 0.846 | 0.867 | 0.886 |
| clade-biased 30% | 0.836 | 0.915 | 0.894 |

Neither route reaches nominal 0.95 anywhere in the factorial either. The model-based Rphylopars
interval is the one construction that improves under clade-biased missingness, to 0.915 at the
default and 0.894 with lambda estimated, while BACE's posterior predictive interval is at its worst
there, 0.836. Interval scores follow: BACE 5.70 against `freq_lambda` 4.82 under clade-biased
missingness, and BACE ahead only under MCAR 10%, 3.47 against 6.75, where `freq_lambda`'s score is
inflated by its own type-error failures in one cell.

### Failures across the factorial

| arm | MCAR 10% | MCAR 30% | MAR 30% | clade 30% |
|---|---|---|---|---|
| BACE | 17.6% | 8.7% | 20.1% | 20.3% |
| freq | 1.6% | 0.3% | 1.7% | 1.6% |
| freq_lambda | 3.6% | 0.3% | 1.6% | 1.7% |

Percentages are of replicates, errored or divergent, both scored at the floor. BACE's failures are
confined to lambda = 1 and are the fixed-threshold monomorphic regime described above; every BACE
failure at lambda = 0.3 count is zero. The frequentist failures are the same monomorphic `castor`
case, plus `freq_lambda`'s Rphylopars type error in one cell.

### Covariate sensitivity

The 18 core cells were re-run with two continuous covariates supplied to every arm that can take them
(ncov = 2; 3,599 of 3,600 replicates at the time of writing, the last one computing). The frequentist
stack receives them through `phylolm` and `phyloglm`; pigauto's GNN receives them through its covariate
input; pigauto with the GNN off ignores them by construction, and its change is exactly zero, which is
the internal check that the two runs share their seeds. Paired on identical replicates, the covariates
improve the frequentist stack's z-RMSE by 0.051 (MCSE 0.0022) at lambda = 0.3, 0.024 (0.0023) at
lambda = 0.7 and nothing at lambda = 1, with `freq_lambda` gaining 0.032, 0.012 and nothing; coverage
moves by less than 0.002 everywhere. No ranking in the core slice changes. BACE was not part of this
run because its cost at n = 1000 put it outside the approved budget, so nothing here compares BACE with
or without covariates.

### What this study does not cover

These results describe one tree shape, one set of seven traits, and two synthetic covariates at most. Nothing here speaks to trees larger
than 1000 tips, to multiple observations per species, or to missingness that depends on the missing
value itself.
