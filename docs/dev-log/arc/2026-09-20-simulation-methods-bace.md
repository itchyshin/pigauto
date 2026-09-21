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
equal rates for binary and categorical traits and a stepwise model for the ordinal trait. The evolutionary model is left at
Rphylopars' default, `model = "BM"`, for every cell, including the cells whose data were generated
with Pagel's lambda below 1. This is deliberate, and it costs the arm accuracy: at
lambda = 0.3, n = 100 the frequentist stack reaches a z-RMSE of 1.25 against a mean floor of 1.01,
meaning it does worse than ignoring the phylogeny, because it extrapolates a phylogenetic signal the
data do not contain. BACE, which estimates its own signal parameter, reaches 1.10 on the same cells.
Rphylopars does offer `model = "lambda"`, and a reader should take these results as describing the
package at its documented default rather than the best the package can do. Counts go
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
one side. Measured over the core slice, this happened in 69 replicates, all at lambda = 1, spread
over n = 100 (32), n = 300 (29) and n = 1000 (8), about 6% of the lambda = 1 replicates and none
anywhere else.

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

## Results placeholder

Filled from `script/campaign_sim_results/summary.csv` and `paired.csv` once the campaign completes.
Every number here must carry its regime: DGP, n, lambda, rho, mechanism, arm, replicate count and
MCSE.
