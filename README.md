# pigauto: Phylogenetic Imputation via Graph AUTO-encoders <img src="https://raw.githubusercontent.com/itchyshin/pigauto/main/man/figures/logo.png" align="right" height="139" alt="pigauto logo"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/itchyshin/pigauto/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/itchyshin/pigauto/actions/workflows/R-CMD-check.yaml)
[![Tests](https://img.shields.io/github/actions/workflow/status/itchyshin/pigauto/test-coverage.yaml?branch=main&label=tests&logo=github)](https://github.com/itchyshin/pigauto/actions/workflows/test-coverage.yaml)
[![Coverage](https://img.shields.io/codecov/c/github/itchyshin/pigauto?label=coverage&logo=codecov)](https://app.codecov.io/gh/itchyshin/pigauto)
[![pkgdown](https://github.com/itchyshin/pigauto/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/itchyshin/pigauto/actions/workflows/pkgdown.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

> **Warning:** pigauto is experimental; use at your own risk. It needs
> further validation.

pigauto fills missing species trait values using a phylogenetic tree, other
traits, and optional covariates. It supports continuous, count, binary,
categorical, ordinal, proportion, zero-inflated-count, and compositional traits.
Traits are encoded on numeric latent scales. The phylogenetic baseline borrows
information across relatives and, on eligible joint routes, across traits. An
optional graph neural network (`gnn = TRUE`) predicts an alternative that is
blended with the baseline using weights calibrated on held-out cells.

## Installation

The workflow below requires pigauto 0.11.1 or later. Install the CRAN
release and check the version:

```r
install.packages("pigauto")
packageVersion("pigauto")
```

If CRAN still provides an older version, or you want the development version,
install from GitHub:

```r
# install.packages("pak") # once, if needed
pak::pak("itchyshin/pigauto")
```

The optional neural network also requires `torch::install_torch()` once.
The default `gnn = FALSE` workflow does not require the torch runtime.

### How pigauto works

The default fit uses a statistical model of trait variation across the tree.
The optional neural network adds another source of predictions. In the
`impute()` workflow:

1. Traits are converted to numbers the model can work with. Measurements are
   standardized, sometimes after a log transformation; ordered categories use
   ranks; unordered categories use one column per category. The package keeps
   the information needed to convert predictions back to the original units
   and labels.
2. Some known values are deliberately hidden. The remaining values are used
   to fit the phylogenetic model, which predicts from relatives and, where
   the joint model is used, from other traits. The hidden values supply
   separate validation and test sets: answers against which predictions can
   be checked.
3. With `gnn = TRUE`, a graph neural network also predicts the traits, using
   connections between species and any supplied covariates. On the validation
   set, pigauto chooses how much weight to give its predictions for each
   trait. The neural weight can be zero. With the default `gnn = FALSE`, no
   network is fitted. If `safety_floor = TRUE`, validation may also choose
   a mean or class-frequency fallback; this option is off by default.
4. The differences between predictions and hidden values measure prediction
   error. For supported traits, validation errors also set the width of
   *conformal prediction bounds*: lower and upper limits intended to contain
   a missing value at a nominal rate, usually 95%. This requires those errors
   to represent the errors at the cells being predicted. Missing values
   concentrated in a poorly represented clade can violate that condition.
   Reusing validation data for model selection also limits the usual
   split-conformal guarantee. Treat the bounds as diagnostics; package-wide
   95% coverage has not been established.
5. Production prediction uses a baseline fitted to all available observations
   and converts predictions back to the original trait scales. `impute()`
   preserves measured values and fills modeled missing cells in one table.
   The separate held-out evaluation uses the baseline and inputs that kept
   the validation/test cells hidden.

`multi_impute()` adds a choice of how to generate several completed tables.
Its `"conformal"` and `"mc_dropout"` routes use the fitted prediction workflow
above. Its `"posterior"` route fits a separate Bayesian phylogenetic model
for continuous traits and uses no neural network. The function map below
explains which output to use for each task.

## Start here

Start with this short workflow on the dataset and tree shipped with pigauto
(300 bird species from AVONET and a matching BirdTree phylogeny):

```r
library(pigauto)

# Load the data and tree shipped with pigauto and assign to common objects
traits <- avonet300
tree   <- tree300

# Check the traits and tree match and are ready to impute
check_pigauto(traits, tree, species_col = "Species_Key")

# Fill in missing trait values m times, using the tree and the other
# traits, to inspect variation among completed datasets
mi <- multi_impute(traits, tree, species_col = "Species_Key", m = 20)

# Extract all m datasets. Before fitting downstream models, check the
# draw method and inference limitations described below.
datasets <- get_imputed_datasets(mi)

# Write an HTML report on the fitted imputation model behind the
# datasets (mi$fit): how well it predicted each trait on held-out data
pigauto_report(mi$fit)
```

Inspect the route used for these datasets:

```r
mi$draws_method
# "conformal" for this example
```

Here the data contain mixed trait types, so the returned datasets are
prediction diagnostics. Compare them to see how the sampled missing values
vary; `with_imputations()` and `pool_mi()` refuse these draws for downstream
inference. Increasing `m` supplies more draws from the same procedure; it
does not make that procedure suitable for inference.

For your own files with one row per species, `read_traits()` moves the
species column into row names. Check and impute those data without a
`species_col` argument:

```r
traits <- read_traits("traits.csv", species_col = "species")
tree <- read_tree("tree.nwk")
check_pigauto(traits, tree)
result <- impute(traits, tree)
```

For several completed tables, use `multi_impute(traits, tree, m = 20)` and
inspect its `draws_method` as described below. A posterior result does not
contain a fitted prediction model in `mi$fit`; use an `impute()` result for
the prediction report in that case.

### Which function should I use?

| Your task | Function | What you receive or do next |
|---|---|---|
| Read trait and tree files | `read_traits()`, `read_tree()` | A trait data frame and a phylogenetic tree. |
| Check species names, trait types and missing values | `check_pigauto(traits, tree)` | An input-check report; no imputation model is fitted. |
| Inspect where values are missing | `preprocess_traits()`, then `plot_missingness()` | Trait summaries, matrices, tree-ordered rows or co-missingness fractions. |
| Obtain one point completion for reporting or exploration | `impute()`, then `completed_data(result)` | One table with observed values preserved and modeled missing cells filled. |
| Inspect variation across several completed tables | `multi_impute(m = 20)`, then `get_imputed_datasets(mi)` | A list of 20 tables. Check `mi$draws_method` and the restrictions below before fitting analyses. |
| Check predictions against hidden known values | `evaluate(fit)` or `cross_validate()`; `pigauto_report(fit)` for an HTML report | Held-out metrics. For lower-level scoring, use `predict(fit, evaluation = TRUE)`. |
| Fit an analysis separately to eligible imputations and combine fixed-effect estimates | `with_imputations()` then `pool_mi()` | A pooled analysis under the documented conditions, not certification that the imputation model suits your analysis. |
| Impute one incomplete continuous predictor while accounting for a specified analysis | `multi_impute_analysis()` | A separate analysis-aware route for its documented `lm`, binomial-logit `glm`, or Gaussian random-intercept `lmer` settings. |

`impute(n_imputations = 20)` still returns **one** completed table. That
legacy argument requests stochastic prediction-diagnostic draws and combines
them into a point prediction; it does not select the posterior sampler.
Use `multi_impute(m = 20)` when you need the separate completed tables and
an explicit record of how their missing values were drawn.

For `multi_impute()`, the default `draws_method = "auto"` selects
`"posterior"` only for all-continuous traits with one row per species,
species names in the row names (`species_col = NULL`), and no covariates or
compositional groups. Otherwise it selects `"conformal"` and explains why.

- **`"conformal"` or `"mc_dropout"`:** use the outputs for prediction
  diagnostics. Do not fit and pool downstream analyses from these draws.
- **`"posterior"`:** inspect `mi$posterior$converged` and
  `mi$posterior$diagnostics`, retain parameter uncertainty, and check the
  analysis conditions in `?multi_impute`. The plug-in setting
  `posterior_control$param_uncertainty = "none"` is diagnostic only and is
  refused by the pooling workflow. Even with convergence and full parameter
  uncertainty, current validation has not established calibrated downstream
  inference; see the results below.

The [multiple-imputation guide](https://itchyshin.github.io/pigauto/articles/multiple-imputation.html)
explains these routes. The [visualisation guide](https://itchyshin.github.io/pigauto/articles/visualisation.html)
explains the report's metrics and missingness plots. The
[function reference](https://itchyshin.github.io/pigauto/reference/index.html)
lists the remaining controls and plotting functions.

BirdTree distributes its tree files in Newick format. Download one from
[BirdTree](https://birdtree.org/downloads/) and read it with `read_tree()`;
`read_tree()` also accepts NEXUS files from other sources. For a downloaded
BirdTree tree, cite Jetz et al. (2012) and BirdTree.org. pigauto uses
the tree you provide and
does not download BirdTree data automatically. The bundled `tree300` is a
teaching example, not a consensus or recommended analysis tree. It was built
specifically for this package and documented with the data help. Jetz et al. describe
the BirdTree source collection; cite BirdTree.org when you use that site.

Users can also run other, tiny installed examples. These are available at
`system.file("extdata", "novice_traits.csv", package = "pigauto")` and
`system.file("extdata", "novice_tree.nwk", package = "pigauto")`. The
matching installed six-expression script is
`system.file("examples", "novice-workflow.R", package = "pigauto")`.

## Example Imputations with Different Traits

### Ordered ecological states

Declare the scientific order of the states explicitly: alphabetical order
(CR, EN, LC, NT, VU) would be wrong for this example. For CSV input, declare
it when reading:

```r
# traits <- read_traits("traits.csv", ordinal_levels = list(
#   threat_status = c("LC", "NT", "VU", "EN", "CR")))
```

For a data frame already in memory, pass the same declaration to `impute()`
(`multi_impute()` and `check_pigauto()` accept it too):

```r
result <- impute(traits, tree, ordinal_levels = list(
  threat_status = c("LC", "NT", "VU", "EN", "CR")))
```

Constructing an ordered factor yourself is equivalent:

```r
traits$threat_status <- ordered(
  traits$threat_status,
  levels = c("LC", "NT", "VU", "EN", "CR")
)
result <- impute(traits, tree)
```

Two mistakes stop with an error rather than guess:

- Declaring the type without the order. `trait_types = c(threat_status =
  "ordinal")` on a plain text column errors and asks for `ordinal_levels`,
  because the type label does not say which state is lowest and which is
  highest.
- Declaring levels that miss a value found in the data. A misspelled level
  errors instead of silently turning that value into missing data.

### Integer-valued continuous measurements

An integer column is otherwise read as a count. Declare a measurement stored
as whole numbers explicitly.

```r
result <- impute(
  traits, tree,
  trait_types = c(body_length_mm = "continuous")
)
```

### Proportions and zero-inflated counts

Declare proportions only when their numeric 0–1 scale is semantically a
proportion. Declare `zi_count` only when zero is a separate process from the
positive count magnitude.

```r
result <- impute(
  traits, tree,
  trait_types = c(habitat_cover = "proportion", parasite_load = "zi_count")
)
```

### Compositions

Each declared composition is one trait: every row must have all components or
all components missing, and observed rows must sum to one.

```r
result <- impute(
  traits, tree,
  multi_proportion_groups = list(diet = c("diet_insect", "diet_fruit", "diet_seed"))
)
```

### Repeated species observations

`read_traits()` requires unique species names. If your file has several
measurements for a species, keep one measurement per row and include a species
ID column that matches the tree tips. Read the file directly, then give the
ID column name to both the input check and `impute()`:

```r
# Read repeated measurements, with one observation per row
traits_obs <- read.csv("traits_repeated.csv", stringsAsFactors = FALSE)

# Check that the data frame is ready for use with pigauto
check_pigauto(traits_obs, tree, species_col = "species")

# Do single imputation
result <- impute(traits_obs, tree, species_col = "species")
```

Read and resolve input errors or relevant warnings before fitting.

Each row is one *observation*: one measurement of that species. The
phylogenetic baseline combines available observations into species summaries,
using rules that depend on the trait type and baseline route:

- numeric traits (continuous, count, ordinal, proportion): the mean on the
  encoded scale, after any log, log1p, rank or logit transformation and
  standardization; this need not equal the raw-data arithmetic mean;
- binary and categorical traits: the joint threshold and one-versus-rest
  routes use the most common class by default. Label-propagation fallback
  routes instead retain within-species class proportions.

In the joint routes, taking the most common class can hide real variation:
a bird recorded as "migratory" in 6 of 10 rows is treated like one recorded
as "migratory" in all 10. For these routes, set
`multi_obs_aggregation = "soft"` to retain the 6-out-of-10 proportion as
partial evidence instead of rounding it to a single class. Label propagation
already uses observed proportions in its species summaries.

With the default baseline fit and a single prediction draw, missing cells
for the same trait and species receive the same species-level prediction.
To let point predictions reflect each observation's covariates, supply
per-row covariates (for example, measurement temperature) and use
`gnn = TRUE`. The network and a separate learned linear covariate term can
both affect predictions. The calibrated weight controls the network's
prediction blend, not every covariate contribution.

With `gnn = FALSE`, covariate use depends on the baseline route. A selected
per-column numeric Brownian-motion fit can use covariates averaged within
species. This can change a species-level prediction when that species has
no remaining observations of the trait. Joint MVN, threshold-joint and
one-versus-rest categorical baselines omit covariates; a numeric fallback
within the same mixed fit may still use them. These baseline routes do not
use per-row covariates to distinguish observations within a species. Read
the route-specific warnings and `?fit_baseline` when supplying covariates.
Stochastic completion draws can differ between rows even without covariates;
that draw variation is not an observation-specific covariate effect.

### When data and tree species names do not match

The trait table and the tree usually come from different sources, so the same
species can be spelled differently in each: `Parus_major` in the tree but
`Parus major`, a synonym, or a typo in the data. pigauto matches names
exactly, letter for letter, and only a species whose name appears in **both**
can be imputed. If some species match, unmatched data rows are retained
with a warning and their missing values stay missing. If no species match,
the input check reports an error and imputation stops.

So before fixing anything by hand, ask the input check which names fall where:

```r
check <- check_pigauto(traits, tree)
check$species$matched    # in both data and tree: modeled, missing values filled
check$species$data_only  # in the data only: row kept in the output, nothing filled
check$species$tree_only  # in the tree only: used while fitting, add no output rows
```

If `data_only` lists species you expected to be imputed, correct their names
in the trait table to match the tree's tip labels and run the check again.

## Defaults, and when to change them

Choose the function for your intended use, then check its assumptions and
validation limits.

> **Planning to analyse the imputed data afterwards? Read this first.**
> Downstream estimates can be biased or overconfident when the imputation
> model does not suit the analysis or fails to carry its uncertainty forward.
> In the registered study of regression slopes after posterior imputation,
> five of eight settings failed the prespecified criteria and three were
> indeterminate. None passed. These results do not validate that route for
> downstream inference. The
> [registered results](https://itchyshin.github.io/pigauto/articles/validation-results.html#inferential-results)
> describe the tested settings and criteria.

| You want | Use | Default? |
|---|---|---|
| Filled-in traits | `impute(traits, tree)` | Yes |
| A lower and upper bound around one imputed value | the conformal bounds in the prediction output, available for continuous, count, ordinal and proportion traits. For zero-inflated counts the bounds describe how large a *nonzero* count would be, and must be read together with the predicted probability that the count is zero | When scores are available for the trait |
| To fit a linear-Gaussian model using imputed continuous traits | `multi_impute(traits, tree)`, then `with_imputations()` and `pool_mi()` | Under the documented input restrictions, for variables included in the imputation model and analysed on their imputation scale; `"auto"` then selects posterior draws |
| To fit `lm`, binomial-logit `glm`, or Gaussian random-intercept `lmer` where one predictor has missing values | `multi_impute_analysis()` | No |

**Why only some draws may be combined.** Rubin's rules combine estimates
and their uncertainty across analyses of the completed datasets. That
requires imputations that represent uncertainty under a model suitable for
the analysis. The posterior route samples missing values and model
parameters under its specified Bayesian phylogenetic model; that alone does
not establish calibrated inference for a downstream analysis.

The `"conformal"` route converts a residual-error quantile into a sampling
scale, while `"mc_dropout"` uses baseline perturbations and, when enabled,
neural dropout. These procedures do not supply the joint model and parameter
uncertainty required by the documented inferential workflow. In the
phylogenetic regression simulations described in `?multi_impute`, conformal
draws biased the pooled slope towards zero. This is evidence from those
settings, not a universal direction or size of bias. `with_imputations()`
refuses conformal and MC-dropout draws.

**Combining Bayesian analyses.** For suitable imputations and compatible
analyses, Bayesian fits such as brms and MCMCglmm can be combined by stacking
posterior draws from the separately fitted datasets, for example with
`brms::brm_multiple()`. This does not make diagnostic imputations valid for
inference. See the [multiple-imputation article](https://itchyshin.github.io/pigauto/articles/multiple-imputation.html)
for conditions and examples using nlme, lme4, glmmTMB, drmTMB, gllvmTMB and
brms.

**Why the neural network is off by default.** pigauto can add a graph
neural network on top of the phylogenetic baseline (`gnn = TRUE`), and the
default is off. In the earlier single-observation benchmarks it did not predict held-out values better
than the baseline and took longer to fit. One later study with repeated
observations and covariates did find gains in two prespecified scenarios,
but none in the no-signal control, and a separately fitted GLS had lower
error in the interactive scenario. Those two studies used different data and
different comparisons, so together they do not show the network is generally
better, and the default stays off. See the
[validation results and their limits](https://itchyshin.github.io/pigauto/articles/validation-results.html#when-does-the-neural-component-add-predictive-value)
for the repeated-observation study.

## Prediction outputs and uncertainty

- **`completed` is the filled-in table.** It retains measured values and
  fills only modeled missing cells. Use it for point summaries and
  exploration. `result$prediction` contains predictions across cells,
  including observed ones. For an honest prediction check, use
  `evaluate(result$fit)` or `predict(result$fit, evaluation = TRUE)`;
  comparing production predictions with their observed inputs is not a
  held-out test.

- **`se` depends on the trait type.** For continuous, count, ordinal and
  proportion traits, it is a model-based uncertainty estimate, with
  back-transformation approximations where needed. For binary and
  categorical traits it is an uncertainty score, not a Gaussian standard
  error. Inspect `result$prediction$probabilities` for these class traits;
  do not calculate their intervals as estimate plus or minus `1.96 * se`.

- **Conformal bounds describe prediction uncertainty under conditions.**
  Where scores are available, continuous, count, ordinal and proportion
  traits have lower and upper prediction bounds. Calibration errors must
  represent errors at the predicted cells, the exchangeability condition
  behind split conformal prediction. Missing an entire clade can violate
  this condition when its errors differ from those in the calibration set.
  Validation data also inform model selection and gate calibration, so the
  standard split-conformal guarantee does not automatically apply. Wider
  intervals are not proof of adequate coverage. For zero-inflated counts,
  bounds describe positive-count magnitude and must be read with the
  probability of a nonzero count. Compositional (`multi_proportion`) traits
  have no conformal bounds. See the [interval limitations](https://itchyshin.github.io/pigauto/articles/gnn-architecture.html#uncertainty-quantification).

- **One completion does not carry imputation uncertainty into an analysis.**
  A model that treats filled values as measured ignores their uncertainty
  and can give biased estimates or intervals that are too narrow. The
  inferential workflows begin with eligible posterior imputations or
  `multi_impute_analysis()` in its documented narrow setting. Generating
  several tables is not enough: the imputation model must also suit the
  analysis, and the validation limits above still apply.

## Advanced controls

The default journey above is the supported starting point. These controls are
for an explicit technical workflow and are kept separate from it.

```r
# Solver, EM, and refinement controls
result <- impute(
  traits, tree,
  joint_solver = "rphylopars",
  em_iterations = 2L,
  joint_refine_iter = 1L,
  conformal_split_val = TRUE
)

# Prediction route: the default "auto" chooses per trait between the
# cross-trait "exact" conditional and "per_column"; force one with
result <- impute(traits, tree, predict_method = "exact")
```

See `?impute`, `?fit_baseline`, and the articles for all calibration,
refinement, and prediction controls.

## More help

- [Getting started](https://itchyshin.github.io/pigauto/articles/getting-started.html)
- [Mixed-type traits](https://itchyshin.github.io/pigauto/articles/mixed-types.html)
- [Multiple imputation and downstream models](https://itchyshin.github.io/pigauto/articles/multiple-imputation.html)
- [Common pitfalls](https://itchyshin.github.io/pigauto/articles/common-pitfalls.html)
- [Tree prediction sensitivity](https://itchyshin.github.io/pigauto/articles/tree-uncertainty.html)

## Citation

Use `citation("pigauto")` in R to get the current package citation.
