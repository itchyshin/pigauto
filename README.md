# pigauto: Phylogenetic Imputation via Graph AUTO-encoders <img src="https://raw.githubusercontent.com/itchyshin/pigauto/main/man/figures/logo.png" align="right" height="139" alt="pigauto logo"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/itchyshin/pigauto/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/itchyshin/pigauto/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/itchyshin/pigauto/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/itchyshin/pigauto/actions/workflows/pkgdown.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

> [!WARNING] **pigauto is experimental; use at your own risk.** It needs
> further validation. Point estimates are the supported claim. Prediction
> intervals are nominal held-out diagnostics, not package-wide certification;
> covariance routes have focused-test evidence only.

pigauto fills missing species trait values using a phylogenetic tree, other
traits, and optional covariates. It supports continuous, count, binary,
categorical, ordinal, proportion, zero-inflated-count, and compositional traits.

## Installation

```r
install.packages("pigauto")
# torch::install_torch() # once, if you opt into gnn = TRUE
```

## Start here

Put a CSV trait table (with a unique `species` column) and a Newick or NEXUS
tree beside your R script. Then use this six-line journey:

```r
library(pigauto)

traits <- read_traits("traits.csv")
tree <- read_tree("tree.nwk")
check_pigauto(traits, tree)
result <- impute(traits, tree)
completed <- completed_data(result)
pigauto_report(result)
```

`impute()` now uses the phylogenetic baseline by default (`gnn = FALSE`),
estimates Pagel's lambda for eligible traits, and leaves the mean safety
floor and phylogenetic-signal gate off. To train the optional graph neural
network, install the torch runtime once and call
`impute(traits, tree, gnn = TRUE)`. The former GNN-on behaviour is an
explicit choice in 0.11.0.

`check_pigauto()` and `completed_data()` were added in 0.11.0. If your
installed version is 0.10.0, call `impute(traits, tree)` directly; that
version does not include `check_pigauto()`.

`check_pigauto()` runs before fitting: it reports input errors, species/tree
matching, trait declarations, runtime availability, and size. Resolve an
`"error"` status before calling `impute()`; a fully observed target table
does not need fitting, so use `cross_validate()` instead.

The output roles are deliberately separate:

- `completed` is the input-shaped trait data with observed values retained and
  modeled missing cells filled.
- `result$prediction` contains all-cell diagnostic predictions, not a second
  completed dataset.
- `result$prediction$se` is type-dependent uncertainty; discrete values are
  not Gaussian standard errors. Conformal bounds are nominal held-out
  diagnostics.
- Split conformal intervals assume the calibration cells resemble the
  scored cells. Masking a whole clade breaks that assumption because
  calibration then comes from other clades.
- A `pigauto_result` does not itself authorize downstream inference. The
  supported pooling entries begin with `multi_impute_analysis()` in its
  documented narrow regime or, for continuous traits,
  `multi_impute(draws_method = "posterior")`. `multi_impute()`'s default,
  `draws_method = "auto"`, picks `"posterior"` whenever the data allow it
  and otherwise says why it fell back to non-poolable conformal draws.

Tiny installed example inputs are available at
`system.file("extdata", "novice_traits.csv", package = "pigauto")` and
`system.file("extdata", "novice_tree.nwk", package = "pigauto")`. The
matching installed six-expression script is
`system.file("examples", "novice-workflow.R", package = "pigauto")`.

## Copy-ready data recipes

### Ordered ecological states

Use an ordered factor when states have a known order.

```r
traits$threat_status <- ordered(
  traits$threat_status,
  levels = c("LC", "NT", "VU", "EN", "CR")
)
result <- impute(traits, tree)
```

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

### Repeated observations

`read_traits()` requires unique species names. For repeated observations,
read the CSV directly and supply `species_col`.

```r
traits_obs <- read.csv("traits_repeated.csv", stringsAsFactors = FALSE)
check_pigauto(traits_obs, tree, species_col = "species")
result <- impute(traits_obs, tree, species_col = "species")
```

### Reconciling species names

Inspect the structured species report before changing data:

```r
check <- check_pigauto(traits, tree)
check$species$matched
check$species$data_only
check$species$tree_only
```

`data_only` rows remain in `completed` but are not modeled or filled.
`tree_only` tips are internal all-missing rows used only while fitting.

## Defaults, and when to change them

`impute(traits, tree)` with its defaults is the starting point. The table
says which route to take for each need; reach for another route only for
the need in the left-hand column.

| You want | Use | Default? |
|---|---|---|
| Filled-in traits | `impute(traits, tree)` | Yes |
| Uncertainty for one imputed value | the per-cell conformal interval in `result$prediction` | Yes |
| A linear downstream analysis of imputed continuous traits on their imputation scale | `multi_impute(traits, tree)`, then `with_imputations()` and `pool_mi()` | Yes, when every trait is continuous, with one row per species and no covariates: the default `draws_method = "auto"` then uses posterior draws |
| One incomplete covariate in `lm`, binomial-logit `glm`, or Gaussian random-intercept `lmer` | `multi_impute_analysis()` | No |

For other data, `draws_method = "auto"` falls back to `"conformal"` draws
and says so. Those draws, like `"mc_dropout"`, describe the spread of
plausible values but are not proper multiple imputations: in a simulation
of a phylogenetic regression they biased the pooled slope towards zero, and
`with_imputations()` refuses them. Bayesian downstream fits (brms, MCMCglmm) are
combined by stacking posterior draws, for example with
`brms::brm_multiple()`, not with `pool_mi()`. Worked examples for nlme,
lme4, glmmTMB, drmTMB, gllvmTMB and brms are in the
[multiple imputation article](https://itchyshin.github.io/pigauto/articles/multiple-imputation.html).
The ability to fit and extract coefficients from a model class does not by
itself validate the imputation model for that analysis; check the article's
conditions before interpreting pooled estimates.

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

Nakagawa S (2026). *pigauto: Phylogenetic Imputation via Graph Autoencoder*.
