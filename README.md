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

This page describes the GitHub development version, 0.11.0.9000. CRAN
0.11.0 lacks some functions shown here, including `get_imputed_datasets()`
and the `predict(..., evaluation = TRUE)` argument. Use `evaluate(fit)`
for honest held-out evaluation on CRAN 0.11.0. A matching version number
alone did not distinguish earlier GitHub fixes from that CRAN build.

For the released version:

```r
install.packages("pigauto")
packageVersion("pigauto")
```

Install the GitHub version to use the development workflow below:

```r
# install.packages("pak") # once, if needed
pak::pak("itchyshin/pigauto")
# torch::install_torch() # once, if you opt into gnn = TRUE
```

### How pigauto works?

1. Traits are encoded on numeric scales, with a map for decoding them back to
   their original types. The tree supplies the relationships among species.
2. Known cells are held out for validation and testing. Baseline models use
   the remaining observations to predict these cells. Joint routes can borrow
   information from other traits as well as relatives; fallback routes differ.
3. With `gnn = TRUE`, a graph neural network is trained and its contribution
   is calibrated by trait. A zero neural weight keeps the baseline contribution.
   The default, `gnn = FALSE`, fits no neural network.
4. Held-out prediction errors provide diagnostics and, for supported traits,
   conformal prediction bounds. These bounds have assumptions and are not a
   package-wide coverage guarantee.
5. Predictions are decoded to the original trait scales. `impute()` returns
   one completed table; `multi_impute()` returns several completed tables and
   records how the draws were generated. For downstream inference, check software eligibility and fit the analysis
   separately to each dataset. The posterior route uses no GNN; eligibility
   does not establish inferential calibration.

See [Multiple imputation](https://itchyshin.github.io/pigauto/articles/multiple-imputation.html)
for the difference between diagnostic draws and inferential imputations, and
[Visualisation](https://itchyshin.github.io/pigauto/articles/visualisation.html)
for missingness patterns and held-out checks.

## Start here

Start with this short workflow on the dataset and tree shipped with pigauto
(300 bird species from AVONET and a matching BirdTree phylogeny):

```r
library(pigauto)

# Load the data and tree shipped with pigauto and assign to common objects
traits <- avonet300
tree   <- tree300

# Alternatively, read your own trait table (a CSV with a unique `species`
# column) and a Newick or NEXUS phylogenetic tree:
# traits <- read_traits("traits.csv")
# tree   <- read_tree("tree.nwk")

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

Check `mi$draws_method` before analysing the datasets: for mixed-type data
like `avonet300` the draws describe the spread of plausible imputed values
(method `"conformal"`), and Rubin's-rules pooling via `with_imputations()` and
`pool_mi()` is reserved for the proper-imputation routes documented in
`?multi_impute` and `?multi_impute_analysis`.

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

Declare the scientific order of the states — pigauto never guesses it, because
alphabetical order (CR, EN, LC, NT, VU) would be wrong. For CSV input, declare
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
# Read the data file with repeated species names, where each species is a separate row
traits_obs <- read.csv("traits_repeated.csv", stringsAsFactors = FALSE)

# Check that the data frame is ready for use with pigauto
check_pigauto(traits_obs, tree, species_col = "species")

# Do single imputation
result <- impute(traits_obs, tree, species_col = "species")
```

Read and resolve input errors or relevant warnings before fitting.

With repeated rows, pigauto models species, not individual rows. Each row is
one *observation*: one measurement of that species. Before fitting, the rows
for each species are combined into a single species summary, one rule per
trait type:

- numeric traits (continuous, count, ordinal, proportion): the average of
  that species' observed rows;
- class traits (binary, categorical): the most common class across the rows.

Taking the most common class can hide real variation within a species: a bird
recorded as "migratory" in 6 of 10 rows is treated exactly like one recorded
as "migratory" in all 10. Set `multi_obs_aggregation = "soft"` to keep that
6-out-of-10 as partial evidence instead of rounding it to a single class.

Because the model works at the species level, every missing cell of a species
is filled with the same species-level prediction, whichever row it sits in.
Rows of the same species can only receive different predictions when you both
supply per-row covariates (for example, the temperature at which each
measurement was taken) and turn on the neural network with `gnn = TRUE`.
Under the default `gnn = FALSE`, covariates are ignored with a warning.

### When data and tree species names do not match

The trait table and the tree usually come from different sources, so the same
species can be spelled differently in each — `Parus_major` in the tree but
`Parus major`, a synonym, or a typo in the data. pigauto matches names
exactly, letter for letter, and only a species whose name appears in **both**
gets imputed. A name that matches nothing is not an error: `impute()` warns,
keeps that row in the output, and fills none of its missing values.

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

Start with `impute(traits, tree)` and its defaults. The table below maps what
you want the function to use; the paragraphs after it explain *why* the
caveats exist, so you can judge when they apply to you.

> **Planning to analyse the imputed data afterwards? Read this first.**
> Filling in values and then fitting a model to them is the riskiest way to
> use any imputation package: the analysis inherits every error of the
> imputation model, and nothing in the output warns you. So the bar is not
> "can pigauto produce the datasets and combine the fits" — it can — but
> "are the combined estimates right". We tested that directly in a
> registered study of regression slopes estimated after imputation: of the
> eight settings tested, five failed our criteria and three were
> indeterminate. None passed. Until that changes, treat analyses built on
> pigauto imputations as unvalidated, including on the posterior route. The
> [registered results](https://itchyshin.github.io/pigauto/articles/validation-results.html#inferential-results)
> list the settings and the pass/fail criteria.

| You want | Use | Default? |
|---|---|---|
| Filled-in traits | `impute(traits, tree)` | Yes |
| A lower and upper bound around one imputed value | the conformal bounds in the prediction output, available for continuous, count, ordinal and proportion traits. For zero-inflated counts the bounds describe how large a *nonzero* count would be, and must be read together with the predicted probability that the count is zero | Yes, for those trait types |
| To fit a linear model to imputed continuous traits | `multi_impute(traits, tree)`, then `with_imputations()` and `pool_mi()` | Yes, when every trait is continuous with one row per species and no covariates — the default `draws_method = "auto"` then uses posterior draws |
| To fit `lm`, binomial-logit `glm`, or Gaussian random-intercept `lmer` where one predictor has missing values | `multi_impute_analysis()` | No |

**Why only some draws may be combined.** Averaging results across `m`
imputed datasets (Rubin's rules, what `pool_mi()` does) is only valid when
the datasets differ by the imputation model's real uncertainty — so-called
*proper* imputations, which pigauto can currently guarantee only in the
all-continuous setting above. For any other data, `multi_impute()` tells you
it is falling back to `"conformal"` draws. Those draws (and `"mc_dropout"`
draws) are centred on the model's best guess and only *describe* the spread
of plausible values; they vary too little between datasets. Combine them
anyway and your analysis looks more certain than it is and its estimates
shrink: in our simulation of a phylogenetic regression, the pooled slope was
biased towards zero. That is why `with_imputations()` refuses these draws
outright instead of letting an analysis silently inherit the bias. Check
`mi$draws_method` to see which kind you have.

**Why Bayesian models are combined differently.** Bayesian fits (brms,
MCMCglmm) already produce posterior draws, so imputation uncertainty is
carried by fitting the model to each dataset and stacking those draws — for
example with `brms::brm_multiple()` — not by `pool_mi()`. Worked examples
for nlme, lme4, glmmTMB, drmTMB, gllvmTMB and brms are in the
[multiple imputation article](https://itchyshin.github.io/pigauto/articles/multiple-imputation.html).
And throughout: that a model class *runs* and returns coefficients does not
validate the imputations for that analysis — check the article's conditions
before interpreting pooled estimates.

**Why the neural network is off by default.** pigauto can add a graph
neural network on top of the phylogenetic baseline (`gnn = TRUE`), and the
default is off because the network has to earn its place with evidence. In
our single-observation benchmarks it did not predict held-out values better
than the baseline and took longer to fit. One later study with repeated
observations and covariates did find gains in two prespecified scenarios —
but none in the no-signal control, and a separately fitted GLS had lower
error in the interactive scenario. Those two studies used different data and
different comparisons, so together they do not show the network is generally
better, and the default stays off. See the
[validation results and their limits](https://itchyshin.github.io/pigauto/articles/validation-results.html#when-does-the-neural-component-add-predictive-value)
for the repeated-observation study.

## Prediction outputs and uncertainty

The result object holds several outputs that look interchangeable but are
not. Each point below guards against a specific mistake that nothing in the
output would flag — the numbers would simply look more trustworthy than they
are.

- **Analyse `completed`, not `prediction`.** `completed` is your trait table
  with every observed value kept and only the modeled missing cells filled —
  this is the dataset to use. `result$prediction` holds the model's
  prediction for *every* cell, including the ones you actually measured. It
  exists so you can judge the model (how close does it come on known
  values?); analyse it and you have quietly replaced your measurements with
  model output.

- **`se` is only a standard error for numeric traits.** For continuous,
  count, ordinal and proportion traits, `result$prediction$se` is the
  model's standard error for that cell. For binary and categorical traits
  the same column holds an uncertainty *score* — how far the classification
  is from certain — so "estimate ± 1.96 × se" arithmetic is meaningless
  there. Read the class probabilities in `result$prediction$probabilities`
  instead.

- **Conformal bounds are a calibrated diagnostic, not a guarantee.** For
  continuous, count, ordinal and proportion traits, each imputed value comes
  with a lower and upper bound, tuned on held-out cells to aim for 95%
  coverage. That tuning assumes the held-out cells resemble the cells being
  predicted. Missing an entire clade breaks the assumption: the bounds were
  tuned on other clades' errors, which may not represent the missing clade.
  The same held-out cells also help choose the model, which nudges the
  bounds optimistic. So read them as a useful diagnostic, not a certified
  95% interval. For zero-inflated counts the bounds describe how large a
  *nonzero* count would be and must be read with the probability that the
  count is nonzero;
  compositional (`multi_proportion`) traits have no bounds yet. Details:
  [interval limitations](https://itchyshin.github.io/pigauto/articles/gnn-architecture.html#uncertainty-quantification).

- **Do not do inference on one completed dataset!** However rich the
  uncertainty columns look, a single filled-in dataset cannot carry the
  uncertainty of the imputation itself, so a model fitted to it overstates
  certainty. Analyses that need that uncertainty go through `multi_impute()`
  or, in its documented narrow setting, `multi_impute_analysis()`.
  `multi_impute()` defaults to `draws_method = "auto"`, which uses posterior
  draws when the data allow; with any other data it falls back to conformal
  draws and prints the reason. Those fallback draws cannot be pooled — see
  "Why only some draws may be combined" above. And as the warning box above
  explains, being *able* to generate and pool draws still does not establish
  that the resulting inference is calibrated.

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
