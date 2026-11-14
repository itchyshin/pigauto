---
title: "Mixed-Type Trait Imputation"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Mixed-Type Trait Imputation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Overview

Real comparative datasets contain many kinds of traits: body mass
(continuous), clutch size (count), migratory status (yes/no), diet type
(carnivore/herbivore/omnivore), threat status (LC < VU < EN). pigauto
handles all of them in a single model — you do not need separate
imputation runs for different column types.

| Type | R class | Example | Notes |
|------|---------|---------|-------|
| Continuous | `numeric` | Body mass, wing length | Auto-detected |
| Count | `integer` | Clutch size, litter size | Auto-detected |
| Binary | `factor` (2 levels) | Migratory yes/no | Auto-detected |
| Categorical | `factor` (>2 levels) | Diet type, lifestyle | Auto-detected |
| Ordinal | `ordered` | Threat status (LC < VU < EN) | Auto-detected |
| Proportion | `numeric` (0–1) | Habitat cover, diet fraction | Requires `trait_types = "proportion"` override |
| ZI count | `integer` (zero-inflated) | Parasite load, rare behaviour | Requires `trait_types = "zi_count"` override; experimental — accuracy more variable than other types |
| Multi-proportion | K `numeric` columns summing to 1 | Diet composition, plumage-colour fractions, microbiome relative abundances | Requires `multi_proportion_groups = list(<name> = c("col1", ..., "colK"))`; encoded via centred log-ratio (CLR) + per-component z-score |

The first five rows in this table are auto-detected from R column class;
`proportion`, `zi_count`, and `multi_proportion` must be declared explicitly
(`trait_types` or
`multi_proportion_groups`). All eight share the same latent space — the
phylogenetic baseline and GNN correction both operate in this space,
and type-specific logic appears only at encoding, loss computation,
and decoding.

Before fitting a mixed table, use `check_pigauto(traits, tree)`. Its trait
table records the detected type, declaration source, usable matched rows, and
copy-ready corrections; its species report distinguishes matched, data-only,
and tree-only identities.

## Synthetic example


``` r
library(pigauto)
library(ape)

set.seed(42)
n <- 60
tree <- rtree(n)

traits <- data.frame(
  row.names = tree$tip.label,
  mass      = exp(rnorm(n, 3, 0.5)),
  clutch    = as.integer(rpois(n, 3) + 1L),
  migr      = factor(sample(c("no", "yes"), n, replace = TRUE)),
  diet      = factor(sample(c("herb", "carn", "omni"), n, replace = TRUE)),
  threat    = ordered(sample(c("LC", "VU", "EN"), n, replace = TRUE),
                      levels = c("LC", "VU", "EN"))
)
```

## Preprocessing

`preprocess_traits()` auto-detects column types from R classes:


``` r
pd <- preprocess_traits(traits, tree)
#> Auto-detected trait type "count" (log1p-encoded) for integer column(s): clutch. If any is a continuous measurement stored as whole numbers (common after read.csv), override it: trait_types = c(clutch = "continuous").
print(pd)
#> pigauto_data
#>   Species: 60 
#>   Traits:  5 
#>   Types:   binary=1, categorical=1, continuous=1, count=1, ordinal=1 
#>   Latent columns: 7
```

The `trait_map` records each trait's type, levels, latent column
indices, and normalisation parameters:


``` r
str(pd$trait_map, max.level = 1)
#> List of 5
#>  $ mass  :List of 9
#>  $ clutch:List of 9
#>  $ migr  :List of 8
#>  $ diet  :List of 8
#>  $ threat:List of 8
pd$trait_map$diet
#> $name
#> [1] "diet"
#> 
#> $type
#> [1] "categorical"
#> 
#> $n_latent
#> [1] 3
#> 
#> $latent_cols
#> [1] 4 5 6
#> 
#> $levels
#> [1] "carn" "herb" "omni"
#> 
#> $log_transform
#> [1] FALSE
#> 
#> $mean
#> [1] NA
#> 
#> $sd
#> [1] NA
```

## Creating splits

When `trait_map` is supplied, `make_missing_splits()` operates at the
original-trait level.  For categorical traits, all K one-hot columns
are held out together:


``` r
spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.20,
                           seed = 1, trait_map = pd$trait_map)
cat("Val cells (latent):", length(spl$val_idx), "\n")
#> Val cells (latent): 21
cat("Test cells (latent):", length(spl$test_idx), "\n")
#> Test cells (latent): 57
```

## Baseline fitting

The baseline uses phylogenetic conditional MVN for continuous-family latent
columns, and label-propagation or threshold/liability candidates for
discrete-family columns:


``` r
bl <- fit_baseline(pd, tree, splits = spl)
dim(bl$mu)
#> [1] 60  7
```

## Fitting

`fit_pigauto()` calibrates the phylogenetic baseline for every trait
type on the validation cells and computes the conformal intervals. With
the default `gnn = FALSE` no network is trained. With `gnn = TRUE` it also
trains the GNN correction, using type-specific losses and trait-level
corruption masking when a `trait_map` is present.


``` r
fit <- fit_pigauto(
  pd, tree,
  splits = spl,
  verbose = FALSE,
  seed = 1
)
#> Warning: Small validation set for 5 trait(s): mass (n=2), clutch (n=4), migr
#> (n=3), diet (n=3), threat (n=3). Calibrated gate and conformal scores will be
#> noisy for these trait(s). 95% split-conformal coverage is NOT achievable for 3
#> trait(s) with fewer than 19 validation cells (mass (n=2), clutch (n=4), threat
#> (n=3)): the achievable ceiling is n_val / (n_val + 1), which only reaches 0.95
#> at n_val >= 19. See `?fit_pigauto` under 'Calibration at small n' for smoothing
#> options.
print(fit)
#> pigauto_fit
#>   GNN     : off (baseline only)
#>   Species : 60 
#>   Traits  : 5 -- mass, clutch, migr, diet, threat 
#>   Types   : binary=1, categorical=1, continuous=1, count=1, ordinal=1 
#>   Architecture: hidden_dim = 64 | k_eigen = 4 
#>   Best val loss : 1.169 
#>   Test loss     : 1.1277 
#>   Gate calibration: yes
#>   Conformal scores: 3 traits
```

## Prediction and decoding

`predict()` decodes latent predictions back to original types:


``` r
pred <- predict(fit, return_se = TRUE)
head(pred$imputed)
#>          mass clutch migr diet threat
#> t58  9.518026      4  yes carn     LC
#> t8   9.629026      5  yes carn     LC
#> t36 21.377761      3  yes herb     VU
#> t4  12.202983      8  yes carn     LC
#> t22 20.067241      4   no herb     VU
#> t18 16.213927      4   no omni     EN
```

For binary and categorical traits, class probabilities are available:


``` r
# Binary: probability of "yes"
head(pred$probabilities$migr)
#>       t58        t8       t36        t4       t22       t18 
#> 0.7875313 0.7875313 0.7875313 0.7875313 0.2124687 0.2124687

# Categorical: probability of each diet class
head(pred$probabilities$diet)
#>          carn      herb      omni
#> t58 0.6495271 0.1752365 0.1752365
#> t8  0.6495271 0.1752365 0.1752365
#> t36 0.1752365 0.6495271 0.1752365
#> t4  0.6495271 0.1752365 0.1752365
#> t22 0.1752365 0.6495271 0.1752365
#> t18 0.1752365 0.1752365 0.6495271
```

The SE matrix provides type-appropriate uncertainty (see the
uncertainty-mechanisms table in the GNN architecture vignette for what
each column's SE does and does not guarantee):


``` r
head(pred$se)
#>     mass clutch      migr      diet threat
#> t58    0      0 0.2124687 0.3504729      0
#> t8     0      0 0.2124687 0.3504729      0
#> t36    0      0 0.2124687 0.3504729      0
#> t4     0      0 0.2124687 0.3504729      0
#> t22    0      0 0.2124687 0.3504729      0
#> t18    0      0 0.2124687 0.3504729      0
```

## Evaluation

`evaluate_imputation()` dispatches type-specific metrics:


``` r
ev <- evaluate_imputation(pred, pd$X_scaled, spl)
print(ev)
#>    split  trait        type  n rmse pearson_r coverage_95 mae spearman_rho
#> 1    val   mass  continuous  2    0         1           1  NA           NA
#> 2    val clutch       count  4    0         1          NA   0           NA
#> 3    val   migr      binary  3   NA        NA          NA  NA           NA
#> 4    val   diet categorical  3   NA        NA          NA  NA           NA
#> 5    val threat     ordinal  3    0        NA          NA  NA            1
#> 6   test   mass  continuous 13    0         1           1  NA           NA
#> 7   test clutch       count 10    0         1          NA   0           NA
#> 8   test   migr      binary 10   NA        NA          NA  NA           NA
#> 9   test   diet categorical  6   NA        NA          NA  NA           NA
#> 10  test threat     ordinal  6    0        NA          NA  NA            1
#>    accuracy      brier zero_accuracy aitchison rmse_clr simplex_mae
#> 1        NA         NA            NA        NA       NA          NA
#> 2        NA         NA            NA        NA       NA          NA
#> 3         1 0.04514297            NA        NA       NA          NA
#> 4         1         NA            NA        NA       NA          NA
#> 5        NA         NA            NA        NA       NA          NA
#> 6        NA         NA            NA        NA       NA          NA
#> 7        NA         NA            NA        NA       NA          NA
#> 8         1 0.04514297            NA        NA       NA          NA
#> 9         1         NA            NA        NA       NA          NA
#> 10       NA         NA            NA        NA       NA          NA
```

## Experimental analysis-aware multiple imputation

`multi_impute_analysis()` requires the downstream formula before imputation.
The initial backend supports exactly one incomplete continuous covariate under
MAR. It does not yet impute the mixed discrete and compositional outcomes shown
earlier in this vignette for downstream inference.


``` r
# Precompute nonlinear auxiliaries explicitly when required.
analysis_data$z_sq <- analysis_data$z^2
mi <- multi_impute_analysis(
  data = analysis_data, formula = y ~ x + z,
  missing = "x", model = "glm", m = 50L,
  auxiliary = "z_sq", seed = 1L
)

# Fit a downstream model to each
fits <- with_imputations(mi, function(d) {
  glm(y ~ x + z, data = d, family = binomial)
})

# Pool with Rubin's rules
pool_mi(fits)
```

`pool_mi()` returns a tidy data.frame with `estimate`, `std.error`, `p.value`,
`df` (classical Rubin degrees of freedom by default, or Barnard-Rubin when
`df_fun` supplies complete-data residual degrees of freedom), `fmi` (fraction of missing information),
and `riv` (relative increase in variance) per coefficient.

The `"lm"` branch uses proper Bayesian normal-regression MI, `"glm"` uses
`smcfcs` for a binomial-logit model, and `"lmer"` uses `jomo.smc` for one
Gaussian random intercept. Fixed effects only: multiple incomplete columns,
missing discrete covariates, MNAR, interactions/nonlinear terms involving the
incomplete covariate, random slopes, variance components, correlations,
BLUPs/conditional modes, and latent loadings are unsupported.

For prediction diagnostics, the lower-level
`predict(fit, n_imputations = 5L)` interface returns stochastic datasets.
These are not supported for downstream inference. The analysis-aware
`multi_impute_analysis()` route admits its own narrow class of objects
to `with_imputations()` and `pool_mi()`; the posterior route from
`multi_impute()` supports the separately documented continuous-trait
analyses.
