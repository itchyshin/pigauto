# pigauto: Phylogenetic Imputation via Graph AUTO-encoders <img src="man/figures/logo.png" align="right" height="139" alt="pigauto logo"/>

**Missing trait data should not stop a comparative analysis.**

Most trait databases are incomplete. Existing imputation tools handle one
trait type at a time, ignore cross-trait correlations, or lack uncertainty
quantification suitable for downstream inference. pigauto addresses all three
gaps by combining three sources of information:
(1) the **phylogenetic tree**, which tells us that closely related
species tend to share similar traits, (2) **correlations among the
traits themselves**, which let observed traits inform predictions of
missing ones, and (3) **environmental covariates** (climate,
geography, habitat), which capture trait variation that phylogeny alone
cannot explain. The package handles continuous measurements, counts,
binary variables, ordered categories, and unordered categories — all
in a single call.

Under the hood pigauto blends a phylogenetic baseline (Brownian-motion
imputation for continuous traits, phylogenetic label propagation for
discrete traits) with a graph neural network that learns additional
patterns from tree topology, inter-trait correlations, and any
user-supplied covariates. A validation-calibrated gate controls how
much the neural network contributes: when the phylogenetic baseline is
already good enough, the gate closes and the network stays out of the
way. The same gated safety applies to covariates — they only
contribute when they demonstrably improve imputation accuracy.

## Key features

- **Multiple trait types**: continuous, binary, categorical, ordinal,
  count, and proportion — all in a single model, no separate analyses
- **Uses the phylogeny**: closely related species inform predictions,
  exactly as you would expect from a comparative method
- **Learns cross-trait patterns**: if body mass predicts beak length,
  observed masses help impute missing beak lengths
- **Environmental covariates**: supply climate, habitat, or geographic
  variables to improve imputation when trait variation has a strong
  environmental component beyond phylogenetic signal
- **Safe by default**: a per-trait gate prevents the neural network
  from degrading traits the phylogenetic baseline already handles well;
  the same gated safety applies to covariates
- **Uncertainty quantification**: conformal prediction intervals give
  95% coverage for continuous, count, and ordinal traits
- **Multiple imputation**: `multi_impute()` → `with_imputations()` →
  `pool_mi()` implements Rubin's rules for downstream inference
- **Tree uncertainty**: `multi_impute_trees()` imputes across a
  posterior sample of trees so that phylogenetic uncertainty propagates
  into standard errors via Rubin's rules
- **Scales to large trees**: tested up to 10,000 species on a laptop
- **Multiple observations per species**: handles intraspecific data

## Documentation

- **Live site**: <https://itchyshin.github.io/pigauto> — tutorials,
  function reference, and changelog as clickable web pages, rebuilt
  on every push to `main`.
- **Documentation index**: [`DOCS.md`](DOCS.md) — the full map of
  every README / tutorial / vignette / man page / benchmark report
  in this repo, with in-repo paths and live-site URLs side by side.
- **First-stop tutorial**: `vignettes/getting-started.Rmd` (rendered
  at [articles/getting-started.html](https://itchyshin.github.io/pigauto/articles/getting-started.html)).
- **Full PCM workflow** on mixed-type traits (continuous + categorical
  + ordinal): live site at
  [pigauto_workflow_mixed.html](https://itchyshin.github.io/pigauto/pigauto_workflow_mixed.html),
  source at [`script/make_workflow_mixed_html.R`](script/make_workflow_mixed_html.R).
- **Architecture notes** for contributors: [`CLAUDE.md`](CLAUDE.md).

## Installation

```r
pak::pak("itchyshin/pigauto")

# First-time torch setup (required)
torch::install_torch()
```

## Quick start

```r
library(pigauto)
data(avonet300, tree300)

# Set species as rownames
df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

# One-call imputation (attention + calibration + conformal by default)
result <- impute(df, tree300)

# Access results
result$prediction$imputed$Mass              # imputed values
result$prediction$se[, "Mass"]             # standard errors
result$prediction$conformal_lower[, "Mass"] # 95% CI lower
result$prediction$conformal_upper[, "Mass"] # 95% CI upper
result$prediction$probabilities$Trophic.Level  # categorical probabilities

# Auto-generate an HTML report
pigauto_report(result)
```

## Using environmental covariates

When trait variation has a strong environmental component, supplying
covariates can improve imputation accuracy. The same covariates
typically serve as predictors in the downstream regression — pigauto
uses them for both purposes:

```r
data(delhey5809, tree_delhey)
df <- delhey5809
rownames(df) <- df$Species_Key

# Traits to impute (some lightness values missing)
traits <- df[, c("lightness_male", "lightness_female")]

# Covariates: fully observed, help imputation AND are regression predictors
covs <- df[, c("annual_mean_temperature", "annual_precipitation",
               "percent_tree_cover", "midLatitude")]

# Impute with all three information sources
result <- impute(traits, tree_delhey, covariates = covs)
```

Covariates must be fully observed (no NAs) and numeric. The gated
safety ensures they never degrade imputation — when phylogenetic
signal already explains the data well, the covariate pathway closes
automatically. See the [comparative study walkthrough](https://itchyshin.github.io/pigauto/pigauto_walkthrough_covariates.html)
for a complete example with downstream regression and Rubin's rules.

## Multiple imputation for downstream analysis

Imputation is a *means*, not an end. The reason to fill in missing cells
is almost always to run a downstream analysis — a phylogenetic
regression, a comparative study of trait evolution, a trait-environment
test. Plugging point-estimate imputations directly into a regression
underestimates standard errors, because it treats imputed cells as if
they were observed. The consequences are well documented for ecological
and phylogenetic datasets in Nakagawa & Freckleton
([2008, *TREE*](https://doi.org/10.1016/j.tree.2008.06.014)) and
([2011, *Behav Ecol Sociobiol*](https://doi.org/10.1007/s00265-010-1044-7)).

The remedy, due to Rubin (1987), is **multiple imputation**: generate
`M` stochastic completions of the trait matrix, fit the downstream
model on each, and pool the `M` outputs via Rubin's rules. pigauto
exposes this workflow as three functions:

```r
library(pigauto)
library(glmmTMB)   # attach so propto() is visible to the formula parser
library(ape)

data(avonet300, tree300)
df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL

# 1. Generate 100 complete datasets (MC-dropout sampling from the GNN)
mi <- multi_impute(df, tree300, m = 100)

# 2. Fit a phylogenetic mixed model on each imputation.
#    Vphy must be a CORRELATION matrix, not a covariance: glmmTMB's
#    propto() estimates sigma^2 freely, so passing the raw vcv(tree)
#    (diagonal = tree height) would rescale the variance component.
Vphy <- cov2cor(ape::vcv(tree300))
fits <- with_imputations(mi, function(d) {
  d$species <- factor(rownames(d), levels = rownames(Vphy))
  d$dummy   <- factor(1)
  glmmTMB(
    log(Mass) ~ log(Wing.Length) + Trophic.Level +
      propto(0 + species | dummy, Vphy),
    data = d
  )
})

# 3. Pool with Rubin's rules
pool_mi(
  fits,
  coef_fun = function(f) fixef(f)$cond,
  vcov_fun = function(f) vcov(f)$cond
)
#>               term estimate std.error   df statistic p.value conf.low conf.high    fmi   riv
#>        (Intercept)   ...
#>   log(Wing.Length)   ...
#>     Trophic.LevelC   ...
```

The pooled table reports `fmi` (fraction of missing information) and
`riv` (relative increase in variance due to non-response) for each
coefficient. A rule of thumb from Rubin is `M ≥ 100 × fmi`: small `fmi`
(< 0.1) means `M = 20` is plenty, but `fmi > 0.3` on any coefficient
means you should push `M` into the hundreds before the pooled standard
errors stop drifting.

`pool_mi()` uses Rubin's rules, which are the correct tool for
frequentist fits (`lm`, `nlme::gls`, `glmmTMB`, `phylolm`, etc.). For
a fully Bayesian workflow you want to concatenate posterior samples
across imputations instead; see the companion **BACE** package's
`pool_posteriors()` on top of `MCMCglmm`. `pool_mi()` deliberately
rejects `MCMCglmm` fits to keep the two paradigms from being mixed.

## Pipeline functions

For fine-grained control:

```r
pd       <- preprocess_traits(df, tree300, log_transform = TRUE)
splits   <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 42)
graph    <- build_phylo_graph(tree300)  # adaptive k_eigen
baseline <- fit_baseline(pd, tree300, splits = splits)
fit      <- fit_pigauto(pd, tree300, splits = splits, graph = graph,
                        baseline = baseline, epochs = 2000L)
pred     <- predict(fit, return_se = TRUE, n_imputations = 10L)

# Evaluate on held-out test set
eval_df <- evaluate(fit, data = pd, splits = splits)

# Generate report
pigauto_report(fit, data = pd, splits = splits)
```

## Evaluation and cross-validation

```r
# Evaluate on test set
eval_df <- evaluate(fit, data = pd, splits = splits)
#>   trait              type        metric    value
#>   Mass               continuous  rmse      0.222
#>   Mass               continuous  pearson_r 0.977
#>   Trophic.Level      categorical accuracy  0.535

# k-fold cross-validation
cv <- cross_validate(pd, tree300, k = 5, seeds = 1:3)
summary(cv)

# Compare methods
comp <- compare_methods(pd, tree300, seeds = 1:3, epochs = 1000)
plot_comparison(comp)

# Simulation benchmark (BM, OU, regime shift, non-linear, mixed types)
bench <- simulate_benchmark(n_species = 100, n_reps = 3, epochs = 500)
summary(bench)
plot(bench)
```

## Plotting

```r
# Training history
plot(fit, type = "history")

# Calibrated gate values
plot(fit, type = "gates")

# Observed vs predicted scatter
plot(pred, data = pd, splits = splits, type = "scatter")

# Conformal interval visualisation
plot(pred, data = pd, splits = splits, type = "intervals", trait = "Mass")
```

## Multiple observations per species

When trait data has multiple rows per species (e.g. individual-level
measurements at different conditions), use `species_col` to identify
the species column.  Observation-level covariates are handled
naturally: the GNN's refinement layer re-injects covariate information
after phylogenetic message passing so that each observation gets its
own prediction:

```r
data(ctmax_sim, tree300)

# ctmax_sim has columns: species, acclim_temp, CTmax
# Multiple observations per species at different acclimation temperatures
# Observation-level covariate: acclim_temp varies within species

traits <- ctmax_sim[, c("species", "CTmax")]
covs   <- data.frame(acclim_temp = ctmax_sim$acclim_temp)

result <- impute(traits, tree300, species_col = "species",
                 covariates = covs)

# Each observation gets a covariate-conditional prediction:
# CTmax at 20°C acclimation ≠ CTmax at 30°C acclimation, even
# for the same species.
```

The GNN aggregates observations to species level before phylogenetic
message passing, then broadcasts back to observation level.  When
covariates are supplied, the observation-level refinement MLP
re-injects covariate values after the broadcast, enabling
covariate-conditional predictions within species.

## Architecture

pigauto is a **gated ensemble of two predictors**:

1. A **phylogenetic baseline** with one branch per trait type:
   Brownian motion (conditional on the phylogenetic correlation matrix)
   for continuous, count, and ordinal traits; Gaussian-kernel
   phylogenetic label propagation for binary and categorical traits.
2. A **graph neural network correction** — an internal torch module
   (`ResidualPhyloDAE`: encoder → 2 attention-based message-passing
   layers with layer norm and ResNet-style skip connections → decoder)
   that produces a full per-cell prediction from spectral node features,
   the corrupted latent matrix, any user-supplied covariates, and an
   adjacency-biased attention over the phylogeny.

Prediction is the blend

    pred = (1 - r_cal) * baseline + r_cal * delta_GNN

where `r_cal` is a **per-trait** gate. The GNN and baseline are not
competitors: the GNN is trained end-to-end by minimising the
type-appropriate loss (MSE for continuous/count/ordinal, BCE for
binary, cross-entropy for categorical) on the blend, with a shrinkage
penalty on `delta - baseline` and an L2 penalty on the gate that keeps
the GNN close to the baseline unless it demonstrably helps. After
training, `r_cal` is grid-searched per trait on a held-out validation
split; for discrete traits a split-validation cross-check prevents the
GNN from degrading baseline accuracy.

> **Note.** The internal module name `ResidualPhyloDAE` retains the
> word "residual" because its GNN layers use ResNet-style residual skip
> connections. The GNN's output `delta` is *not* a statistical residual
> `y - baseline`: it is a full prediction that is blended with the
> baseline at inference time.

Post-training pipeline:
1. Gate values are **calibrated** on the validation set (grid search per trait)
2. **Conformal scores** are computed from validation residuals
3. At test time, conformal intervals provide **guaranteed coverage**

## Trait types

| R class    | pigauto type | Encoding         | Loss          | Baseline                   | Notes                      |
|------------|-------------|------------------|---------------|----------------------------|----------------------------|
| numeric    | continuous  | optional log + z | MSE           | Phylogenetic BM            | Auto-detected              |
| integer    | count       | log1p + z        | MSE           | Phylogenetic BM            | Auto-detected              |
| factor(2)  | binary      | 0/1              | BCE           | Phylo label propagation    | Auto-detected              |
| factor(>2) | categorical | one-hot (K cols) | cross-entropy | Phylo label propagation    | Auto-detected              |
| ordered    | ordinal     | integer + z      | MSE           | Phylogenetic BM            | Auto-detected              |
| numeric(0–1) | proportion | logit + z       | MSE           | Phylogenetic BM            | Requires `trait_types` override |
| integer(ZI) | zi_count   | gate + log1p + z | BCE + MSE    | LP + BM                    | Experimental; requires `trait_types` override |

## Benchmark results

**Simulation** (15 scenarios, 3 replicates, 100–300 species):

| Scenario                    | Continuous RMSE | Count RMSE | Discrete Acc |
|-----------------------------|-----------------|------------|-------------|
| BM, low signal              | +1.3%           | —          | tied         |
| BM, high signal             | +1.4%           | —          | tied         |
| Mixed types, high signal    | +1.8%           | +1.8%      | tied         |
| Variable phylo signal       | **+2.5%**       | +2.7%      | tied         |
| Large tree (300 spp)        | +2.0%           | +1.7%      | tied         |
| OU moderate α=2             | +1.1%           | —          | tied         |
| Non-linear correlations     | +0.8%           | —          | tied         |

**AVONET 300** (7 traits, 5 replicates) with calibration:

| Trait             | Type        | BM     | GNN    | Conformal Coverage |
|-------------------|-------------|--------|--------|--------------------|
| Mass              | continuous  | 0.222  | 0.222  | 94.4%              |
| Beak.Length_Culmen| continuous  | 0.347  | 0.347  | 95.8%              |
| Tarsus.Length     | continuous  | 0.280  | 0.280  | 94.9%              |
| Wing.Length       | continuous  | 0.235  | 0.235  | 90.6%              |
| Trophic.Level     | categorical | 53.5%  | 53.5%  | —                  |
| Primary.Lifestyle | categorical | 58.4%  | 58.4%  | —                  |
| Migration         | ordinal     | 0.861  | 0.861  | 95.3%              |

Gate calibration ensures the GNN never degrades performance: when the
BM baseline is optimal, calibrated gates → 0.

## Bundled data

- `avonet300`: 300 bird species with 4 continuous morphometric traits
  (Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length) plus 3
  ecological traits (Trophic.Level [categorical], Primary.Lifestyle
  [categorical], Migration [ordinal])
- `tree300`: pruned BirdTree phylogeny matching `avonet300`
- `trees300`: 10 posterior trees from the BirdTree Hackett backbone,
  pruned to the same 300 tips — for `multi_impute_trees()` examples
- `avonet_full`: the same 7-trait schema extended to all 9,993 bird
  species for which AVONET3 and the BirdTree Stage2 Hackett MCC tree
  agree. Native missingness is preserved (24 NA cells) so users see
  the real-world missingness pattern. Use this for scale benchmarks;
  use `avonet300` for quick examples and unit tests.
- `tree_full`: the matching 9,993-tip BirdTree phylogeny.
- `delhey5809`: plumage lightness and 4 environmental covariates for
  5,809 bird species (Delhey 2019) — for covariate benchmarks
- `tree_delhey`: matching BirdTree phylogeny for `delhey5809`
- `ctmax_sim`: simulated multi-observation CTmax data (1,464 obs across
  300 species, with acclimation temperature as an observation-level
  covariate) — for multi-obs + covariates examples; uses `tree300`

## Citation

Nakagawa S (2026). pigauto: Phylogenetic Imputation via Graph
Autoencoder. R package version 0.6.1.
