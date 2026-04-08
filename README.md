# pigauto

Phylogenetic trait imputation via graph neural network. Combines a
phylogenetic baseline (Brownian motion via Rphylopars) with a residual
graph autoencoder that learns corrections from tree topology and
inter-trait correlations. Supports continuous, binary, categorical,
ordinal, and count traits in a unified latent space.

## Key features

- **Mixed-type traits**: continuous, binary, categorical, ordinal, and count
  in a single model
- **Attention-based message passing**: learned phylogenetic attention with
  adjacency-prior bias — the model attends to informative neighbours
  rather than treating all equally
- **Phylogenetic label propagation**: species-specific baselines for
  discrete traits using phylogenetic similarity, not flat frequencies
- **Validation-calibrated gates**: post-training gate optimisation prevents
  the GNN from adding noise when the baseline is already strong
- **Conformal prediction intervals**: distribution-free 95% coverage
  guarantees for continuous/count/ordinal traits
- **Adaptive spectral encoding**: `k_eigen` scales automatically with tree
  size (4 for 30 tips → 15 for 300 → 32 for 640+)
- **Multiple observations per species**: aggregates individual-level data
  before phylogenetic message passing
- **Auto-generated HTML reports**: `pigauto_report()` produces interactive
  benchmark reports with Chart.js

## Installation

```r
# Install from local source
devtools::install()

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
measurements), use `species_col` to identify the species column:

```r
# traits_df has columns: species, mass, wing_length, ...
# Multiple rows per species are allowed
result <- impute(traits_df, tree, species_col = "species")

# Or with the lower-level API
pd <- preprocess_traits(traits_df, tree, species_col = "species")
pd$n_obs      # number of observations
pd$n_species  # number of unique species
```

The GNN aggregates observations to species level before phylogenetic
message passing, then broadcasts back to observation level.
Rphylopars handles within-species replication natively in the
baseline.

## Architecture

**ResidualPhyloDAE**: encoder (2 layers) → attention-based graph message
passing (2 layers with layer normalisation, phylogenetic attention bias)
→ decoder (2 layers) → per-column gated blending with phylogenetic
baseline.

Prediction: `pred = (1 - r_cal) * BM_baseline + r_cal * GNN_delta`,
where `r_cal` is a per-trait gate optimised on the validation set.

Post-training pipeline:
1. Gate values are **calibrated** on the validation set (grid search per trait)
2. **Conformal scores** are computed from validation residuals
3. At test time, conformal intervals provide **guaranteed coverage**

## Trait types

| R class    | pigauto type | Encoding         | Loss          | Baseline                   |
|------------|-------------|------------------|---------------|----------------------------|
| numeric    | continuous  | optional log + z | MSE           | Rphylopars BM              |
| integer    | count       | log1p + z        | MSE           | Rphylopars BM              |
| factor(2)  | binary      | 0/1              | BCE           | Phylo label propagation    |
| factor(>2) | categorical | one-hot (K cols) | cross-entropy | Phylo label propagation    |
| ordered    | ordinal     | integer + z      | MSE           | Rphylopars BM              |

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

## Citation

Nakagawa S (2026). pigauto: Phylogenetic Imputation via Graph
Autoencoder. R package version 0.3.0.
