# pigauto

Phylogenetic trait imputation via graph neural network. Combines a
phylogenetic baseline (Brownian motion via Rphylopars) with a
residual graph autoencoder that learns corrections from tree topology
and inter-trait correlations. Supports continuous, binary, categorical,
ordinal, and count traits in a unified latent space.

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

# One-call imputation
result <- impute(df, tree300)

# Access imputed values and standard errors
result$prediction$imputed$Mass
result$prediction$se[, "Mass"]
```

For fine-grained control, use the individual pipeline functions:

```r
pd       <- preprocess_traits(df, tree300, log_transform = TRUE)
splits   <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 42)
graph    <- build_phylo_graph(tree300, k_eigen = 8)
baseline <- fit_baseline(pd, tree300, splits = splits)
fit      <- fit_pigauto(pd, tree300, splits = splits, graph = graph,
                        baseline = baseline, epochs = 2000L)
pred     <- predict(fit, return_se = TRUE, n_imputations = 5L)
```

## Architecture

**ResidualPhyloDAE**: encoder (2 layers) -> multi-hop graph message
passing (2 layers with layer normalisation) -> decoder (2 layers) ->
per-column gated blending with phylogenetic baseline.

Prediction formula: `pred = (1-r) * BM_baseline + r * GNN_prediction`,
where `r` is a learnable per-trait gate bounded in (0, 0.8). Under pure
Brownian motion, `r` shrinks to near zero and the baseline is preserved.
Where the BM model is misspecified, `r` opens to allow GNN corrections.

## Trait types

| R class   | pigauto type | Encoding          | Loss           |
|-----------|-------------|-------------------|----------------|
| numeric   | continuous  | optional log + z  | MSE            |
| integer   | count       | log1p + z         | MSE            |
| factor(2) | binary      | 0/1               | BCE            |
| factor(>2)| categorical | one-hot (K cols)  | cross-entropy  |
| ordered   | ordinal     | integer + z       | MSE            |

## Benchmark results

**Simulation** (12 scenarios, 3 replicates, 150 species):

| Trait type   | GNN vs BM baseline |
|-------------|-------------------|
| Continuous  | +0.9% RMSE improvement |
| Count       | +1.6% RMSE improvement |
| Binary      | 0.0% accuracy change |
| Categorical | 0.0% accuracy change |

**AVONET 2000 bird species** (8 continuous morphometric/geographic traits):

- Geographic traits (latitude, longitude, range): +0.7% to +2.0% RMSE improvement
- Morphometric traits (mass, wing, beak): near-neutral to slight overfitting
- 95% CI coverage: 78-96%

## Vignettes

- `vignette("getting-started")` -- full pipeline walkthrough
- `vignette("mixed-types")` -- mixed trait type example

## Citation

Nakagawa S (2026). pigauto: Phylogenetic Imputation via Graph
Autoencoder. R package version 0.2.0.
