# pigauto: Phylogenetic Imputation via Graph AUTO-encoders <img src="https://raw.githubusercontent.com/itchyshin/pigauto/main/man/figures/logo.png" align="right" height="139" alt="pigauto logo"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/itchyshin/pigauto/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/itchyshin/pigauto/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/itchyshin/pigauto/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/itchyshin/pigauto/actions/workflows/pkgdown.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**Missing trait data should not stop a comparative analysis.**

pigauto helps fill missing values in species trait data. Give it a trait table
and a phylogenetic tree, and it uses related species, other traits, and optional
environmental data to estimate the gaps. It works with common continuous and
categorical trait types, and provides prediction intervals where available.

## The workflow

```
Trait data + phylogenetic tree
              ↓
          impute()
              ↓
Completed trait data + prediction intervals
```

## Installation

```r
# CRAN release
install.packages("pigauto")

# Development version
pak::pak("itchyshin/pigauto")

# First-time torch setup (required once)
torch::install_torch()
```

## Quick start

```r
library(pigauto)
library(ape)

data(avonet300, tree300)
df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

# ── Step 1: point imputation ──────────────────────────────────────────────
# avonet300 is fully observed, so there is nothing to impute as-is.
# Demonstrate the workflow by hiding 30 Mass values:
set.seed(1L)
hide <- sample(which(!is.na(df$Mass)), 30L)
df_obs <- df
df_obs$Mass[hide] <- NA

result <- impute(df_obs, tree300)

# `completed` keeps observed values and fills the NAs:
result$completed$Mass[hide]                          # pigauto's imputations
df$Mass[hide]                                        # held-out truth, for comparison

# `imputed_mask` flags which cells were filled:
sum(result$imputed_mask[, "Mass"])                   # 30

# Conformal 95% intervals are stored on the prediction object:
result$prediction$conformal_lower[hide, "Mass"]
result$prediction$conformal_upper[hide, "Mass"]

# `result$prediction$imputed$Mass` contains predictions for ALL species
# (observed + missing) -- intended for diagnostics, not as the imputed
# output.  See `?impute` ("What gets imputed (read this first)") for the
# distinction, and for the `n_imputations = 20, pool_method = "mode"`
# recommendation on imbalanced ordinal traits like Migration.

```

## Using environmental covariates

When trait variation has a strong environmental component, supplying
covariates can improve imputation. The same covariates may also serve as
predictors in the downstream model:

```r
data(ctmax_sim, tree300)
traits <- ctmax_sim[, c("species", "CTmax")]
covs   <- ctmax_sim["acclim_temp"]

result <- impute(traits, tree300, species_col = "species",
                 covariates = covs)
```

Covariates must be fully observed. Numeric columns are z-scored;
factor/ordered columns are one-hot encoded automatically. The calibrated gate
can close when the phylogenetic baseline already explains the validation cells,
so covariates should be treated as an evidence-backed addition rather than an
automatic improvement.

## Phylogenetic tree uncertainty

`multi_impute_trees()` remains an experimental prediction-sensitivity tool.
Its stochastic completions are not intended for downstream inference.

## Trait types

### Automatic type detection

pigauto infers each trait's type from its R class — no `trait_types`
argument needed for most data:

| R class | pigauto type | How to set in R |
|---|---|---|
| `numeric` | continuous | default for `read.csv()` numeric columns |
| `integer` | count | `as.integer(x)` |
| `factor` (2 levels) | binary | `factor(x)` |
| `factor` (>2 levels) | categorical | `factor(x)` |
| `ordered` factor | ordinal | `ordered(x, levels = c("low","mid","high"))` |
| `character` | → factor → binary/categorical | auto-converted |
| `logical` | binary | `as.logical(x)` |

**Two types require explicit override** (indistinguishable from R class):

```r
result <- impute(df, tree,
                 trait_types = c(Survival  = "proportion",
                                 Parasites = "zi_count"))
```

**Compositional (multi-proportion) data** — K columns per row summing to 1
(e.g. plumage-colour proportions, diet composition, microbiome relative
abundances). Declare the group separately because these columns belong
to ONE trait, not K independent ones:

```r
# df has columns black, blue, red, ..., yellow (12 colours that sum to 1)
result <- impute(df, tree,
                 multi_proportion_groups = list(
                   colour = c("black", "blue", "red", "rufous",
                              "white", "yellow", ...)))

# imputed compositions sum to 1 per row:
result$prediction$probabilities$colour   # n_species x K matrix
```

Encoding is centred log-ratio (CLR) + per-component z-score; baseline is
Brownian motion on CLR space; decode is softmax back to the simplex.

### Full type reference

| R class | pigauto type | Encoding | Loss | Baseline |
|---|---|---|---|---|
| numeric | continuous | log (optional) + z | MSE | Phylogenetic BM |
| integer | count | log1p + z | MSE | Phylogenetic BM |
| factor(2) | binary | 0/1 | BCE | Phylo label propagation |
| factor(>2) | categorical | one-hot | cross-entropy | Phylo label propagation |
| ordered | ordinal | integer + z | MSE | Phylogenetic BM |
| numeric(0–1) | proportion | logit + z | MSE | Phylogenetic BM |
| integer(ZI) | zi_count | gate + log1p + z | BCE + MSE | LP + BM |
| K numeric cols summing to 1 | multi_proportion | CLR + per-component z | MSE (CLR) | Per-component BM |

> **Trait vs covariate.** A trait is something you want to impute (NAs
> allowed). A covariate is something that helps imputation (must be fully
> observed). The same variable can be either depending on your question:
> IUCN status unknown for Data Deficient species → put it in `traits` as
> an `ordered` factor; IUCN status fully known → pass as a covariate.

## Architecture

pigauto blends two predictors per trait:

    pred = (1 - r_cal) × baseline + r_cal × GNN_delta

The **baseline** is Brownian-motion conditional imputation for
continuous/count/ordinal/proportion traits, and phylogenetic label
propagation for binary/categorical traits. The **GNN** is an
attention-based graph neural network trained on the phylogenetic
topology, cross-trait correlations, and any user covariates.
`r_cal` is a per-trait gate calibrated on a held-out validation split:
when the baseline is already optimal on the validation cells, the gate can
close to zero so the GNN contributes nothing. This is a validation-calibrated
safety mechanism, not a universal performance guarantee.

## Benchmarks

Current benchmark pages live under the pkgdown Methodology menu and should be
read as versioned evidence: each page states its data-generating regime, sample
size, missingness, seed, and comparison set. We keep headline numbers out of
this README because both pigauto and the reference implementations are still
under active development, and the numbers move with each release.

### Caveats from multi-seed evidence

Three findings worth knowing before running pigauto on AVONET-scale
data:

1. **Continuous-trait results have non-trivial cross-seed
   variance** — multi-seed evidence on AVONET n=1500 (3 seeds at
   N_IMP=20, see `useful/MEMO_2026-05-01_multiseed_n20_and_default_flip.md`)
   shows lift bands of ±10 % to ±20 % around the mean for traits
   like Tarsus.Length and Beak.Length_Culmen. We recommend
   `multi_impute(m = 20, ...)` and reporting mean ± SD across
   ≥ 3 seeds for any quantitative claim. Single-seed RMSE numbers
   should not be interpreted as point estimates of expected
   performance.

2. **Mass on AVONET is currently unstable.** At N_IMP=20 the BM
   joint baseline emits singular-Σ warnings on Mass, and one of
   three seeds produces a tail outlier with RMSE ≈ 11× the
   column-mean baseline (others beat it). Cause is under
   investigation; expected to be a back-transform tail issue. If
   you are imputing Mass, run multiple seeds and inspect the
   distribution of imputed values for outliers.

3. **Migration (3-level ordinal) has needed special care** in the
   AVONET multi-seed checks. The current baseline includes
   per-trait ordinal path selection with an LP-via-OVR candidate for
   low-K ordinal traits, but re-run the multi-seed evidence before
   treating this regime as resolved for a quantitative claim.

### The gate protects against regression

For high-signal discrete traits the calibrated gate closes completely
(correct behaviour — phylogenetic label propagation is already optimal).
The GNN contributes meaningfully for OU-misspecified continuous traits,
zero-inflated counts at high zero-fraction, and low-signal traits when
strong covariates are supplied.

## Multiple observations per species

```r
data(ctmax_sim, tree300)
traits <- ctmax_sim[, c("species", "CTmax")]
covs   <- data.frame(acclim_temp = ctmax_sim$acclim_temp)

result <- impute(traits, tree300, species_col = "species",
                 covariates = covs)
```

Each observation gets a covariate-conditional prediction: the GNN
aggregates to species level for phylogenetic message passing, then
re-injects covariate values at observation level so predictions within
a species differ by covariate context.

## Documentation

- **Live site**: <https://itchyshin.github.io/pigauto/>
- **Getting started**: `vignettes/getting-started.Rmd`
  ([rendered](https://itchyshin.github.io/pigauto/articles/getting-started.html))
- **Mixed types**:
  [mixed-types.html](https://itchyshin.github.io/pigauto/articles/mixed-types.html)
- **Posterior-tree prediction sensitivity**:
  [tree-uncertainty.html](https://itchyshin.github.io/pigauto/articles/tree-uncertainty.html)
- **Common pitfalls**:
  [common-pitfalls.html](https://itchyshin.github.io/pigauto/articles/common-pitfalls.html)

## Acknowledgements

Russell Dinnage ([@rdinnager](https://github.com/rdinnager)) contributed to
graph-scaling work and foundation-model/TabPFN exploration. Bhavya Jain
([@b1805](https://github.com/b1805)) contributed correctness fixes,
safeguards and tests, and a large attention ablation. They are acknowledged
as package contributors; Shinichi Nakagawa is the package author and
maintainer for this release.

## Citation

Nakagawa S (2026). *pigauto: Phylogenetic Imputation via Graph
Autoencoder*. R package version 0.10.0.
