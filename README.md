# pigauto: Phylogenetic Imputation via Graph AUTO-encoders <img src="man/figures/logo.png" align="right" height="139" alt="pigauto logo"/>

**Missing trait data should not stop a comparative analysis.**

pigauto fills gaps in species trait matrices by combining the phylogenetic
tree, cross-trait correlations, and optional environmental covariates — then
propagates imputation uncertainty through to your downstream model via
multiple imputation and Rubin's rules.

## The workflow

```
Raw data (with NAs)
       ↓
   impute()          # point estimates + conformal uncertainty intervals
       ↓
multi_impute(m=50)   # 50 stochastic complete datasets
       ↓
with_imputations()   # fit your model on each dataset
       ↓
   pool_mi()         # pool with Rubin's rules → correct SEs
```

## Installation

```r
pak::pak("itchyshin/pigauto")

# First-time torch setup (required once)
torch::install_torch()
```

## Quick start

```r
library(pigauto)
library(glmmTMB)
library(ape)

data(avonet300, tree300)
df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

# ── Step 1: point imputation ──────────────────────────────────────────────
result <- impute(df, tree300)

result$prediction$imputed$Mass               # imputed values
result$prediction$conformal_lower[, "Mass"]  # 95% CI lower
result$prediction$conformal_upper[, "Mass"]  # 95% CI upper

# ── Step 2: generate 50 complete datasets ────────────────────────────────
mi <- multi_impute(df, tree300, m = 50L)

# ── Step 3: fit downstream model on each ────────────────────────────────
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

# ── Step 4: pool with Rubin's rules ──────────────────────────────────────
pool_mi(
  fits,
  coef_fun = function(f) fixef(f)$cond,
  vcov_fun = function(f) vcov(f)$cond
)
#>               term  estimate std.error    df  fmi
#>        (Intercept)     ...
#>   log(Wing.Length)     ...
#>     Trophic.LevelC     ...
```

The pooled table reports `fmi` (fraction of missing information) per
coefficient. When `fmi > 0.1` on any term, increase `m` until the
standard errors stop changing. For most comparative datasets `m = 50`
is sufficient.

## Using environmental covariates

When trait variation has a strong environmental component, supplying
covariates improves imputation. The same covariates typically serve as
predictors in the downstream model:

```r
data(delhey5809, tree_delhey)
df <- delhey5809
rownames(df) <- df$Species_Key

traits <- df[, c("lightness_male", "lightness_female")]
covs   <- df[, c("annual_mean_temperature", "annual_precipitation",
                 "percent_tree_cover", "midLatitude")]

result <- impute(traits, tree_delhey, covariates = covs)
```

Covariates must be fully observed. Numeric columns are z-scored;
factor/ordered columns are one-hot encoded automatically (v0.6.1+).
A per-trait gated safety prevents covariates from degrading imputation:
when phylogenetic signal explains the data well, the covariate pathway
closes automatically. See the
[covariate walkthrough](https://itchyshin.github.io/pigauto/pigauto_walkthrough_covariates.html)
for a complete example.

## Phylogenetic tree uncertainty

If you have a posterior sample of trees (e.g. 50 trees from BirdTree),
tree uncertainty enters the analysis at **two** places — and pigauto
handles them separately because they are conceptually distinct.

**Step 1 — imputation** (pigauto's job).
`multi_impute_trees()` runs a full pigauto fit on each posterior tree,
so each completed dataset is conditional on a *different* tree:

```r
data(avonet300, trees300)   # trees300 = 50 posterior trees
df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL

mi <- multi_impute_trees(df, trees = trees300, m_per_tree = 5L)
#> 50 trees × 5 imputations = 250 completed datasets
#> mi$tree_index[i]  →  which posterior tree produced dataset i
```

**Step 2 — analysis** (Nakagawa & de Villemereuil 2019).
For each completed dataset, fit the downstream comparative model
using the *same* tree that produced that dataset, then pool the
T × M fits with Rubin's rules:

```r
fits <- Map(
  function(dat, t_idx) {
    dat$species <- rownames(dat)
    nlme::gls(
      log(Mass) ~ log(Wing.Length),
      correlation = ape::corBrownian(phy = trees300[[t_idx]], form = ~species),
      data = dat, method = "ML"
    )
  },
  mi$datasets,
  mi$tree_index
)
pool_mi(fits)
```

The pooled standard errors now reflect imputation uncertainty,
phylogenetic tree uncertainty, and their interaction — all in one
Rubin's-rules step. Reference: Nakagawa S, de Villemereuil P (2019).
*Systematic Biology* 68(4):632–641. doi:
[10.1093/sysbio/syy089](https://doi.org/10.1093/sysbio/syy089).

### Compute cost scales linearly with T

Each posterior tree requires a fresh pigauto fit (no caching possible —
the tree *is* the model). Rough budget on a modern CPU laptop:

| Species n | 1 fit | T = 50 | T = 10 |
|---:|---:|---:|---:|
| 300 | ~30–60 s | 25–50 min | 5–10 min |
| 5,000 | ~5–10 min | 4–8 hr | ~1 hr |
| 10,000 | ~20–40 min | 17–33 hr | 3–7 hr |

**Guidance for large trees.** The 2019 paper notes that 10–20 posterior
trees are usually enough (use the "relative efficiency" index in
`pool_mi()$fmi` to check convergence). For 10,000-species datasets we
recommend T = 10–20, or parallelising the T fits across machines
(each fit is independent — good HPC / cloud use case).
If you have no posterior tree sample, use a single MCC tree via
`impute()` or `multi_impute()` and note the tree-uncertainty caveat in
your paper.

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
when the baseline is already optimal the gate closes to zero and the
GNN contributes nothing, guaranteeing the network never degrades accuracy.

## Benchmarks

Current per-trait-type numbers, scaling curves, AVONET 9,993 results,
tree-uncertainty comparisons, and the covariate-simulation sweep all live
in the
[validation suite](https://itchyshin.github.io/pigauto/dev/validation_suite.html).
We keep the headline numbers out of this README because both pigauto and
the reference implementations (including BACE) are still under active
development, and the numbers move with each release. The validation
suite is regenerated from the actual benchmark scripts every release, so
it is the single source of truth.

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

3. **Migration (3-level ordinal) regresses by ≈ 11 pp** vs the
   per-class-mode baseline on AVONET, consistent across all
   tested seeds. The ordinal threshold-joint baseline (Phase B3)
   is under-determined for K=3 and the calibrated gate cannot
   route around it because LP is not yet a corner in the
   safety-floor simplex grid. Phase F (LP corner for ordinal)
   addresses this and is queued for v0.9.2.

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

- **Live site**: <https://itchyshin.github.io/pigauto>
- **Getting started**: `vignettes/getting-started.Rmd`
  ([rendered](https://itchyshin.github.io/pigauto/articles/getting-started.html))
- **Full PCM workflow** (mixed types + downstream inference):
  [pigauto_workflow_mixed.html](https://itchyshin.github.io/pigauto/pigauto_workflow_mixed.html)
- **Architecture notes**: [`CLAUDE.md`](CLAUDE.md)

## Citation

Nakagawa S (2026). *pigauto: Phylogenetic Imputation via Graph
Autoencoder*. R package version 0.9.1.9009.
