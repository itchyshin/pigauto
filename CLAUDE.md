# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`pigauto` is an R package (version 0.5.0) for phylogenetic trait imputation. It fits a gated ensemble of a phylogenetic baseline and an attention-based graph neural network correction. For continuous/count/ordinal traits the baseline is Brownian motion (via an internal conditional-MVN implementation using the phylogenetic correlation matrix `R = cov2cor(vcv(tree))`; see `R/bm_internal.R`); for binary/categorical traits it is phylogenetic label propagation. Optional environmental covariates are threaded through the GNN with gated safety. Prediction is the per-trait blend `(1 - r_cal) * baseline + r_cal * delta_GNN`, with `r_cal` calibrated on a held-out validation split. See `README.md` for the user-facing API; this file documents the internals.

**A note on "residual"**: the internal torch class is named `ResidualPhyloDAE` because its GNN layers use ResNet-style residual skip connections. The GNN output `delta` is **not** a statistical residual `y - baseline` — it is a full per-cell prediction trained end-to-end via type-appropriate loss (MSE / BCE / cross-entropy) on the blend, not on `y - baseline`. Do not re-introduce user-facing prose that describes the GNN as "learning a residual from the baseline".

## Common commands

```r
# First-time setup (after a fresh clone or torch upgrade)
devtools::install()
torch::install_torch()

# Normal dev loop
devtools::load_all()
devtools::test()                                  # all tests
devtools::test(filter = "fit-predict")            # one file
testthat::test_file("tests/testthat/test-mixed-types.R")  # single file, no reload
devtools::check()                                 # R CMD check before commits
devtools::document()                              # regenerate NAMESPACE/man after roxygen edits
```

Tests live in `tests/testthat/` (testthat 3rd edition). There is no linter or formatter configured — match surrounding style.

## Architecture

### The core formulation

Every prediction is

    pred = (1 - r_cal) * BM_baseline + r_cal * GNN_delta

where `r_cal` is a **per-latent-column gate** (vector of length `p_latent`) that is:

1. **Learned** during training as `sigmoid(res_raw) * gate_cap` with an L2 penalty `lambda_gate` pushing it toward zero.
2. **Re-calibrated** on the validation set after training via grid search per trait (`fit_pigauto.R`, look for the v4 calibration block). For discrete traits the calibration uses 0-1 loss with an absolute cell floor (min 2 cells on half A, 1 on half B) and a split-validation cross-check to prevent the GNN from degrading baseline accuracy.
3. Enforced so that `r_cal = 0` is always a valid fallback. When the baseline is already optimal the calibrated gate closes and the GNN becomes a no-op. **Do not remove this safety without understanding why it exists** — without it the GNN degrades discrete accuracy on high-phylo-signal data.

### The model (`R/model_residual_dae.R`)

`ResidualPhyloDAE` is an internal `torch::nn_module`. The "Residual" in the class name refers to the ResNet-style skip connections inside the GNN layers, not to a statistical residual.

- **Input**: `x` (`n_obs × p_latent` with corrupted cells replaced by a learnable mask token), `coords` (`n_species × k_eigen` spectral features from the Laplacian), `covs` (`n_obs × cov_dim` — currently baseline mean + NA-mask indicator; no user covariates yet), `adj` (`n_species × n_species` symmetric-normalised Gaussian kernel).
- **Encoder**: 2 linear layers + ReLU + dropout, concatenating `x`, `coords`, `covs`.
- **Message passing**: `n_gnn_layers` (default 2) graph layers. `use_attention = TRUE` (default) uses scaled dot-product attention with a learnable log-adjacency bias, so the model starts close to the phylogenetic prior and learns deviations. Each layer has its own learnable alpha gate, LayerNorm, and ResNet-style skip connection.
- **Decoder**: 2 linear layers → `delta` (`n_obs × p_latent`). Note: `delta` is a full per-cell prediction, not `y - baseline`.
- **Output**: `(1 - rs) * baseline_mu + rs * delta`, where `rs` is the per-column sigmoid gate bounded in `(0, gate_cap)`.

### Data flow

**Primary user entry point**: `impute(traits, tree)` in `R/impute.R`. Everything below is what it calls internally.

```
preprocess_traits()  →  build_phylo_graph()  →  fit_baseline()  →  fit_pigauto()  →  predict()  →  evaluate()
   (trait_map)          (coords, adj)          (mu, se)            (model + cal)      (pigauto_pred)
```

`impute()` chains all of these in one call and returns a `pigauto_result` containing `completed` (the user's dataframe with NAs filled), `imputed_mask`, `prediction`, and the intermediate objects.

### Key S3 classes

| Class | Produced by | Contains |
|---|---|---|
| `pigauto_data` | `preprocess_traits()` | `X_scaled`, `trait_map`, `obs_to_species`, species/trait names |
| `pigauto_fit` | `fit_pigauto()` | `model_state`, `model_config`, calibrated gates, conformal scores, splits |
| `pigauto_pred` | `predict.pigauto_fit()` | `imputed`, `imputed_latent`, `se`, `probabilities`, conformal intervals |
| `pigauto_result` | `impute()` | wraps the above plus `completed` and `imputed_mask` |
| `pigauto_benchmark` | `simulate_benchmark()` | multi-scenario results data.frame |
| `pigauto_cv` | `cross_validate()` | per-fold metrics |
| `pigauto_mi` | `multi_impute()` | `datasets` (list of M), `m`, `pooled_point`, `se`, `imputed_mask`, `fit`, `tree` |
| `pigauto_mi_fits` | `with_imputations()` | list of M downstream fits (or `pigauto_mi_error` for failed draws) |
| `pigauto_pooled` | `pool_mi()` | tidy data.frame with Rubin's-rules pooled coefficients (`estimate`, `std.error`, `df`, `fmi`, `riv`, ...) |

The pipeline functions are for fine-grained control and for writing benchmarks. For the multiple-imputation → downstream-inference workflow, see `R/multi_impute.R` → `R/with_imputations.R` → `R/pool_mi.R` (Rubin 1987; Barnard & Rubin 1999; Nakagawa & Freckleton 2008, 2011). User-facing tutorial: `inst/doc/pigauto_workflow_mixed.html` (source: `script/make_workflow_mixed_html.R`) is the self-contained walk-through on the bundled AVONET 300 mixed-type dataset, covering all three analysis paths (Path A = pigauto+glmmTMB+Rubin; Path B = pigauto+MCMCglmm+posterior concatenation; Path C = BACE integrated). There is no continuous-only sibling tutorial — it was removed in the v0.3.0 docs cleanup because the mixed-type walk-through strictly covers it.

### Trait-type handling

All seven types live in a single latent matrix:

| R class    | pigauto type | Latent cols | Loss          | Baseline                                  |
|------------|-------------|-------------|---------------|-------------------------------------------|
| numeric    | continuous  | 1 (z-score, optional log) | MSE | Internal phylogenetic BM                  |
| integer    | count       | 1 (log1p + z) | MSE          | Internal phylogenetic BM                  |
| ordered    | ordinal     | 1 (integer + z) | MSE        | Internal phylogenetic BM                  |
| factor(2)  | binary      | 1 (0/1)     | BCE           | Phylogenetic label propagation            |
| factor(>2) | categorical | K (one-hot) | cross-entropy | Phylogenetic label propagation            |
| numeric(0–1) | proportion | 1 (logit + z) | MSE       | Internal phylogenetic BM                  |
| integer(ZI) | zi_count   | 2 (gate 0/1 + log1p-z) | BCE + cond. MSE | LP (gate) + BM (magnitude) |

Proportion and zi_count require explicit `trait_types` overrides — they are not auto-detected. ZI count uses two latent columns: column 1 is a binary gate (0=zero, 1=non-zero), column 2 is the log1p-z magnitude (NA when observed value is zero). Expected value prediction: `E[X] = P(non-zero) * E[count | non-zero]`.

`trait_map` (list of descriptors) is the single source of truth for type, latent-column range, levels, mean, sd, and log-transform flag. `preprocess_traits()` builds it; `fit_baseline()`, `fit_pigauto()`, `predict()`, and `evaluate()` all use it. **When touching encoding or decoding, change the encoder and decoder in sync and update every consumer of `trait_map`**.

Critical: because categorical traits expand to K latent columns, `length(trait_names)` and `ncol(X_scaled)` differ. Always use `colnames(data$X_scaled)` (latent names) when dealing with latent matrices, not `data$trait_names`.

### Multi-observation per species

When `preprocess_traits(traits, tree, species_col = "species")` is used, `data$X_scaled` has `n_obs` rows and `data$obs_to_species` maps each row to a species index. The GNN operates at species level: `scatter_mean(h, obs_to_species, n_species)` aggregates before message passing, and `index_select(1L, obs_to_species)` broadcasts back after. The baseline `mu` is species-level and is expanded to observation level in `fit_pigauto` before blending. **The infrastructure is there but end-to-end behaviour is not exercised by a benchmark** — add one before claiming it works on real multi-obs data.

### Post-training (`fit_pigauto.R`)

After the training loop:
1. Val-set grid search picks a per-trait calibrated gate (continuous: MSE; discrete: 0-1 loss with absolute cell floor and half-A/half-B cross-check).
2. Conformal scores are computed from validation residuals (continuous/count/ordinal) → stored in the fit and used at prediction time for 95% coverage intervals.

## Repository layout

- `R/` — package source (~7100 lines). Everything with an `@export` tag is user-facing.
- `tests/testthat/` — testthat 3rd edition. One test file per broad area: preprocess, graph, masking, fit-predict, mixed-types, multi-impute, new-features.
- `BACE/` — a **separate, self-contained R package** (Bayesian phylogenetic imputation via MCMCglmm) kept in-tree as a reference implementation and comparison baseline. It has its own `R/`, `tests/`, `vignettes/`, and `DESCRIPTION`. `Grep` and `Glob` results for generic terms (`impute`, `phylo`, `trait`) will include BACE files — always check the path prefix. Pigauto wraps BACE only in `R/fit_baseline_bace.R`. BACE is `Suggests:`-only, `^BACE$` is in `.Rbuildignore`, and BACE's own tests are not part of pigauto's test suite. Do not modify BACE as part of pigauto work.
- `script/` — benchmark drivers, logs, and HTML/RDS outputs. Ignored by `R CMD build`. Key entries: `validate_avonet_full.{R,log,md,rds}` (full-scale validation), `bench_scaling_v031.{R,log,rds}` (scaling benchmark), `bench_avonet_missingness.{R,rds,md}` + `make_avonet_missingness_html.R` (missingness sweep), and the per-type benchmark suite: `bench_{continuous,binary,ordinal,count,categorical,proportion,zi_count,missingness_mechanism}.R` (drivers) + `make_bench_*_html.R` (HTML generators). Each driver outputs `.rds` + `.md`; each HTML generator outputs to both `script/` and `pkgdown/assets/dev/`. Anything named `bench_v2.*`, `bench_v3.*`, `bench_v4.*`, or `benchmark_*` is a stale snapshot from earlier phases — do not treat them as reference implementations.
- `dev/` — scratch experiments. Ignored by `R CMD build`.
- `avonet/`, `data/`, `data-raw/` — the bundled AVONET 300-species dataset and its build scripts.

## Non-obvious gotchas

### `.Rbuildignore` is wide

`BACE/`, `script/`, `dev/`, `useful/`, `avonet/`, `checkpoints*`, `data-raw/`, `NEWS.md`, and a few others are excluded from the built package. New files placed in those directories will not ship. If you need something in the installed package, put it in `R/`, `inst/`, `data/`, `man/`, `vignettes/`, or `tests/`.

### `k_eigen = "auto"` scales with tree size

`build_phylo_graph()` picks `k_eigen` adaptively: 4 for ~30 tips, 15 for 300, up to 32 for 640+. When writing tests with tiny trees (≤40 tips), pass a small explicit `k_eigen` or let `auto` handle it, but expect different numeric results across tree sizes.

### Gate initialisation depends on trait type

Continuous/count/ordinal/proportion columns initialise `res_raw` at `-1` (effective gate ≈ 0.135 × gate_cap). Binary/categorical columns initialise near 0 (fully closed). ZI count: column 1 (gate) initialises like binary (near 0), column 2 (magnitude) like continuous (-1). This is deliberate: discrete traits should rely entirely on phylogenetic label propagation until the GNN earns the right to contribute. Preserve this when editing `ResidualPhyloDAE$initialize`.

### Multi-obs code path in `fit_pigauto.R`

Several spots switch behaviour on `multi_obs`: baseline expansion (`MU <- MU_species[obs_to_sp, ]`), the `torch_tensor(obs_to_sp, ..., torch_long())` creation, and the `obs_to_species` argument passed to the model forward. If you add a new input tensor that is computed at species level, remember to expand it, and vice versa for observation-level tensors going into the GNN.

## Host-specific notes (optional)

On the primary author's machine, two persistent memory notes live under `~/.claude/projects/-Users-z3437171-Dropbox-Github-Local-pigauto/memory/`: `user_profile.md` (author priorities) and `project_bace.md` (BACE internals). They are **not** portable — ignore this section if the path does not exist on the current host.
