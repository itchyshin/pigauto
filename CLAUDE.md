# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`pigauto` is an R package (version 0.7.0) for phylogenetic trait imputation. It fits a gated ensemble of a phylogenetic baseline and an attention-based graph neural network correction. For continuous/count/ordinal traits the baseline is Brownian motion (via an internal conditional-MVN implementation using the phylogenetic correlation matrix `R = cov2cor(vcv(tree))`; see `R/bm_internal.R`); for binary/categorical traits it is phylogenetic label propagation. Optional environmental covariates are threaded through the GNN with gated safety. Prediction is the per-trait blend `(1 - r_cal) * baseline + r_cal * delta_GNN`, with `r_cal` calibrated on a held-out validation split. See `README.md` for the user-facing API; this file documents the internals.

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

`evaluate_imputation()` returns a wide-format data.frame. Most types populate `rmse` / `pearson_r` / `mae` / `accuracy` / `brier` / `spearman_rho`. `multi_proportion` additionally populates three compositional-specific columns: `aitchison` (Euclidean distance in CLR space — the natural compositional metric), `rmse_clr` (RMSE on z-scored CLR latent, comparable to continuous RMSE), and `simplex_mae` (mean abs error on the decoded proportions). Other types carry `NA` in these columns.

The pipeline functions are for fine-grained control and for writing benchmarks. For the multiple-imputation → downstream-inference workflow, see `R/multi_impute.R` → `R/with_imputations.R` → `R/pool_mi.R` (Rubin 1987; Barnard & Rubin 1999; Nakagawa & Freckleton 2008, 2011). User-facing tutorial: `inst/doc/pigauto_workflow_mixed.html` (source: `script/make_workflow_mixed_html.R`) is the self-contained walk-through on the bundled AVONET 300 mixed-type dataset, covering all three analysis paths (Path A = pigauto+glmmTMB+Rubin; Path B = pigauto+MCMCglmm+posterior concatenation; Path C = BACE integrated). There is no continuous-only sibling tutorial — it was removed in the v0.3.0 docs cleanup because the mixed-type walk-through strictly covers it.

### Trait-type handling

All eight types live in a single latent matrix:

| R class    | pigauto type | Latent cols | Loss          | Baseline                                  |
|------------|-------------|-------------|---------------|-------------------------------------------|
| numeric    | continuous  | 1 (z-score, optional log) | MSE | Internal phylogenetic BM                  |
| integer    | count       | 1 (log1p + z) | MSE          | Internal phylogenetic BM                  |
| ordered    | ordinal     | 1 (integer + z) | MSE        | Internal phylogenetic BM                  |
| factor(2)  | binary      | 1 (0/1)     | BCE           | Phylogenetic label propagation            |
| factor(>2) | categorical | K (one-hot) | cross-entropy | Phylogenetic label propagation            |
| numeric(0–1) | proportion | 1 (logit + z) | MSE       | Internal phylogenetic BM                  |
| integer(ZI) | zi_count   | 2 (gate 0/1 + log1p-z) | BCE + cond. MSE | LP (gate) + BM (magnitude) |
| K numeric cols summing to 1 | multi_proportion | K (CLR + per-component z) | MSE (CLR) | Per-component BM on CLR |

Proportion and zi_count require explicit `trait_types` overrides — they are not auto-detected. ZI count uses two latent columns: column 1 is a binary gate (0=zero, 1=non-zero), column 2 is the log1p-z magnitude (NA when observed value is zero). Expected value prediction: `E[X] = P(non-zero) * E[count | non-zero]`.

**multi_proportion is declared differently from the others.** Not via `trait_types`, but via a separate `multi_proportion_groups` argument on `preprocess_traits()`, `impute()`, and `multi_impute()`:

```r
impute(df, tree, multi_proportion_groups = list(
  colour = c("black", "blue", "red", ..., "yellow")))
```

Zero handling: small `epsilon = 1e-6` added before `log()` so `log(0)` is safe; rows are re-normalised so the composition is preserved. Decoding: CLR values projected to sum-zero, then softmax → simplex. Row-level masking: the entire K-column composition is either fully observed or fully missing (you can't observe 3 of K components), so DAE corruption and val/test splits operate on whole rows. The K input columns are excluded from per-column type detection and replaced by ONE group entry in `trait_map` (key = group name, `input_cols = c(...)`, `levels = c(...)`, `mean`/`sd` as K-vectors).

`trait_map` (list of descriptors) is the single source of truth for type, latent-column range, levels, mean, sd, and log-transform flag. `preprocess_traits()` builds it; `fit_baseline()`, `fit_pigauto()`, `predict()`, and `evaluate()` all use it. **When touching encoding or decoding, change the encoder and decoder in sync and update every consumer of `trait_map`**.

Critical: because categorical traits expand to K latent columns, `length(trait_names)` and `ncol(X_scaled)` differ. Always use `colnames(data$X_scaled)` (latent names) when dealing with latent matrices, not `data$trait_names`.

### Multi-observation per species

When `preprocess_traits(traits, tree, species_col = "species")` is used, `data$X_scaled` has `n_obs` rows and `data$obs_to_species` maps each row to a species index. The GNN operates at species level: `scatter_mean(h, obs_to_species, n_species)` aggregates before message passing, and `index_select(1L, obs_to_species)` broadcasts back after. The baseline `mu` is species-level and is expanded to observation level in `fit_pigauto` before blending.

**Observation-level refinement** (v0.6.0+): when `covariates` are supplied in multi-obs mode, the model creates an `obs_refine` MLP (`nn_sequential(linear(hidden + n_user_cov, hidden), relu, linear(hidden, hidden))`) that re-injects user covariates after the species-level broadcast via a residual connection: `h = h + obs_refine(cat(h, user_covs))`. This allows different observations of the same species to receive different predictions based on their covariate context (e.g., CTmax at 20°C vs 30°C acclimation). The refinement is controlled by `n_user_cov` in model_config; when 0, no refinement MLP is created. The bundled `ctmax_sim` dataset + `tree300` exercises this code path. Benchmark: `script/bench_multi_obs.R`.

### Liability encoding (Phase 1 of Level C — unreleased / v0.8.0-alpha)

`R/liability.R` formalises the liability interpretation of each trait type
as the foundation for a future joint multivariate-BM baseline. Continuous
types ARE their liability (identity / log1p / logit / CLR on the
preprocessed latent scale). Binary uses a 1D threshold at 0. Ordinal uses
K-1 thresholds. Categorical uses K liabilities with an argmax constraint
(plug-in approximation; exact Gibbs in Phase 6 EM).

`estep_liability(tm, observed, mu_prior, sd_prior)` is the single dispatcher
that all Phase 2+ code should call to compute per-cell posterior mean/var
of a liability given an observation. One non-obvious detail: for `ordinal`,
`observed` comes z-scored from `X_scaled`, so the dispatcher un-z-scores
before recovering the integer class index. Do not bypass this step when
calling `estep_liability_ordinal()` directly.

Nothing in `fit_baseline()` / `fit_pigauto()` consumes this module yet —
Phase 2 will wire the joint multivariate-BM baseline on top.

### Post-training (`fit_pigauto.R`)

After the training loop:
1. Val-set grid search picks a per-trait calibrated gate (continuous: MSE; discrete: 0-1 loss with absolute cell floor and half-A/half-B cross-check).
2. Conformal scores are computed from validation residuals (continuous/count/ordinal) → stored in the fit and used at prediction time for 95% coverage intervals.

## Repository layout

- `R/` — package source. Everything with an `@export` tag is user-facing.
- `tests/testthat/` — testthat 3rd edition, 558 tests total. One test file per broad area: preprocess, graph, masking, fit-predict, mixed-types, multi-impute, multi-proportion, new-features.
- `BACE/` — a **separate, self-contained R package** (Bayesian phylogenetic imputation via MCMCglmm) kept in-tree as a reference implementation and comparison baseline. It has its own `R/`, `tests/`, `vignettes/`, and `DESCRIPTION`. `Grep` and `Glob` results for generic terms (`impute`, `phylo`, `trait`) will include BACE files — always check the path prefix. Pigauto wraps BACE only in `R/fit_baseline_bace.R`. BACE is `Suggests:`-only, `^BACE$` is in `.Rbuildignore`, and BACE's own tests are not part of pigauto's test suite. Do not modify BACE as part of pigauto work.
- `script/` — benchmark drivers, logs, and HTML/RDS outputs. Ignored by `R CMD build`. Key entries: `validate_avonet_full.{R,log,md,rds}` (full-scale validation), `bench_scaling_v031.{R,log,rds}` (scaling benchmark), `bench_avonet_missingness.{R,rds,md}` + `make_avonet_missingness_html.R` (missingness sweep), and the per-type benchmark suite: `bench_{continuous,binary,ordinal,count,categorical,proportion,zi_count,multi_proportion,missingness_mechanism}.R` (drivers) + `make_bench_*_html.R` (HTML generators). Each driver outputs `.rds` + `.md`; each HTML generator outputs to both `script/` and `pkgdown/assets/dev/`. Anything named `bench_v2.*`, `bench_v3.*`, `bench_v4.*`, or `benchmark_*` is a stale snapshot from earlier phases — do not treat them as reference implementations.
- `dev/` — scratch experiments. Ignored by `R CMD build`.
- `avonet/`, `data/`, `data-raw/` — the bundled AVONET 300-species dataset and its build scripts.

## Uncertainty quantification design

pigauto uses **three distinct uncertainty mechanisms** — do not conflate them:

### 1. Baseline SE (analytic, BM conditional MVN)
Source: `R/bm_internal.R` → `bm_impute_col()`.  
For continuous/count/ordinal/proportion traits: standard conditional-MVN formula,
`SE = sqrt(σ² (1 - h_i))` where `h_i = diag(R_mo R_oo⁻¹ R_om)`.  
**Validity**: exact under Brownian motion; model-dependent (assumes BM is the true process).  
**Used in**: `pred$se` for these trait types (after delta-method back-transformation).

### 2. Conformal prediction intervals (distribution-free, empirically calibrated)
Source: `R/fit_helpers.R` → `compute_conformal_scores()`.  
After training, on the held-out validation set: `score_j = quantile(|truth - blended_pred|, ⌈(1−α)(1+1/n)⌉/n)`.  
**Validity**: split conformal guarantee — exactly ≥95% marginal coverage regardless of model assumptions or trait distribution. No Gaussianity needed.  
**Used in**: `pred$conformal_lower`, `pred$conformal_upper`. This is the primary 95% CI.

### 3. Multiple-imputation draws (`multi_impute()`)
Two methods selectable via `draws_method`:
- **`"conformal"` (default)**: Single pass; missing cells sampled from N(μ, conformal_score/1.96) on the appropriate transformed scale. Falls back to BM-SE-based Normal sampling when conformal scores are missing, and to Bernoulli/Categorical for discrete types. Preferred because conformal scores are calibrated against actual held-out residuals regardless of gate value.
- **`"mc_dropout"`**: M GNN forward passes in training mode (dropout active). Each imputation `m` draws `t_BM_draw ~ N(BM_mu, BM_se)` on the latent scale (held fixed for all refine steps of that imputation), then blends `pred = (1 - r_cal) * t_BM_draw + r_cal * GNN_dropout(t_BM_draw)`. When `r_cal = 0` (gate closed — BM dominates): `pred = t_BM_draw` → between-imputation variance = BM posterior variance, non-zero ✓. When `r_cal > 0`: both BM draws and GNN dropout contribute variance. `BM_se = 0` for observed cells so they are never perturbed. **Note on conservatism**: BM-draw MI is wider than conformal MI (AVONET300: Mass MC SD ≈ 290 vs conformal/1.96 ≈ 23) because BM SE reflects prior uncertainty while conformal reflects actual prediction error. For downstream Rubin's rules, conformal is better calibrated. Implementation: `predict_pigauto.R` lines 168–211.

### 4. `pred$se` for discrete types — uncertainty scores, not SEs
Binary: `min(p, 1-p)` — probability of being wrong (0 = certain, 0.5 = maximally uncertain).  
Categorical: `1 - max(p_k)` — margin from certainty (0 = certain, (K-1)/K = maximally uncertain).  
**These are not standard errors in the Gaussian sense.** Do not use them in Rubin's rules arithmetic. They are uncertainty scores for ranking/reporting. For MI draws, binary/categorical always use Bernoulli/Categorical sampling from the probability vector, not Normal draws.

### What NOT to do
- Do not use `pred$se` for binary/categorical as if it were a Normal SE.
- Do not back-calculate a "±1.96×SE" interval for discrete types.
- Do not use BM SE alone as a 95% CI — use the conformal interval instead (it is wider and better calibrated when the GNN adds prediction error beyond BM).

## Tree uncertainty — two-step workflow

Tree uncertainty enters the analysis at TWO distinct places. Do not conflate them.

**Step 1 — imputation (pigauto's job).** `multi_impute_trees(traits, trees = trees300, m_per_tree = 5L)` runs a full pigauto fit once per posterior tree, producing `T × m_per_tree` completed datasets. Each dataset is conditional on a specific tree; `mi$tree_index[i]` records which tree produced dataset `i`. No caching across trees is possible — the tree IS the model (graph, baseline, GNN weights all change per tree).

**Step 2 — analysis (user's responsibility).** For each completed dataset, refit the downstream comparative model using the SAME tree that produced it, then pool all T × M fits via `pool_mi()`:

```r
fits <- Map(function(dat, t_idx) {
  dat$species <- rownames(dat)
  nlme::gls(y ~ x, correlation = ape::corBrownian(phy = trees[[t_idx]], form = ~species),
            data = dat, method = "ML")
}, mi$datasets, mi$tree_index)
pool_mi(fits)
```

This is the Nakagawa & de Villemereuil (2019, *Syst. Biol.* 68:632–641) algorithm — trees as missing data, pooled via Rubin's rules. pigauto handles step 1; step 2 stays with the user because the downstream model is user-chosen.

**Compute cost is linear in T.** Rough budget: n=300 × T=50 ≈ 25–50 min, n=5,000 × T=50 ≈ 4–8 hr, n=10,000 × T=50 ≈ 17–33 hr. At n ≥ 5,000 reduce T to 10–20 (the 2019 paper's "relative efficiency" index typically converges well before T=50) or parallelise across machines.

**Future work (not built): `share_gnn = TRUE`.** Train the GNN once on the MCC tree and reuse it across posterior trees while recomputing only the cheap baseline per tree. ~10–15× speedup for large n. Gate safety means worst case is "GNN is useless → per-tree BM fallback" — still tree-uncertainty-aware.

## Non-obvious gotchas

### `.Rbuildignore` is wide

`BACE/`, `script/`, `dev/`, `useful/`, `avonet/`, `checkpoints*`, `data-raw/`, `NEWS.md`, and a few others are excluded from the built package. New files placed in those directories will not ship. If you need something in the installed package, put it in `R/`, `inst/`, `data/`, `man/`, `vignettes/`, or `tests/`.

### `k_eigen = "auto"` scales with tree size

`build_phylo_graph()` picks `k_eigen` adaptively: 4 for ~30 tips, 15 for 300, up to 32 for 640+. When writing tests with tiny trees (≤40 tips), pass a small explicit `k_eigen` or let `auto` handle it, but expect different numeric results across tree sizes.

### Gate initialisation depends on trait type

Continuous/count/ordinal/proportion columns initialise `res_raw` at `-1` (effective gate ≈ 0.135 × gate_cap). Binary/categorical columns initialise near 0 (fully closed). ZI count: column 1 (gate) initialises like binary (near 0), column 2 (magnitude) like continuous (-1). This is deliberate: discrete traits should rely entirely on phylogenetic label propagation until the GNN earns the right to contribute. Preserve this when editing `ResidualPhyloDAE$initialize`.

### Multi-obs code path in `fit_pigauto.R`

Several spots switch behaviour on `multi_obs`: baseline expansion (`MU <- MU_species[obs_to_sp, ]`), the `torch_tensor(obs_to_sp, ..., torch_long())` creation, and the `obs_to_species` argument passed to the model forward. If you add a new input tensor that is computed at species level, remember to expand it, and vice versa for observation-level tensors going into the GNN.

### pkgdown GitHub Actions on pull requests

`.github/workflows/pkgdown.yaml` has a job-level `if: github.event_name != 'pull_request'` that skips the entire pkgdown job on PR events. This is deliberate. The job attaches to the `github-pages` environment at job level (for the deploy step), and that environment has a protection rule that only allows `main` to deploy. The environment check fires BEFORE any `step: if` filter, so PR runs would always fail at the env gate even though the Upload/Deploy steps were themselves already conditionally skipped on PRs. Skipping the whole job on PRs eliminates the spurious red X on PR checks. Do not remove the job-level `if` without splitting into a separate build-on-PR job that does not attach to the `github-pages` environment.

## Host-specific notes (optional)

On the primary author's machine, two persistent memory notes live under `~/.claude/projects/-Users-z3437171-Dropbox-Github-Local-pigauto/memory/`: `user_profile.md` (author priorities) and `project_bace.md` (BACE internals). They are **not** portable — ignore this section if the path does not exist on the current host.
