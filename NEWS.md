# pigauto unreleased (v0.9.0-alpha — Graph Transformer GNN + multi-obs Level-C)

## Phase 10: Level-C baselines activate for multi-observation datasets

- `fit_baseline()` no longer gates out multi-obs data from the joint MVN
  (Phase 2), threshold-joint (Phase 3), and OVR categorical (Phase 6)
  paths. Multi-obs input is aggregated to species level inside each
  Level-C helper (mean for continuous-family, threshold-at-0.5 for
  binary/ZI-gate, argmax-one-hot for categorical) and then processed
  identically to single-obs data.
- New internal helper: `aggregate_to_species()` in
  `R/joint_threshold_baseline.R`. Mask-then-aggregate ordering prevents
  val/test leakage; single-obs data passes through unchanged with
  masking preserved for back-compat.
- All three Level-C entry points (`fit_joint_mvn_baseline`,
  `fit_joint_threshold_baseline`, `fit_ovr_categorical_fits`) accept
  multi-obs via aggregation at their top. `build_liability_matrix`
  likewise.
- Four `!multi_obs` guards removed from `fit_baseline()` dispatch
  (threshold-joint condition, continuous-joint condition, OVR loop,
  and the associated `stopifnot` in `build_liability_matrix`).
- Aggregation caveat: binary/categorical aggregation is lossy. A
  species with mixed observations (6/10 class-1) becomes a single
  class-1 observation at species level. For datasets with strong
  within-species variability in discrete traits, prefer single-obs
  mode with a user-computed species-level summary.

---

## Phase 9: Graph Transformer backbone for the GNN

- New internal module `R/graph_transformer_block.R`. `GraphTransformerBlock`
  is a standard transformer-encoder block adapted for graph inputs:
  multi-head attention with a learnable per-head `log(adjacency + eps)`
  bias (the phylo prior carries through to each head), FFN with
  configurable width multiplier, pre-norm, residual skips. The FFN
  output projection initialises to zero so each block is near-identity
  at step 0 — the stacked composition is still near-identity, preserving
  pigauto's gate-closed-at-init safety property.

- `ResidualPhyloDAE` now accepts `use_transformer_blocks = TRUE` (new
  default), `n_heads = 4L`, and `ffn_mult = 4L`. When the flag is on,
  the 2-layer single-head attention stack is replaced by `n_gnn_layers`
  instances of `GraphTransformerBlock`. The legacy architecture is fully
  preserved behind `use_transformer_blocks = FALSE` and reproduces
  pre-Phase-9 numbers exactly.

- `fit_pigauto()` gains the same three arguments and threads them into
  the model constructor. `predict.pigauto_fit()` reads
  `use_transformer_blocks` from `model_config` with an implicit
  `FALSE` fallback, so pre-Phase-9 saved fits reconstruct via the legacy
  path without any user-visible breakage.

- Bench (`script/bench_discriminative_phase9.R`): on 4 mixed-signal
  simulation scenarios (n=200, reps=2, epochs=100), transformer and
  legacy GNN produce identical predictions in high-signal scenarios
  (gate closes to zero — the architecture below the gate is irrelevant)
  and the transformer shows small RMSE lift in moderate/low signal
  scenarios where the gate stays partly open. Runtime ~25–30% slower
  per fit (not 2–4× as initially feared).

- Larger transformer gains likely need bigger data plus self-supervised
  pretraining (strategic roadmap option B). Phase 9 ships the
  architecture and the plumbing; the ceiling lift is future work.

---

# pigauto unreleased (v0.8.0-alpha — Level C baseline)

## Infrastructure

- New internal module `R/liability.R` with the liability-encoding contract for
  all 8 trait types. `liability_info()` returns the per-type schema
  (`n_liability`, `kind`, `thresholds`). `estep_liability()` computes the
  posterior mean and variance of a Gaussian liability given an observed value,
  dispatching over continuous (point-mass), binary (truncated Gaussian),
  ordinal (interval-truncated), and categorical (plug-in argmax boost).
  Exact Gibbs E-step for categorical deferred to Phase 6 EM.

- The ordinal dispatcher correctly un-z-scores observations from `X_scaled`
  before recovering the integer class index, and clamps to `[1..K]`. A
  regression test (`test-liability.R`) covers this path for both mid-range
  and edge-of-distribution classes.

## Internal (not user-facing yet)

- These functions are the foundation for the Phase 2+ joint multivariate
  baseline. No behavioural change in this release — nothing new is called
  from `fit_baseline()` or `fit_pigauto()` yet.

- New `R/joint_mvn_baseline.R`: when `Rphylopars` is installed and the data
  has 2+ BM-eligible latent columns, `fit_baseline()` delegates to a joint
  multivariate-BM baseline via `Rphylopars::phylopars()`. Captures cross-trait
  phylogenetic correlation that the per-column path missed. Falls back
  gracefully to per-column BM when `Rphylopars` is unavailable, when any
  `multi_proportion` trait is present, or in multi-obs mode. Single-trait
  case reproduces current output within 1e-3 (tested). Benchmark on
  correlated simulated BM data: **33.7% RMSE lift** vs per-column BM
  (`script/bench_joint_baseline.R`, n=100 tips × 10 reps).

### Level C Phase 3 — BINARY via threshold joint baseline

- `fit_baseline()` now dispatches to `fit_joint_threshold_baseline()`
  when BINARY cols are present alongside continuous cols and
  `Rphylopars` is available. Observed binary cells are transformed to
  truncated-Gaussian posterior means via `estep_liability_binary`
  (Phase 1), combined with continuous latents, and jointly fit by
  `Rphylopars::phylopars()`. Binary output is decoded back to
  `logit(P(y=1))` so downstream BCE and gate math are unchanged.
- New file: `R/joint_threshold_baseline.R`
  (`build_liability_matrix`, `fit_joint_threshold_baseline`,
  `decode_binary_liability`).
- New test file: `tests/testthat/test-joint-threshold-baseline.R`
  (6 new test blocks, 56 expectations covering builder, fitter,
  decoder, dispatcher, graceful fallback, and the <2-obs edge case).
- New benchmark: `script/bench_binary_joint.R` — A/B on correlated
  continuous+binary BM simulations. **Result: 3/3 scenarios with lift
  >= 2pp**, peaking at 16pp at rho = 0.6 (moderate cross-trait
  correlation where LP has the least same-trait phylo signal to
  exploit).
- Graceful fallback: when `Rphylopars` is unavailable the LP baseline
  runs unchanged (verified by `test-joint-threshold-baseline.R`).
  Per-column fallback when a binary col has <2 observations.
- Ordinal threshold decoding is deferred to a later phase (Phase 6 EM).
  Ordinal still rides in the Phase-2 MVN via its z-scored integer
  encoding.

### Level C Phase 6 — CATEGORICAL via OVR K-fits (stable, default)

- `fit_baseline()` now dispatches categorical traits through K independent
  `fit_joint_threshold_baseline` calls — one per class, treating each as
  a "is_class_k vs rest" binary threshold problem. This is the same OVR
  strategy BACE uses, which lifted their AVONET categorical accuracy
  from ~42% to 72%.
- Each individual fit has only 1 categorical-related col, so Rphylopars
  stays well-conditioned — avoiding the rank-(K-1) numerical instability
  that blocked the single-fit approach in Phase 4.
- New file: `R/ovr_categorical.R` (`fit_ovr_categorical_fits()` and
  `decode_ovr_categorical()`).
- `include_categorical` and `cat_encoding` arguments removed (the
  single-fit joint_K / OVR-in-one-matrix paths are retired — the
  K-independent-fits path is strictly better and stable).
- Cost: K times the phylopars cost per categorical trait. For n <= 1000
  species and K <= 10, this is still fast relative to GNN training.
- Benchmark: OVR K-fits beats LP on correlated simulated data.

### Level C Phase 5 — ZI gate joins the threshold-joint path

- Zero-inflated count traits have two latent cols: a binary gate
  (observed-is-zero vs observed-is-nonzero) and a continuous magnitude
  (log1p-z). The magnitude col already rode Phase 2's joint MVN
  baseline (as a BM-eligible continuous). Phase 5 routes the gate col
  through the Phase 3 threshold-joint path alongside binary traits,
  so ZI gates now benefit from cross-trait correlation with other
  continuous/binary traits just like regular binaries.
- `build_liability_matrix()` now includes the ZI gate col as a
  binary-like liability (treated identically to binary by the
  truncated-Gaussian E-step).
- `fit_baseline()` dispatcher collects gate cols into `binary_cols`
  for threshold-joint dispatch; the ZI gate LP loop skips cols
  already handled.
- Count, proportion, ordinal already ride Phase 2's joint MVN via
  the `bm_cols` set — no Phase 5 changes needed for them. They were
  already cross-trait correlation-aware.

**Deferred to Phase 6**: multi_proportion (K CLR cols hit the same
phylopars rank-deficiency instability as categorical K-col groups).
multi_proportion continues to run per-column BM via its existing path.

# pigauto 0.7.0

## New trait type: multi_proportion (compositional data)

- **`multi_proportion`**: a new trait type for compositional data — K columns
  per row that sum to 1 (e.g. plumage-colour proportions, diet composition,
  microbiome relative abundances, allele frequencies). Declared via a new
  `multi_proportion_groups = list(colour = c("black", ..., "yellow"))`
  argument to `impute()`, `preprocess_traits()`, and `multi_impute()`.

- **Encoding**: centred log-ratio (CLR) + per-component z-score. Small
  epsilon (`1e-6`) is applied to zero components before log to keep the
  transform safe; rows are re-normalised so the composition is preserved.

- **Baseline**: K independent Brownian-motion imputations on the CLR
  columns. BM is well-defined on CLR space (Euclidean), unlike the raw
  simplex.

- **GNN output**: K logits projected to sum-zero, then softmax at decode
  time — guaranteed to lie on the simplex.

- **Loss**: MSE on CLR columns (averaged across K), applied group-wise so
  the whole composition is masked together during DAE corruption.

- **Metrics**: Aitchison distance (the natural compositional-data metric),
  RMSE on CLR space (comparable to continuous traits), and simplex MAE
  (mean absolute error on decoded proportions).

- **Benchmark**: `script/bench_multi_proportion.R` with 5-scenario
  signal sweep and K ∈ {3, 5, 8, 12} secondary sweep. Report rendered
  by `script/make_bench_multi_proportion_html.R`.

- **Tests**: `tests/testthat/test-multi-proportion.R` (5 test blocks,
  20 expectations) — encoding, validation, baseline, end-to-end `impute()`.

## Internal

- `evaluate_imputation()` output adds three columns: `aitchison`,
  `rmse_clr`, `simplex_mae`. Existing metrics columns are unchanged.

- `impute()` gains explicit `trait_types` and `multi_proportion_groups`
  arguments (previously `trait_types` was documented as passed through
  `...` but was not actually forwarded; this fixes that).

# pigauto 0.6.2

## Documentation

- **Tree uncertainty workflow clarified.** `multi_impute_trees()`,
  `trees300`, the README, and `vignette("getting-started")` now all
  describe the two-step workflow explicitly: step 1 is tree-aware
  imputation (pigauto's job — already done correctly); step 2 is
  tree-aware downstream analysis (the user's responsibility,
  following Nakagawa & de Villemereuil 2019, *Syst. Biol.* 68:632–641).
  Every doc now includes a complete `Map()`-over-`mi$tree_index` code
  example showing how to run the corresponding tree in the downstream
  model.
- **Compute-cost table added** to the tree-uncertainty docs: wall-clock
  budgets at `n` = 300 / 5,000 / 10,000 species with `T` = 10 / 50
  trees, plus guidance on reducing `T` for large trees and
  parallelising across machines.

## Internal

- No API changes; no numerical changes to any fitted model. Doc-only
  patch release.

# pigauto 0.6.1

## Bug fixes

- **MC dropout zero variance fixed** (`R/predict_pigauto.R`): when
  `n_imputations > 1`, the blend equation now uses a BM posterior draw
  `t_BM_draw ~ N(BM_mu, BM_se)` instead of the deterministic `t_MU`.
  When the calibrated gate is zero (BM dominates), each imputation now
  samples from the correct BM posterior rather than returning the same
  point estimate, giving non-zero between-imputation variance.

- **Conformal draw scale bug fixed** (`R/multi_impute.R`): for
  log-transformed traits (Mass, Wing length, etc.), draws are now made
  on the latent z-score scale and back-transformed, avoiding near-zero
  variance caused by dividing a log-scale SE by an original-scale value.
  The same fix applies to log1p-scaled counts and logit-scaled
  proportions.

## New features

- `draws_method = "conformal"` is now the default for `multi_impute()`.
  Draws sample from the split-conformal calibrated uncertainty
  distribution (better calibrated than MC dropout for downstream
  Rubin's-rules pooling).
- `draws_method = "mc_dropout"` remains available and is now correct
  (see bug fix above).

## Hex sticker

- New `man/figures/logo.png`: pink pig driving a car with a 7-tip
  ultrametric pectinate cladogram above the pig's head. Forest green /
  mint palette (v3) selected as the official sticker.

# pigauto 0.6.0

## Observation-level covariate refinement for multi-obs data

When datasets have multiple observations per species measured under
different conditions (e.g. CTmax at different acclimation temperatures),
the GNN now produces **covariate-conditional predictions within species**.
A new refinement MLP (`obs_refine`) re-injects user-supplied covariates
after species-level phylogenetic message passing, so that different
observations of the same species receive different predicted values
based on their covariate context. This is controlled by a residual
connection that preserves the phylogenetic signal from message passing.

- `R/model_residual_dae.R`: new `n_user_cov` parameter and `obs_refine`
  MLP (activated only when `species_col` + `covariates` are both supplied)
- `R/fit_pigauto.R`: `n_user_cov` stored in model config
- `R/predict_pigauto.R`: backward-compatible model reconstruction

## New bundled dataset: `ctmax_sim`

Simulated multi-observation-per-species CTmax data (1,464 observations
across 300 species, with acclimation temperature as an observation-level
covariate). 30% of species are entirely unobserved. Uses `tree300`.
See `data-raw/make_ctmax_sim.R` for the generation script.

## New benchmark: multi-observation imputation

`script/bench_multi_obs.R` simulates CTmax-like datasets under varying
phylogenetic signal and within-species acclimation response ratios,
comparing species_mean, pigauto (no covariates), and pigauto (with
covariates). The observation-level refinement gives measurable RMSE
improvement when within-species covariate effects are strong.

---

# pigauto 0.5.0

## Three sources of information

pigauto now explicitly combines three sources of information for
imputation: (1) the phylogenetic tree, (2) cross-trait correlations,
and (3) optional environmental covariates. The `covariates` argument
is accepted by `impute()`, `multi_impute()`, `multi_impute_trees()`,
and the lower-level pipeline functions. Covariates are threaded
through the GNN with the same gated safety that protects against
GNN degradation — they only contribute when they demonstrably
improve accuracy.

## Internal BM baseline (Rphylopars no longer required)

The Brownian-motion baseline for continuous, count, ordinal, and
proportion traits is now computed internally using conditional
multivariate normal imputation on the phylogenetic correlation matrix
`R = cov2cor(vcv(tree))` (Goolsby et al. 2017). This removes the
hard dependency on `Rphylopars`, which is now in `Suggests` only.
The internal implementation (`R/bm_internal.R`) uses univariate
per-column imputation with Cholesky decomposition and nugget
regularisation for near-singular submatrices. Validation against
`Rphylopars` shows r = 0.97–0.98 (expected given univariate vs
multivariate difference).

## Tree-uncertainty MI via `multi_impute_trees()`

New function `multi_impute_trees()` performs multiple imputation
across a posterior sample of phylogenetic trees, so that phylogenetic
uncertainty propagates into downstream standard errors via Rubin's
rules. Demonstrated with the bundled `trees300` dataset (10 BirdTree
posterior trees). Benchmark results show SE inflation of 1.1–2.1x
and fraction of missing information (FMI) rising from ~0.03
(single-tree) to 0.22–0.79 (multi-tree) depending on missingness
level.

## New bundled datasets

- `trees300`: 10 posterior trees from the BirdTree Hackett backbone,
  pruned to the same 300 tips as `tree300`
- `delhey5809`: plumage lightness and 4 environmental covariates
  (absolute latitude, mean temperature, mean precipitation,
  tree cover) for 5,809 bird species (Delhey 2019)
- `tree_delhey`: matching BirdTree phylogeny for `delhey5809`

## Benchmark suite

Thirteen benchmark reports are now accessible from the pkgdown site
under Articles > Benchmarks:

- Per-type benchmarks: continuous, binary, ordinal, count,
  categorical, proportion, zero-inflated counts
- Missingness mechanisms: MCAR vs MAR (trait/phylo) vs MNAR
- Tree uncertainty: single-tree vs multi-tree MI
- Environmental covariates: Delhey 5,809-species real data
  (high phylo signal, covariates give ~zero lift) and simulation
  sweep (4 lambda x beta scenarios; 8–15% RMSE reduction when
  env effects are strong, ~0% when phylo signal dominates)

Each report is a self-contained HTML page generated from
`script/make_bench_*_html.R` and backed by reproducible
`script/bench_*.R` drivers.

## Documentation

- All user-facing prose updated from "two sources" to "three
  sources" of information
- Rphylopars references replaced with "phylogenetic BM" or
  "internal conditional-MVN baseline" throughout
- `DOCS.md` updated with new benchmarks, datasets, and functions

---

# pigauto 0.4.0

## Breaking: TabPFN baseline removed

`fit_baseline_tabpfn()` and `setup_tabpfn()` are gone, along with the
`reticulate` suggested dependency, the `script/bench_tabpfn.*` files,
`script/run_tabpfn.py`, `tests/testthat/test-tabpfn.R`, and the
`dev/bench_tabpfn.html` benchmark page. The TabPFN wrapper was
originally added as a curiosity — a test of whether a general-purpose
tabular transformer could compete with pigauto as a drop-in baseline.
In practice it is a poor fit for what pigauto is trying to do:

1. **It is not phylogenetic.** TabPFN treats rows as exchangeable and
   cannot use the tree at all. The whole point of pigauto is that the
   tree carries information about trait evolution; feeding traits
   through a tree-blind model throws that away.
2. **Quadratic attention does not scale.** TabPFN's attention is
   O(n²) in the number of rows. At `n = 9,993` species (the full
   AVONET dataset) this is already near the edge of what a single
   GPU can handle, and most users of pigauto care precisely about
   scales at which a baseline must be cheap.
3. **The cross-language subprocess architecture was fragile.** R
   `torch` and Python `torch` share `libtorch` and cannot coexist in
   the same process, so TabPFN had to run via `system2()` on a
   separately-installed Python venv. Every time the user's Python,
   `tabimpute`, or `torch` versions drifted, the wrapper broke. This
   is not maintenance the project can afford.
4. **Keeping the benchmark page gave a misleading impression of
   relevance.** Users who hit the "pigauto vs TabPFN" page came away
   thinking TabPFN was a first-class supported baseline in pigauto.
   It never was.

Users who still want to benchmark pigauto against TabPFN can install
`tabimpute` separately and call it directly from Python. The removal
simplifies pigauto's installation, CI, and documentation.

## Docs cleanup (carried over from v0.3.3 work)

- `DOCS.md`: rewritten benchmark and developer-reports sections to
  reflect the v0.3.2+ benchmark set (scaling to 10,000 tips and the
  AVONET missingness sweep). Dropped the obsolete `benchmark_report`
  and `benchmark_final_report` links that v0.3.2 had already pruned
  from `_pkgdown.yml` but had not yet removed from `DOCS.md`.
- `CLAUDE.md`: removed the TabPFN Python-venv setup recipe and the
  "R torch ↔ Python torch cannot share a process" gotcha (no longer
  applicable). Version line bumped to 0.4.0.

# pigauto 0.3.3

## Correctness fix: Vphy must be a correlation matrix in `glmmTMB`'s `propto()`

The multiple-imputation example in `README.md`, the two HTML
tutorials (`pigauto_intro.html`, `pigauto_workflow_mixed.html`), and
`tests/testthat/test-multi-impute.R` all constructed the phylogenetic
random-effect structure with

```r
Vphy <- ape::vcv(tree)            # WRONG: covariance, diag = tree height
glmmTMB(... + propto(0 + species | dummy, Vphy), ...)
```

`propto()` estimates `sigma^2` freely, so passing the raw covariance
matrix (whose diagonal is the tree height) silently rescales the
phylogenetic variance component by tree height. The correct usage
passes a *correlation* matrix:

```r
Vphy <- cov2cor(ape::vcv(tree))   # diag = 1
```

All four sites are now fixed. Thanks to the user who flagged this.

## Honest framing: the GNN is a gated ensemble, not a residual model

The user-facing prose throughout the package historically described
pigauto as a *Residual Phylogenetic Denoising Autoencoder
(ResidualPhyloDAE) that learns a residual correction on top of a BM
baseline*. This wording was misleading for binary and categorical
traits in two ways:

1. The "baseline" for those types is **phylogenetic label propagation**
   (a similarity-weighted average of observed labels using a Gaussian
   kernel on the cophenetic distance) — not a generalised linear model
   that produces residuals on a link scale.
2. The GNN is trained **end-to-end** by minimising the type-appropriate
   loss (MSE for continuous/count/ordinal, BCE for binary,
   cross-entropy for categorical) on the **blended output**, not on
   `y - baseline`. So the GNN's output `delta` is a full per-cell
   prediction, not a statistical residual.

What pigauto actually is, in honest terms, is a **gated ensemble** of
two predictors:

- A phylogenetic **baseline** (Brownian motion via Rphylopars for
  continuous/count/ordinal traits; phylogenetic label propagation for
  binary/categorical traits).
- A **graph neural network correction** (an internal torch module that
  produces a full per-cell prediction `delta` from spectral node
  features and an attention-biased message passing over the
  phylogeny).

Prediction is the per-trait blend
`(1 - r_cal) * baseline + r_cal * delta`, where `r_cal` is calibrated
on a held-out validation split (with a split-validation cross-check
for discrete traits). When the baseline is already optimal the gate
closes and the GNN becomes a no-op.

The internal torch class name `ResidualPhyloDAE` is preserved because
its GNN layers use ResNet-style residual skip connections — that is a
correct, well-defined ML term. But user-facing prose no longer
describes the GNN as "learning a residual on top of the baseline".

Files updated: `README.md`, `DESCRIPTION`, `CLAUDE.md`, `_pkgdown.yml`,
`DOCS.md`, `vignettes/getting-started.Rmd`, `vignettes/mixed-types.Rmd`,
`R/fit_pigauto.R` (roxygen), `R/predict_pigauto.R` (roxygen),
`R/model_residual_dae.R` (clarifying header comment),
`R/simulate_traits.R`, `R/plot.R` (plot titles), and
`script/make_intro_html.R`.

## Tracking issue for v0.4.0 architectural work

The honest fix above is a documentation correction, not a model
change. A genuinely *residual* baseline for binary and ordinal traits
would require a liability-scale model (threshold BM, `phylolm::phyloglm`,
or `MCMCglmm` with `family = "threshold"`) — a design decision that
should be made deliberately rather than rolled into a same-day patch.
A v0.4.0 tracking issue lays out the trade-offs.

# pigauto 0.3.2

## New bundled dataset: `avonet_full` / `tree_full`

The full AVONET3 + BirdTree Stage2 Hackett MCC phylogeny, aligned to
9,993 species with 7 mixed-type traits (4 continuous morphometric +
2 categorical + 1 ordinal) and native missingness preserved.
Complements the 300-species `avonet300` / `tree300` bundled data:

- `avonet300` / `tree300` remain the right choice for quick examples,
  unit tests, and the getting-started vignette.
- `avonet_full` / `tree_full` are the right choice for scale
  benchmarks, realistic-missingness experiments, and anything that
  wants to exercise the v0.3.1 sparse Lanczos / cophenetic caching
  code paths on a real-world phylogeny.

The schemas are identical, so any code that runs on `avonet300` runs
on `avonet_full` with no modification. Combined tarball increment is
~328 KB (xz-9).

## New benchmark: AVONET missingness sweep (20 / 50 / 80%)

*Articles -> Benchmarks -> AVONET missingness sweep* is a new pkgdown
page comparing three methods on the full `avonet_full` dataset at
three MCAR missingness levels:

- **mean / mode** — column mean (continuous, ordinal) or class
  frequency (categorical). No phylogeny.
- **BM baseline** — Rphylopars Brownian-motion for continuous /
  ordinal, phylogenetic label propagation for discrete traits.
- **pigauto** — full pipeline (BM baseline + calibrated GNN delta
  + conformal intervals).

Per-trait RMSE, Pearson r, Spearman rho, and accuracy on held-out
test cells. Single seed per cell, trait-level masking so categorical
traits are held out as whole groups. Hyperparameters match the
`validate_avonet_full.R` scaling run for comparability.

## Docs cleanup

- Regenerated `pigauto_intro` and `pigauto_workflow_mixed` HTML
  tutorials against the v0.3.2 package state.
- Removed stale v0.3.0-era benchmark pages from the navbar
  (`benchmark_report.html`, `benchmark_final_report.html`,
  `test_report.html`) and from `pkgdown/assets/dev/`.
- Pruned 45 orphan files from `script/` (trial artefacts from
  pre-v0.3.1 benchmarking phases: `bench_v2/v3/v4.*`,
  `benchmark_avonet2000.*`, `benchmark_final.*`,
  `benchmark_simulation.*`, scaling duplicates, and Apr-6 demo
  artefacts).
- `README.md` citation now points at 0.3.2.
- `DESCRIPTION` updated to mention both bundled datasets.

# pigauto 0.3.1

## Scaling to 10,000 tips

`pigauto` now runs end-to-end on real 10,000-tip phylogenies. Before this
release the pipeline was tested and tuned around the 300-species bundled
AVONET data, and at 10,000 tips it wedged on a single dense eigendecomposition
of the graph Laplacian (hours of CPU time). Two targeted fixes land in this
release; no rewrite was needed.

### Fix A: sparse Lanczos eigensolver for `build_phylo_graph()`

`build_phylo_graph()` now uses `RSpectra::eigs_sym()` sparse Lanczos to
compute the `k + 1` smallest Laplacian eigenvectors when the tree has more
than 7,500 tips and **RSpectra** is installed. This replaces
`base::eigen(L, symmetric = TRUE)`, which was O(n^3) in time and discarded
>99% of its work (we only need the bottom ~32 eigenvectors).

- On the real 9,993-species AVONET + BirdTree Stage2 Hackett MCC phylogeny,
  the graph stage now takes under a minute, compared with roughly seven
  minutes of dense eigen() on the same hardware.
- The old dense path is kept as the default for smaller trees and as a
  safety fallback for when RSpectra is unavailable or Lanczos does not
  converge. The threshold is conservative (7,500) because highly symmetric
  ultrametric simulations (`ape::rcoal`) have degenerate spectra that need
  large Krylov subspaces to resolve.
- **RSpectra** is listed in `Suggests`. Installing it is strongly recommended
  for any tree larger than ~7,500 tips.

### Fix B: cache the cophenetic distance matrix across the pipeline

The cophenetic distance matrix `D = ape::cophenetic.phylo(tree)` was previously
computed four separate times in a single `impute()` call (three times inside
`build_phylo_graph()` and once inside `fit_baseline()`). Each call allocates a
dense n*n double matrix and runs an O(n^2) tree traversal -- wasted work at
scale. In 0.3.1:

- `build_phylo_graph()` computes `D` exactly once and returns it as part of
  the result list (`graph$D`).
- `fit_baseline()` gains an optional `graph` argument; when supplied, it
  reuses `graph$D` instead of recomputing the cophenetic matrix.
- `impute()` and `fit_pigauto()` pass the graph through automatically, so
  users get the speedup without any code changes.
- `impute()` also explicitly releases `graph$D` after `fit_baseline()`
  finishes, so the ~800 MB cophenetic matrix does not sit in R memory
  during training (a subtle side-effect that caused a large per-epoch
  regression at n = 10,000 until it was dropped).

### Scaling benchmark

See `pkgdown/assets/dev/scaling.html` (linked from the site as
*Articles -> Benchmarks -> Scaling to 10,000 tips*) for a side-by-side
comparison of v0.3.0 vs v0.3.1 at n in {300, 1k, 2k, 3k, 5k, 7.5k, 10k},
plus an end-to-end validation run on the full AVONET3 + BirdTree dataset
(9,993 species, 7 mixed-type traits). The pkgdown article also walks
through both a quick 300-species example with the bundled `avonet300`
data and a full 10k-species example using the BirdTree Stage2 Hackett
MCC phylogeny.

## Bug fixes

- `build_phylo_graph()`: the N > 2000 memory warning was out of date (the
  real blocker above that size was wall-clock time, not memory). The
  warning is now issued at N > 10000 and its text points at the true
  bottleneck.

# pigauto 0.3.0

## Major new features

### Attention-based message passing
The GNN now uses learned attention weights (with phylogenetic adjacency
as positional bias) instead of fixed adjacency matrix multiplication.
The model learns which phylogenetic neighbours are most informative for
each species, while the log-adjacency bias ensures initialisation
respects evolutionary distances. Controlled via `use_attention = TRUE`
(now the default).

### Validation-calibrated gates
After training, per-trait gate values are optimised on the validation set
via grid search. This prevents the GNN from adding noise when the
Brownian motion baseline is already optimal. On the AVONET benchmark,
this eliminates the RMSE regression seen in earlier versions.

### Conformal prediction intervals
Distribution-free 95% prediction intervals computed via split conformal
prediction on the validation residuals. Coverage is guaranteed under
exchangeability. Available for continuous, count, and ordinal traits via
`pred$conformal_lower` and `pred$conformal_upper`.

### Phylogenetic label propagation
Binary and categorical trait baselines now use phylogenetic similarity-
weighted label propagation instead of flat population-level frequencies.
Each species gets a personalised baseline reflecting its phylogenetic
neighbourhood.

### Adaptive spectral encoding
`k_eigen` now defaults to `"auto"`, scaling with tree size:
`min(max(ceil(n/20), 4), 32)`. This gives 4 features for small trees,
15 for 300 tips, and 32 for 640+ tips.

### Auto-generated HTML reports
`pigauto_report()` produces self-contained HTML reports with interactive
Chart.js visualisations: per-trait metrics, calibrated gate values,
training history, and conformal coverage.

### Evaluation framework
- `evaluate()`: per-trait metrics on held-out data
- `compare_methods()`: automated BM vs GNN comparison
- `cross_validate()`: proper k-fold cross-validation
- `summary.pigauto_fit()`: formatted console summary

### Plotting
- `plot.pigauto_fit()`: training history, gate values, conformal scores
- `plot.pigauto_pred()`: scatter plots, conformal intervals, probability
  distributions
- `plot_comparison()`: forest-plot style method comparison

### Simulation benchmark
`simulate_benchmark()` runs a complete simulation study with configurable
scenarios (BM, OU, regime shift, non-linear, mixed types), returning tidy
results with print/summary/plot methods. Useful for methods papers and
for verifying pigauto's behaviour on data with known properties.

## Minor improvements

- Categorical gate calibration now works at the trait level (all K one-hot
  columns share a single gate) using cross-entropy loss, not per-column MSE
- `build_phylo_graph()` default `k_eigen` changed from `8L` to `"auto"`
- `fit_pigauto()` default `k_eigen` changed from `8L` to `"auto"`
- `all_grads_finite()` now handles non-standard tensor types gracefully
- `impute()` uses adaptive k_eigen by default
- Package passes `R CMD check` with 0 errors and 0 notes
- Bumped version to 0.3.0

# pigauto 0.2.0

- Added support for multiple observations per species (`species_col`)
- Added mixed trait type support (continuous, binary, categorical,
  ordinal, count)
- Added AVONET 300 bundled dataset with mixed types
- Extended simulation benchmark (16 scenarios via BACE)

# pigauto 0.1.0

- Initial release
- ResidualPhyloDAE architecture
- Brownian motion baseline via Rphylopars
- Per-column gated blending
- MC dropout for multiple imputation
