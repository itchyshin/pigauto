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
