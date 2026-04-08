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
