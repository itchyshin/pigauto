# BACE Package - Quick Reference Guide

## Key Functions & Signatures

### Main Entry Point: `bace()`
```r
bace(
  fixformula,           # Formula or list of formulas: "y ~ x1 + x2"
  ran_phylo_form,       # Random effects: "~ 1 | species"
  phylo,               # ape::phylo object
  data,                # data.frame with NAs
  runs = 10,           # Initial convergence checking runs
  nitt = 6000,         # MCMC iterations
  thin = 5,            # Thinning rate
  burnin = 1000,       # Burn-in iterations
  n_final = 10,        # Final pooling runs (after convergence)
  species = FALSE,     # Decompose phylogenetic + non-phylo effects
  verbose = TRUE,
  plot = FALSE,
  max_attempts = 3,    # Retry if not converged
  skip_conv = FALSE,   # Skip convergence checking
  sample_size = NULL   # Optional: subsample posteriors
)
```
**Returns**: `bace_complete` object with `$pooled_models`, `$imputed_datasets`, `$convergence`

---

### Core Imputation: `bace_imp()`
```r
bace_imp(
  fixformula,
  ran_phylo_form,
  phylo,
  data,
  nitt = 6000,
  thin = 5,
  burnin = 1000,
  runs = 10,           # Number of chained iterations
  species = FALSE,
  verbose = TRUE,
  ...
)
```
**Returns**: `bace` object with `$data` (list of imputed datasets), `$models`

---

### Convergence Assessment: `assess_convergence()`
```r
assess_convergence(
  bace_object,                                    # Output from bace_imp()
  method = "summary",                             # or "all", "energy", "wasserstein"
  variables = NULL,                               # Specific variables to check
  use_all_data = FALSE,                          # Use only imputed values
  alpha = 0.05,                                   # Significance level
  min_iterations = 3,
  lag = 1,                                        # ACF lag for tests
  pct_change_threshold = 0.05                     # 5% change threshold
)
```
**Returns**: `bace_convergence` object with `$converged` (logical), `$summary_stats`, diagnostic details

---

### Posterior Pooling: `pool_posteriors()`
```r
pool_posteriors(
  bace_final_object,   # Output from bace_final_imp()
  variable = NULL,     # Specific variable or all
  sample_size = NULL   # Optional: subsample posteriors (e.g., 1000)
)
```
**Returns**: `bace_pooled` object with `$models` (MCMCglmm objects), `$n_imputations`

---

### Data Simulation: `sim_bace()`
```r
sim_bace(
  response_type = "gaussian",                     # or "poisson", "binary"
  predictor_types = c("gaussian", "gaussian"),    # Vector of predictor types
  var_names = NULL,
  phylo_signal = NULL,                            # Vector of phylo signal (0-1)
  n_cases = 200,
  n_species = 75,
  missingness = NULL,                             # Vector of missing proportions
  birth = 0.8,                                    # Tree generation params
  death = 0.4,
  ...
)
```
**Returns**: List with `$data`, `$tree`, `$params`, `$complete_data`

---

## File Locations & Contents

| Function | File | Key Details |
|----------|------|-------------|
| `bace()` | `/R/bace.R` (270 lines) | Wrapper orchestrating entire workflow |
| `bace_imp()` | `/R/bace_imp.R` (400+ lines) | Core chained equations imputation |
| `assess_convergence()` | `/R/convergence_diagnostics.R` (800+ lines) | Multiple convergence methods |
| `pool_posteriors()` | `/R/pool_posteriors.R` (360 lines) | Posterior pooling & uncertainty |
| `sim_bace()` | `/R/simulate_simBACE.R` (400+ lines) | Data simulation with mixed types |
| Variable type detection | `/R/prep_functions.R` | `.summarise_var_types()` |
| Formula building | `/R/build_functions.R` | `.build_formula_string()` |
| Model fitting | `/R/model_functions.R` | `.fit_mcmc_model()` |

---

## Variable Type Codes

BACE automatically detects:
- `gaussian`: numeric (continuous)
- `poisson`: numeric (count-like)
- `binary`: factor with exactly 2 levels
- `threshold`: ordered factor
- `categorical`: unordered factor with >2 levels

Access detected types via: `bace_object$types`

---

## Output Structure Examples

### From `bace()`
```r
result$pooled_models
  $models$y           # MCMCglmm object (fully pooled)
  $models$x1
  $n_imputations      # Number pooled (default 10)
  
result$imputed_datasets  # List of n_final datasets
  [[1]]$y, [[1]]$x1, ...
  [[2]]$y, [[2]]$x1, ...
  
result$convergence     # Details from assess_convergence()
  $converged          # TRUE/FALSE
  $method_results
  $summary_stats      # Statistics across iterations
  
result$converged       # Boolean
result$n_attempts      # Number of retry attempts
```

### From `pool_posteriors()`
```r
pooled$models$y$Sol               # Posterior samples for fixed effects
pooled$models$y$VCV               # Posterior samples for random variance
pooled$models$y$DIC               # Model DIC
pooled$models$y$BACE_pooling$n_imputations    # Metadata
```

---

## Convergence Criteria (Summary Method)

Variable considered converged if:
- **ACF**: Lag-1 autocorrelation < 0.3
- **% Change**: Mean absolute change in last half of runs < 5%
- **Trend**: No significant linear trend (p > 0.05)
- **Overall**: Majority of variables pass majority of tests

---

## Test Files & Data

| File | Type | Purpose |
|------|------|---------|
| `tests/testthat/test-bace_imp.R` | Unit tests | Formula/error handling |
| `tests/testthat/test-convergence_diagnostics.R` | Unit tests | Convergence methods |
| `dev/testing_data/avonet_2000_masked.csv` | Real data | 2000 bird species traits |
| `dev/testing_data/avonet_2000_truth.csv` | Validation | Held-out true values |
| `dev/testing_data/run_bace_vs_rphylopars.R` | Comparison script | Benchmark vs Rphylopars |

---

## Typical Workflow

```r
# 1. Create or load phylogenetic tree
library(ape)
tree <- read.tree("my_tree.nwk")
tree <- compute.brlen(tree, method = "Grafen")

# 2. Prepare data (ENSURE proper variable types!)
data <- read.csv("my_data.csv")
str(data)  # Check types before proceeding

# 3. Run full BACE analysis
result <- bace(
  fixformula = "trait1 ~ trait2 + trait3",
  ran_phylo_form = "~ 1 | species",
  phylo = tree,
  data = data,
  runs = 15,
  n_final = 10,
  verbose = TRUE
)

# 4. Check results
summary(result$pooled_models$models$trait1)
result$pooled_models$models$trait1$Sol  # Posterior samples

# 5. Access imputed datasets for downstream analysis
imputed_data_1 <- result$imputed_datasets[[1]]
```

---

## Key Implementation Details

1. **Chained equations**: Variables fit sequentially, using predictions from previous variables
2. **Phylogenetic covariance**: Integrated via MCMCglmm `ginverse` parameter
3. **Mixed types**: Each variable gets appropriate GLM (Gaussian, Poisson, Binomial, Threshold, Multinomial)
4. **Convergence**: Checked via summary statistics + optional energy/Wasserstein distances
5. **Pooling**: Simple concatenation of MCMC samples (properly weighted by imputation)
6. **Memory**: Optional posterior sampling reduces final object size

---

## Important Notes

1. **Data preparation is critical**: Variables must be properly classified (numeric vs factor) before imputation
2. **Convergence may require tuning**: MCMC parameters (nitt, thin, burnin) may need adjustment
3. **Phylogenetic signal matters**: Results depend on phylogenetic relationships
4. **Replicated observations**: Needed for `species=TRUE` decomposition (>30% species)
5. **Output is MCMCglmm**: Can use all standard MCMCglmm methods on pooled results

