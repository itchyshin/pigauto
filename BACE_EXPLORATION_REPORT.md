# BACE Package Exploration Summary

## Overview
BACE (Bayesian Augmentation by Chained Equations) is an R package for imputing missing data in comparative phylogenetic datasets using Bayesian methods with chained equations. It's designed to handle **mixed variable types** (continuous, binary, categorical, ordinal) while accounting for phylogenetic relationships.

---

## 1. Bundled Datasets
- **No built-in bundled datasets** in `data/` or `data-raw/`
- Package includes **real-world test data** in `/dev/testing_data/`:
  - `AVONET.csv` (2MB+ - 2000 bird species with trait measurements)
  - `avonet_2000_masked.csv` (masked version with artificial missingness)
  - `avonet_2000_truth.csv` (held-out true values for evaluation)
  - `Hackett_tree_2000.tre` (pruned phylogenetic tree matching the AVONET data)
- Provides `sim_bace()` function for generating synthetic datasets with controlled properties

---

## 2. Main User-Facing Functions

### Primary: `bace()` - Complete Workflow
**Location**: `/R/bace.R`
- **Purpose**: High-level wrapper orchestrating the full analysis pipeline
- **Input Parameters**:
  - `fixformula`: Fixed effects formula (string or list of formulas)
  - `ran_phylo_form`: Random effects + phylogenetic structure
  - `phylo`: phylo-class tree (from ape package)
  - `data`: DataFrame with missing values
  - `runs`, `nitt`, `thin`, `burnin`: MCMC controls
  - `n_final`: Number of final imputation runs (default=10)
  - `species`: Decompose phylogenetic + non-phylogenetic effects
  - `verbose`, `plot`: Output controls
  - `sample_size`: Optional posterior sampling to reduce memory
  
- **Output**: `bace_complete` object containing:
  - `pooled_models`: MCMCglmm objects with posterior pooling
  - `imputed_datasets`: List of final imputed datasets
  - `convergence`: Convergence diagnostics
  - `converged`: Boolean convergence status
  - `n_attempts`: Number of attempts to converge

- **Workflow**: 
  1. Initial imputation runs (`bace_imp()`) for convergence checking
  2. Convergence assessment (`assess_convergence()`)
  3. Auto-retry with more runs if not converged (up to `max_attempts`)
  4. Final imputation runs (10 by default)
  5. Posterior pooling across imputations

### Secondary: `bace_imp()` - Core Imputation
**Location**: `/R/bace_imp.R`
- Single-iteration imputation using chained equations
- Can be called directly for more control over the process
- Returns `bace` object with imputed datasets and models

### Utility: `sim_bace()` - Data Simulator
**Location**: `/R/simulate_simBACE.R`
- Generates synthetic phylogenetic comparative datasets
- **Supports multiple variable types**:
  - `gaussian`: Continuous normal
  - `binary`: Binary (0/1)
  - `poisson`: Count data
  - `thresholdK`: Ordinal (K levels)
  - `multinomialK`: Categorical (K levels)
- Parameters:
  - `response_type`: Type of response variable
  - `predictor_types`: Vector of predictor types
  - `phylo_signal`: Strength of phylogenetic signal (0-1)
  - `missingness`: Proportion of missing values per variable
  - `n_cases`, `n_species`: Data dimensions
  - `rr`: Enable random slopes
- Returns: Complete dataset + tree + simulation parameters

---

## 3. Variable Type Handling

BACE handles **mixed types** by automatically detecting and fitting appropriate models:

| Variable Type | BACE Name | Detection | Imputation Model |
|---------------|-----------|-----------|------------------|
| Continuous | `gaussian` | numeric | Linear Gaussian |
| Count | `poisson` | numeric (integer-like) | Poisson GLM |
| Binary | `binary` | factor with 2 levels | Binomial GLM |
| Ordinal | `threshold` | ordered factor | Threshold (ordinal) |
| Categorical | `categorical` | unordered factor | Multinomial logit |

**Key Note from Vignette**: Users MUST ensure proper variable classification before imputation:
```r
# Check variable types
str(data)  # Factors for categorical, numeric for continuous, etc.
```

---

## 4. Imputation Process

### Chained Equations Approach
1. **Iteration structure**: Multiple "runs" of the chaining process
2. **Within each run**:
   - Models fit for each variable with missing data
   - Initial predictor missingness filled with:
     - Mean (continuous variables)
     - Random sample from observed values (categorical variables)
   - Models predict missing values for response
   - Next iteration uses these predictions as predictors
   - Cycle continues until convergence
3. **Final step**: Pool posteriors across multiple imputation runs

### Phylogenetic Integration
- **Random effect formula**: `ran_phylo_form = "~ 1 | species"`
- Phylogenetic variance captured via phylogenetic covariance matrix (`ginverse`)
- Optional species effect decomposition (`species=TRUE`):
  - Both phylogenetic + non-phylogenetic species random effects
  - Requires replicated observations per species (>30% recommended)

---

## 5. Convergence Assessment

**Location**: `/R/convergence_diagnostics.R`

BACE implements **multiple convergence diagnostics**:

### Summary Statistics Method (Primary)
- Tracks summary statistics across iteration runs
- Applies four statistical tests to each variable:
  1. **ACF (Autocorrelation)**: abs(lag-1 correlation) < 0.3 → converged
  2. **Percentage Change**: Mean absolute % change in last half of runs < 5% → converged
  3. **Trend Test**: No significant linear trend (p > 0.05) → converged
  4. **Geweke Test**: First 10% vs last 50% not significantly different → converged (currently disabled)
- Overall: Majority of variables must pass majority of tests

### Energy Distance Method
- Calculates energy distance between consecutive imputed datasets
- Checks for stabilization in later iterations
- Tests: No increasing trend + Low CV in last half

### Wasserstein Distance Method
- Per-variable 1D Wasserstein distances (continuous variables only)
- Tests: Low CV + No increasing trend in last half
- >70% of variables must stabilize

### Overall Decision
- **All methods**: Majority vote (≥2 of 3 methods pass)
- **Summary only** (default): All variables must show convergence
- User can select specific methods

---

## 6. Posterior Pooling

**Location**: `/R/pool_posteriors.R`

### Pooling Strategy
- **Concatenates MCMC samples** from all imputation runs
- Naturally accounts for both within and between-imputation variance
- Sample pooling: Combines `nrow(Sol)` per imputation × n_imputations
- Optional posterior sampling to reduce memory (e.g., draw 1000 from each imputation)

### Output Format
- Returns `bace_pooled` objects containing:
  - `models`: Named list of pooled MCMCglmm objects
  - `n_imputations`: Number of pooled imputations
  - `variables`: Names of imputed variables
- Each pooled model inherits MCMCglmm class → works with standard methods:
  - `summary()`, `print()`, `plot()`
  - Extract posteriors via `model$Sol[, "coefficient"]`

### Metadata Tracking
- Each pooled model includes `BACE_pooling` list:
  - `n_imputations`: Number of datasets pooled
  - `n_samples_per_imputation`: Samples per imputation (after optional sampling)
  - `original_samples_per_imputation`: Before sampling
  - `total_samples`: Final posterior size

---

## 7. Evaluation Approach

### Vignette Examples
BACE documentation focuses on **qualitative evaluation**:
- Summary statistics of pooled models
- Effective sample size (ESS) inspection
- Trace and density plots for mixing diagnostics
- Posterior distributions of fixed effects

### Comparison Script: `run_bace_vs_rphylopars.R`
Real comparison against **Rphylopars** on AVONET (2000 species):

**Evaluation Metrics Used**:
1. **Correlation**: `cor(true_value, imputed_prediction)`
2. **RMSE**: `sqrt(mean((prediction - true)^2))`
3. **MAE**: `mean(abs(prediction - true))`
4. **Method agreement**: `cor(bace_prediction, rphylopars_prediction)`

**Held-Out Validation Design**:
- Original traits masked at random locations
- Ground truth values stored separately
- Predictions evaluated only on masked cells
- Metrics calculated per trait

**Data Preprocessing**:
- Log transformation for skewed traits (user-specified)
- Safe transforms handling zeros/negatives via shifts
- Back-transformation for final evaluation on original scale

---

## 8. Testing & Validation

**Location**: `/tests/testthat/`

### Test Suite Coverage
- `test-bace_imp.R`: Formula validation, error handling
- `test-bace_wrapper.R`: Wrapper function tests
- `test-posterior_sampling.R`: Pooling logic
- `test-convergence_diagnostics.R`: Convergence methods
- `test-species_effects.R`: Phylogenetic + species decomposition
- `test-mcmc_param_lists.R`: Model-specific MCMC parameters

### Test Data
- Minimal synthetic trees (10 species)
- Small datasets (10 observations)
- Multiple variable types included
- Focus on input validation and error messages

---

## 9. Key Design Features

### Strengths
1. **Mixed variable type support**: Handles continuous, count, binary, ordinal, categorical
2. **Phylogenetic integration**: Incorporates evolutionary relationships
3. **Multiple imputation**: Proper uncertainty quantification
4. **Convergence checking**: Automated detection + retry mechanism
5. **Flexible formulas**: List-based formulas per variable
6. **Posterior pooling**: Accounts for imputation uncertainty in final inference
7. **Memory optimization**: Optional posterior sampling for large datasets

### Important Assumptions
1. **Data must be completely classified** (factors vs numeric) before imputation
2. **Phylogenetic relationships matter** for predictive quality
3. **Replicated observations per species** recommended for effect decomposition
4. **Convergence may require iterative tuning** of nitt/thin/burnin

---

## 10. Relevant Files for Comparison with pigauto

| Purpose | File Path | Notes |
|---------|-----------|-------|
| **Main interface** | `R/bace.R` | High-level workflow orchestration |
| **Imputation logic** | `R/bace_imp.R` | Core chained equations algorithm |
| **Convergence** | `R/convergence_diagnostics.R` | Multiple diagnostic methods |
| **Posterior pooling** | `R/pool_posteriors.R` | Uncertainty quantification |
| **Variable type handling** | `R/prep_functions.R` | Type detection + preprocessing |
| **Vignette example** | `vignettes/bace.qmd` | Practical workflow + evaluation |
| **Real comparison** | `dev/testing_data/run_bace_vs_rphylopars.R` | Held-out validation design |
| **Data simulator** | `R/simulate_simBACE.R` | Synthetic data with mixed types |
| **Test data** | `dev/testing_data/avonet_*.csv` | Real phylogenetic trait data |

---

## Summary: BACE Evaluation Strategy

BACE's evaluation emphasizes:
1. **Within-imputation diagnostics**: ESS, mixing (via plots)
2. **Between-imputation convergence**: Stability across runs (statistical tests)
3. **Posterior uncertainty**: Proper pooling across imputations
4. **Held-out validation** (in comparison scripts): Correlation, RMSE, MAE on masked data

For fair pigauto comparison, key metrics should include:
- **Imputation accuracy**: RMSE, MAE, correlation on held-out values
- **Uncertainty quantification**: Posterior SD vs actual error
- **Convergence properties**: Speed to stability across runs
- **Mixed variable handling**: Performance across all data types
- **Scalability**: Runtime + memory on large datasets (BACE uses sampling strategy)
