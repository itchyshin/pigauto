# BACE Comparison Framework for pigauto

## Critical Differences in Evaluation Approach

### What BACE Emphasizes
1. **Convergence of the imputation process** (not parameter estimation)
   - Tracks stability of imputations across multiple chained iterations
   - Uses statistical tests on summary statistics (ACF, trend, % change)
   - Auto-retries if not converged

2. **Posterior uncertainty quantification**
   - Multiple imputations (usually n_final=10)
   - Pooling via simple concatenation of MCMC samples
   - Final posteriors reflect both model + imputation uncertainty

3. **Bayesian inference compatibility**
   - Output is MCMCglmm objects with full posterior distributions
   - Allows credible intervals, posterior predictions, etc.
   - Not just point estimates

4. **Phylogenetic integration**
   - Phylogenetic covariance matrix in random effects
   - Optional decomposition of phylogenetic + non-phylo effects
   - Critical for comparative datasets

### What BACE De-Emphasizes (In Documentation)
- **Quantitative imputation accuracy metrics** (not in main vignette)
- **Runtime/scalability comparisons** (AVONET script only, n=2000)
- **Comparison with non-Bayesian methods** (only vs Rphylopars)
- **Uncertainty calibration** (posterior SD vs actual error)

---

## Fair Comparison Metrics for pigauto vs BACE

### 1. Imputation Accuracy (Quantitative)
Required for any comparison paper:
```r
# On held-out masked data
RMSE <- sqrt(mean((pred - true)^2))
MAE <- mean(abs(pred - true))
Correlation <- cor(pred, true)

# Can calculate per-variable and overall
# Should vary by:
#  - Variable type (continuous vs categorical)
#  - Missingness mechanism (MCAR, MAR, MNAR)
#  - Proportion missing (5%, 10%, 25%, 50%)
```

### 2. Uncertainty Quantification
```r
# BACE returns full posteriors; pigauto should provide comparable uncertainty
# Coverage: proportion of times true value falls in 95% CI
# Width: average credible interval width
# Calibration: {actual error} vs {posterior SD}

# For each imputed value:
true_val <- truth[i]
posterior_samples <- pooled$models$trait$Sol[, "intercept"]
ci_lower <- quantile(posterior_samples, 0.025)
ci_upper <- quantile(posterior_samples, 0.975)
coverage <- mean(true_val >= ci_lower & true_val <= ci_upper)
```

### 3. Convergence Properties
```r
# BACE's advantage: explicit convergence checking
# Compare:
#  - Number of iterations to stability
#  - Sensitivity to MCMC settings
#  - Failure rate under challenging conditions

# Key scenarios:
#  - High missingness (>50%)
#  - Mixed variable types
#  - Strong phylogenetic signal
#  - Weak phylogenetic signal
```

### 4. Mixed Variable Type Handling
BACE's core strength. Test:
```r
# Dataset with:
#  - Continuous (gaussian) ~ 3 variables
#  - Binary ~ 2 variables
#  - Ordinal/threshold (4 levels) ~ 2 variables
#  - Categorical (5 levels) ~ 2 variables

# Metrics per type:
#  - Continuous: RMSE, correlation
#  - Binary: Accuracy, AUC (if using predictions for classification)
#  - Ordinal: Ordered accuracy, Spearman correlation
#  - Categorical: Accuracy on majority class, multinomial AUC
```

### 5. Computational Efficiency
```r
# BACE's potential weakness
# Compare on same hardware:

library(tictoc)

# Dataset sizes:
#  - Small (n=50 species, p=10 variables)
#  - Medium (n=200 species, p=20 variables)
#  - Large (n=2000 species, p=40 variables) [AVONET scale]

# Metrics:
#  - Runtime (seconds)
#  - Memory usage (MB)
#  - Scalability (runtime vs n_species)
```

### 6. Phylogenetic Signal Effect
```r
# BACE integrates phylogeny; pigauto may not
# Test datasets with varying phylogenetic signals:

phylo_signal_levels <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# For each level:
#  - Generate data with sim_bace(phylo_signal = level)
#  - Mask values
#  - Run both methods
#  - Compare accuracy

# Expected: BACE gains advantage at moderate-high signals
#           Both similar when signal is weak (random tree)
```

---

## Test Data Recommendations

### 1. Synthetic (Controlled)
Use BACE's `sim_bace()` to ensure comparability:
```r
# Example configuration
sim <- sim_bace(
  response_type = "gaussian",
  predictor_types = c("gaussian", "gaussian", "binary", "threshold4"),
  phylo_signal = 0.6,
  missingness = c(0.2, 0.2, 0.3, 0.2, 0.2),  # Vary by variable
  n_cases = 200,
  n_species = 50,
  seed = 123
)

# Store both:
complete_data <- sim$complete_data      # Ground truth
masked_data <- sim$data                 # Has NAs
tree <- sim$tree
truth_values <- # Extract the masked locations
```

### 2. Real Data (Reproducibility)
Use AVONET 2000 species from BACE repo:
```
/dev/testing_data/
  - avonet_2000_masked.csv       (input)
  - avonet_2000_truth.csv        (gold standard)
  - Hackett_tree_2000.tre        (phylogeny)
```

---

## Scenario-Specific Comparisons

### Scenario 1: Well-Behaved Data (Baseline)
- No strong missingness patterns
- Moderate sample size (n=200)
- All continuous or well-behaved mixed types
- Moderate phylogenetic signal
- **Expected winner**: Both similar; pigauto may be faster

### Scenario 2: High-Dimensional Mixed Types
- 40+ variables of mixed types
- 10-20% missing per variable
- Strong correlations between variables
- **Expected winner**: BACE (handles all types natively)

### Scenario 3: Extreme Missingness
- 50%+ missing in some variables
- Requires chaining iterations
- **Expected winner**: BACE (convergence designed for this)

### Scenario 4: Large-Scale Phylogenetic
- n=2000+ species
- Phylogenetic relationships critical
- Memory/speed important
- **Expected winner**: pigauto (if memory-efficient)

---

## Output Comparison

### BACE Outputs
```r
result <- bace(...)

# 1. Imputed datasets
result$imputed_datasets[[1]]      # Data frame with imputed values
result$imputed_datasets[[2]]
...                               # n_final datasets

# 2. Posterior distributions
result$pooled_models$models$var_name$Sol[, "fixed_effect"]  # All samples
quantile(posterior_samples, c(0.025, 0.5, 0.975))           # CI

# 3. Convergence info
result$convergence$converged      # TRUE/FALSE
result$convergence$summary_stats  # Statistics across runs
```

### pigauto Should Provide
```r
# Comparable outputs:
pigauto_result <- pigauto_impute(...)

# 1. Imputed datasets (matching BACE)
pigauto_result$imputed_datasets[[1]]

# 2. Uncertainty representation
# Options:
#   a) Multiple imputations (like BACE) + posterior aggregation
#   b) Bayesian posterior samples (like BACE)
#   c) Frequentist confidence intervals
#   d) Bootstrap samples

# 3. Convergence/validation
pigauto_result$convergence_info  # Comparable to BACE
```

---

## Key Metrics Table for Comparison

| Metric | How to Measure | BACE Strength | pigauto Strength |
|--------|---|---|---|
| **Continuous accuracy** | RMSE on held-out | Similar | Similar |
| **Categorical accuracy** | Classification error | Strong (native) | ? |
| **Ordinal accuracy** | Ordered accuracy | Strong (native) | ? |
| **Convergence speed** | Iterations to stable | ? | Potentially faster |
| **Computational time** | Runtime(seconds) | ? (slow) | Potentially better |
| **Memory usage** | RAM required | Uses sampling to manage | ? |
| **Phylogenetic benefit** | Gain from using tree | Strong (built-in) | Depends on design |
| **Uncertainty calibration** | Coverage of CIs | Good (Bayesian) | ? |
| **Mixed type handling** | Performance on each type | Excellent | To be tested |
| **Scalability** | Time vs n, p | Limited (Bayesian MCMC) | ? |

---

## Specific Tests from BACE's Own Comparison

BACE vs Rphylopars script shows:
```r
# Metrics computed per trait:
metrics <- preds %>%
  group_by(trait) %>%
  summarise(
    n = n(),
    cor_bace  = cor(true_value, bace_pred),
    rmse_bace = sqrt(mean((bace_pred - true_value)^2)),
    mae_bace  = mean(abs(bace_pred - true_value)),
    
    cor_rphy  = cor(true_value, rphy_pred),
    rmse_rphy = sqrt(mean((rphy_pred - true_value)^2)),
    mae_rphy  = mean(abs(rphy_pred - true_value))
  )
```

**This is a good template for pigauto comparison**

---

## Checklist for Fair Comparison

- [ ] Same synthetic datasets generated with `sim_bace()` (reproducible seed)
- [ ] Same held-out validation locations (mark before running either method)
- [ ] Same evaluation metrics (RMSE, MAE, correlation per variable type)
- [ ] Same hardware/computational setup
- [ ] Multiple random seeds (5-10 replications)
- [ ] Varied scenarios (small/large n, high/low missingness, varied phylo signal)
- [ ] Results table showing metric ± SD across replicates
- [ ] Statistical testing (e.g., paired t-tests if appropriate)
- [ ] Runtime/memory comparison on realistic scale
- [ ] Reproducible code in supplementary materials

