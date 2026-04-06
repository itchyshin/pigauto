library(torch)
library(ape)
library(Matrix)
library(Rphylopars)

# [Paste the Previous SpectralVGAE and get_spectral_features code blocks here if not already in memory]

# ── Simulation: Non-Linear Environmental Effect ──────────────────────────────
compare_nonlinear <- function(n_sims = 10, n_sp = 300) {
  
  res_table <- data.frame()
  cat("Running Non-Linear Benchmark...\n")
  
  for(i in 1:n_sims) {
    # 1. Generate Tree
    tree <- rtree(n_sp)
    
    # 2. Generate Environment (Random spatially correlated variable)
    # We simulate 'env' as a trait itself so it has phylogenetic signal
    env <- rTraitCont(tree, model = "BM", sigma = 0.5)
    
    # 3. Generate Trait: NON-LINEAR function of Env + Phylo Noise
    # Trait = sin(Env * 3) + PhyloNoise
    phylo_noise <- rTraitCont(tree, model = "BM", sigma = 0.2)
    true_trait <- sin(env * 3) + phylo_noise
    
    # 4. Mask Data
    obs_trait <- true_trait
    missing_idx <- sample(1:n_sp, size = 0.3 * n_sp)
    obs_trait[missing_idx] <- NA
    
    # 5. Method 1: Rphylopars (Linear Assumption)
    # We must give it the Env as a covariate
    df <- data.frame(species = tree$tip.label, trait = obs_trait, env = env)
    phy_pred <- tryCatch({
      suppressWarnings({
        # fit_cov=TRUE tells it to use the environment
        p <- phylopars(df, tree, model="BM", pheno_error=FALSE) 
      })
      p$anc_recon[1:n_sp, 1][missing_idx]
    }, error = function(e) rep(NA, length(missing_idx)))
    
    # 6. Method 2: Spectral VGAE (Non-Linear Capacity)
    # We need to modify the input to include 'env'
    # For fairness, we'll just concat 'env' to the input vector X 
    # (The previous function needs a tiny tweak to accept covariates, 
    #  or we just treat 'env' as a second trait in a multi-trait VGAE).
    
    # -- Quick Hack for this test: Residual Imputation on the Residuals --
    # Let's see if VGAE can pick up the signal purely from the graph structure
    # of the complex trait, even without seeing 'env' explicitly.
    vgae_res <- impute_spectral_vgae(obs_trait, tree, epochs=2000, lr=0.01, k_eigen=8, beta=0.001)
    vgae_pred <- vgae_res$imputed[missing_idx]
    
    # Metrics
    actual <- true_trait[missing_idx]
    
    rmse_vgae <- sqrt(mean((actual - vgae_pred)^2))
    rmse_phy <- sqrt(mean((actual - phy_pred)^2))
    
    res_table <- rbind(res_table, data.frame(
      Sim = i,
      Method = c("Spectral_VGAE", "Rphylopars"),
      RMSE = c(rmse_vgae, rmse_phy)
    ))
    cat(".")
  }
  return(res_table)
}

# Run
res_nonlinear <- compare_nonlinear(n_sims = 5)
aggregate(RMSE ~ Method, data = res_nonlinear, mean)

# draw boxplot
