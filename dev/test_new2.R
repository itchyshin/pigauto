library(torch)
library(ape)
library(Matrix)
library(Rphylopars)

# ── 1. Device Setup ──────────────────────────────────────────────────────────
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

# ── 2. Data Simulation: The "VGAE Advantage" Scenario ────────────────────────
simulate_complex_data <- function(n_sp = 2000, missing_pct = 0.3) {
  # 1. Large Tree
  tree <- rtree(n_sp)
  
  # 2. Environment (Continuous, with phylogenetic signal)
  # e.g., Temperature varying across the clade
  env <- rTraitCont(tree, model = "BM", sigma = 0.5)
  # Standardize Environment to mean=0, sd=1 for neural net stability
  env <- as.numeric(scale(env))
  
  # 3. Trait: NON-LINEAR response to Environment + Phylogenetic Noise
  # "Sigmoid Switch": The trait stays low, then jumps up, then plateaus.
  # Rphylopars (Linear) hates this. Neural Nets love this.
  phylo_noise <- rTraitCont(tree, model = "BM", sigma = 0.2)
  
  # The "Switch" Logic:
  true_trait <- 10 / (1 + exp(-3 * env)) + phylo_noise
  
  # 4. Mask Data
  obs_trait <- true_trait
  missing_idx <- sample(1:n_sp, size = missing_pct * n_sp)
  obs_trait[missing_idx] <- NA
  
  list(
    tree = tree,
    env = env,
    true_trait = true_trait,
    obs_trait = obs_trait,
    missing_idx = missing_idx
  )
}

# ── 3. Helper: Spectral Features (Tree Coordinates) ──────────────────────────
get_spectral_features <- function(tree, k = 8) {
  D <- cophenetic(tree)
  # Heuristic sigma
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 0
  
  # Laplacian
  D_deg <- diag(rowSums(A))
  L <- D_deg - A
  
  # Eigen vectors (smallest k)
  eig <- eigen(L, symmetric = TRUE)
  n <- nrow(A)
  coords <- eig$vectors[, (n - k + 1):n]
  
  torch_tensor(coords, dtype = torch_float(), device = device)
}

# ── 4. Helper: Stable Adjacency ──────────────────────────────────────────────
get_stable_adj <- function(tree) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 0
  
  row_sums <- rowSums(A)
  row_sums[row_sums == 0] <- 1
  A_norm <- diag(1/row_sums) %*% A
  
  torch_tensor(A_norm, dtype = torch_float(), device = device)
}

# ── 5. Robust Benchmark with Correlation ─────────────────────────────────────
run_benchmark <- function(n_sims = 5, n_sp = 2000) {
  
  results <- data.frame()
  cat(sprintf("Running Benchmark (N=%d)...\n", n_sp))
  
  for(i in 1:n_sims) {
    cat(sprintf("\n[Sim %d/%d]\n", i, n_sims))
    
    # 1. Sim Data (Complex Sigmoid)
    tree <- rtree(n_sp)
    # FIX: Add tiny length to branches to stop Rphylopars singular warnings
    tree$edge.length <- tree$edge.length + 1e-5 
    
    env <- rTraitCont(tree, model = "BM", sigma = 0.5)
    env_std <- as.numeric(scale(env))
    phylo_noise <- rTraitCont(tree, model = "BM", sigma = 0.2)
    
    # The Sigmoid Function
    true_trait <- 10 / (1 + exp(-3 * env_std)) + phylo_noise
    
    obs_trait <- true_trait
    missing_idx <- sample(1:n_sp, size = 0.3 * n_sp)
    obs_trait[missing_idx] <- NA
    
    dat <- list(tree=tree, env=env_std, obs_trait=obs_trait, missing_idx=missing_idx)
    
    # 2. Rphylopars
    cat("Running Rphylopars... ")
    t0 <- Sys.time()
    phy_pred <- tryCatch({
      df <- data.frame(species = tree$tip.label, trait = obs_trait, env = env_std)
      suppressWarnings({
        p <- phylopars(df, tree, model="BM", pheno_error=FALSE)
      })
      p$anc_recon[1:n_sp, 1][missing_idx]
    }, error = function(e) {
      cat("Failed. ")
      return(rep(NA, length(missing_idx)))
    })
    cat(sprintf("Done (%.1fs)\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
    
    # 3. Spectral VGAE
    cat("Running VGAE... \n")
    t0 <- Sys.time()
    # Ensure you have the 'impute_vgae_cov' function loaded from previous steps
    gcn_pred_full <- impute_vgae_cov(dat, epochs=1000, lr=0.01) 
    gcn_pred <- gcn_pred_full[missing_idx]
    cat(sprintf("VGAE Done (%.1fs)\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
    
    # 4. Metrics
    truth <- true_trait[missing_idx]
    
    # RMSE
    rmse_phy <- sqrt(mean((truth - phy_pred)^2, na.rm=TRUE))
    rmse_gcn <- sqrt(mean((truth - gcn_pred)^2))
    
    # Correlation (Handling NAs for safety)
    cor_phy <- cor(truth, phy_pred, use = "complete.obs")
    cor_gcn <- cor(truth, gcn_pred, use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Sim = i,
      Method = c("Rphylopars", "VGAE"),
      RMSE = c(rmse_phy, rmse_gcn),
      Cor  = c(cor_phy, cor_gcn)
    ))
  }
  return(results)
}

# ── 6. Execute & Plot ────────────────────────────────────────────────────────
final_res <- run_benchmark(n_sims = 5, n_sp = 2000)

# Print Summary
cat("\n--- FINAL RESULTS (N=2000 Non-Linear) ---\n")
print(aggregate(cbind(RMSE, Cor) ~ Method, data = final_res, mean))

# Visualization
library(ggplot2)
library(patchwork) # Great for combining plots, install if needed

# RMSE Plot
p1 <- ggplot(final_res, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot() +
  geom_jitter(width=0.1, size=3, alpha=0.7) +
  theme_minimal() +
  labs(title = "Imputation Error (RMSE)",
       subtitle = "Lower is Better",
       y = "RMSE") +
  theme(legend.position = "none")

# Correlation Plot
p2 <- ggplot(final_res, aes(x = Method, y = Cor, fill = Method)) +
  geom_boxplot() +
  geom_jitter(width=0.1, size=3, alpha=0.7) +
  theme_minimal() +
  labs(title = "Imputation Correlation",
       subtitle = "Higher is Better",
       y = "Pearson Correlation") +
  theme(legend.position = "none")

# Combine side-by-side
p1 + p2