# ==============================================================================
# PHYLO-DAE: Phylogenetic Denoising Graph Auto-Encoder
# ==============================================================================
# Goal: Impute missing trait values for species on a phylogenetic tree.
# Strategy:
#   1. Input: Trait + Environment + "Tree Coordinates" (Spectral Features).
#   2. Denoising: Randomly set 20% of KNOWN data to 0 during training.
#   3. Graph Mixing: Allow species to "borrow" info from neighbors.
#   4. Loss: Calculate error ONLY on the data we deliberately corrupted.
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars) # For comparison baseline
library(ggplot2)
library(patchwork)  # For arranging plots (install.packages("patchwork"))

# ── 1. Device Configuration ───────────────────────────────────────────────────
# Detects if you have an NVIDIA GPU (cuda) or Mac GPU (mps). Defaults to CPU.
device <- if (cuda_is_available()) {
  torch_device("cuda")
} else if (backends_mps_is_available()) {
  torch_device("mps") # Mac M1/M2/M3 acceleration
} else {
  torch_device("cpu")
}
cat("Running on:", as.character(device), "\n")

# ── 2. Feature Engineering: "Tree Coordinates" ────────────────────────────────
# Why? Neural networks don't understand "trees" natively.
# We turn the tree into "Spectral Coordinates" (Eigenvectors of the Laplacian).
# This gives every species a set of numbers [c1, c2, c3...] representing its
# location in the evolutionary tree. Close relatives have similar coordinates.
get_spectral_features <- function(tree, k = 8) {
  # 1. Get evolutionary distances (Cophenetic distance)
  D <- cophenetic(tree)
  
  # 2. Convert Distance -> Similarity (Gaussian Kernel)
  # "sigma" controls how fast similarity decays. 
  # 0.5 * median ensures the scale matches the tree's size.
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 0 # No self-loops (a species isn't its own neighbor)
  
  # 3. Calculate Graph Laplacian (L = Degree - Adjacency)
  # This matrix mathematically captures how "clumpy" the tree is.
  L <- diag(rowSums(A)) - A
  
  # 4. Eigen-Decomposition
  # The eigenvectors of L are the "coordinates" of the graph.
  eig <- eigen(L, symmetric = TRUE)
  n <- nrow(A)
  
  # 5. Select the 'k' smallest non-zero eigenvectors.
  # These represent the broadest clusters in the tree (e.g., "Mammals vs Birds")
  # rather than tiny noise.
  coords <- eig$vectors[, (n - k + 1):n] 
  
  # Return as Torch Tensor
  torch_tensor(coords, dtype = torch_float(), device = device)
}

# ── 3. Graph Adjacency: The "Message Passing" Map ─────────────────────────────
# This matrix tells the Neural Network who talks to whom.
get_stable_adj <- function(tree) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 0
  
  # Row-Normalize:
  # Crucial! Ensures that a species with 100 neighbors doesn't get
  # 100x more signal than a species with 1 neighbor.
  row_sums <- rowSums(A)
  row_sums[row_sums == 0] <- 1 # Safety to prevent division by zero
  A_norm <- diag(1/row_sums) %*% A
  
  torch_tensor(A_norm, dtype = torch_float(), device = device)
}

# ── 4. The Model Architecture ─────────────────────────────────────────────────
PhyloDAE <- nn_module(
  "PhyloDAE",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim) {
    # We combine 3 inputs:
    # 1. The Trait (maybe missing/corrupted)
    # 2. The Tree Coordinates (Spectral features)
    # 3. The Environment (Temperature, Rainfall, etc.)
    total_input <- input_dim + coord_dim + cov_dim
    
    # ── Encoder (Compression) ──
    # Wide hidden layer (128 units) to capture complex non-linear curves
    self$enc1 <- nn_linear(total_input, hidden_dim)
    self$enc2 <- nn_linear(hidden_dim, hidden_dim)
    
    # ── Graph Mixing Layer (The "Phylo" Magic) ──
    # This layer allows the latent representation of species i
    # to be updated by the weighted average of its relatives.
    self$graph_mix <- nn_linear(hidden_dim, hidden_dim)
    
    # ── Decoder (Reconstruction) ──
    # Expands back from Hidden Dimension -> 1 (The Trait)
    self$dec1 <- nn_linear(hidden_dim, hidden_dim)
    self$dec2 <- nn_linear(hidden_dim, input_dim)
    
    # Activation & Regularization
    self$act  <- nn_relu()       # Allows learning curves (non-linear)
    self$drop <- nn_dropout(0.5) # 50% Dropout: Essential for Denoising
  },
  
  forward = function(x, coords, covs, adj) {
    # 1. Stack all inputs together [N, Total_Features]
    combined <- torch_cat(list(x, coords, covs), dim = 2)
    
    # 2. Encoder Pass
    h <- self$enc1(combined)
    h <- self$act(h)
    h <- self$drop(h) # Randomly kill neurons to force robustness
    
    h <- self$enc2(h)
    h <- self$act(h)
    
    # 3. Graph Mixing (Residual Connection)
    # Formula: H_new = H_old + Weights * (Adjacency * H_old)
    # Interpretation: "Keep my own personality (H_old), but add a 
    # flavor of my family's personality (Adj * H)."
    neighbor_info <- torch_matmul(adj, h)
    h_mix <- self$graph_mix(neighbor_info)
    h <- h + h_mix 
    
    # 4. Decoder Pass
    h <- self$dec1(h)
    h <- self$act(h)
    out <- self$dec2(h) # Final prediction
    
    return(out)
  }
)

# ── 5. The Training Engine ────────────────────────────────────────────────────
impute_phylo_dae <- function(sim_data, epochs=2000, lr=0.001, k_eigen=8, corruption_rate=0.2) {
  
  # Unpack simulation data
  trait_vec <- sim_data$obs_trait
  tree <- sim_data$tree
  env_vec <- sim_data$env
  
  # ── A) Pre-Processing ──
  # 1. Identify where data is TRULY missing (we never train on this)
  is_missing <- is.na(trait_vec)
  
  # 2. Initial Fill: Replace NA with Mean (just as a placeholder)
  mu_global <- mean(trait_vec, na.rm=TRUE)
  X_filled <- trait_vec
  X_filled[is_missing] <- mu_global
  
  # 3. Z-Score Scaling: Neural Nets fail if data isn't ~N(0,1)
  x_mean <- mean(X_filled)
  x_sd <- sd(X_filled)
  if(x_sd == 0) x_sd <- 1
  X_scaled <- (X_filled - x_mean) / x_sd
  
  # ── B) Tensor Conversion ──
  # Note: reshaping to [N, 1] is critical for matrix multiplication
  X_full <- torch_tensor(matrix(X_scaled, ncol=1), dtype=torch_float(), device=device)
  Env_t  <- torch_tensor(matrix(env_vec, ncol=1), dtype=torch_float(), device=device)
  
  # The Mask of "What we actually know"
  M_true_obs <- torch_tensor(matrix(!is_missing, ncol=1), dtype=torch_bool(), device=device)
  
  adj_t <- get_stable_adj(tree)
  coords_t <- get_spectral_features(tree, k = k_eigen)
  
  # ── C) Model Initialization ──
  model <- PhyloDAE(input_dim=1, hidden_dim=128, coord_dim=k_eigen, cov_dim=1)
  model$to(device = device)
  
  # Adam Optimizer: Standard for Deep Learning
  opt <- optim_adam(model$parameters, lr = lr)
  
  # ── D) The Denoising Loop ──
  for(i in 1:epochs) {
    model$train()
    opt$zero_grad()
    
    # 1. DYNAMIC CORRUPTION (The "Denoising" Part)
    # Every epoch, we create a NEW random mask.
    # We select 'corruption_rate' (e.g., 20%) of the data to HIDE.
    noise_mask <- torch_rand_like(X_full) > corruption_rate
    
    # 2. Corrupt the Input
    # multiply by mask: 1=keep, 0=destroy. 
    # The model sees ZEROS where there used to be data.
    X_corrupted <- X_full * noise_mask
    
    # 3. Forward Pass
    # The model tries to predict the trait for everyone
    pred <- model(X_corrupted, coords_t, Env_t, adj_t)
    
    # 4. Calculate Loss
    # CRITICAL: We only penalize the model on the data we JUST HID.
    # We say: "I showed you 0, but I know the answer is 1.5. Learn that connection!"
    # Target Mask = (Truly Observed Data) AND (Data we artificially corrupted)
    target_mask <- M_true_obs & (!noise_mask)
    
    if (target_mask$sum()$item() > 0) {
      # MSE Loss
      loss <- nnf_mse_loss(pred[target_mask], X_full[target_mask])
      loss$backward()
      opt$step()
    }
    
    # Print progress every 500 epochs
    if(i %% 500 == 0) cat(".")
  }
  cat("\n")
  
  # ── E) Final Prediction ──
  model$eval()
  # For the final answer, we feed the FULL uncorrupted data
  final_pred_scaled <- as.numeric(model(X_full, coords_t, Env_t, adj_t)$cpu())
  
  # Un-scale back to original units (e.g., mm, kg)
  final_pred <- (final_pred_scaled * x_sd) + x_mean
  
  return(final_pred)
}

# ── 6. Benchmark: N=200 Species, Non-Linear Sigmoid ──────────────────────────
run_benchmark <- function(n_sims = 5, n_sp = 200) {
  
  results <- data.frame()
  cat(sprintf("Running Phylo-DAE Benchmark (N=%d)...\n", n_sp))
  
  for(i in 1:n_sims) {
    cat(sprintf("[Sim %d/%d] ", i, n_sims))
    
    # 1. Simulate "S-Curve" Data
    # This is non-linear. Rphylopars (Linear) hates this.
    tree <- rtree(n_sp)
    tree$edge.length <- tree$edge.length + 1e-5 # prevent singular matrix crash
    
    env <- rTraitCont(tree, model = "BM", sigma = 0.5)
    env_std <- as.numeric(scale(env))
    phylo_noise <- rTraitCont(tree, model = "BM", sigma = 0.2)
    
    true_trait <- 10 / (1 + exp(-3 * env_std)) + phylo_noise
    
    obs_trait <- true_trait
    missing_idx <- sample(1:n_sp, size = 0.3 * n_sp)
    obs_trait[missing_idx] <- NA
    
    dat <- list(tree=tree, env=env_std, obs_trait=obs_trait, missing_idx=missing_idx)
    
    # 2. Method A: Rphylopars (The Baseline)
    cat("Rphylopars... ")
    phy_pred <- tryCatch({
      df <- data.frame(species = tree$tip.label, trait = obs_trait, env = env_std)
      suppressWarnings({
        p <- phylopars(df, tree, model="BM", pheno_error=FALSE)
      })
      p$anc_recon[1:n_sp, 1][missing_idx]
    }, error = function(e) rep(NA, length(missing_idx)))
    
    # 3. Method B: Phylo-DAE (Our Model)
    cat("Phylo-DAE... ")
    dae_full <- impute_phylo_dae(dat, epochs=1500, lr=0.005, corruption_rate=0.2)
    dae_pred <- dae_full[missing_idx]
    
    # 4. Calculate Metrics (RMSE and Correlation)
    truth <- true_trait[missing_idx]
    
    rmse_phy <- sqrt(mean((truth - phy_pred)^2, na.rm=TRUE))
    rmse_dae <- sqrt(mean((truth - dae_pred)^2))
    
    cor_phy <- cor(truth, phy_pred, use = "complete.obs")
    cor_dae <- cor(truth, dae_pred, use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Sim = i,
      Method = c("Rphylopars", "Phylo-DAE"),
      RMSE = c(rmse_phy, rmse_dae),
      Cor  = c(cor_phy, cor_dae)
    ))
    cat("Done.\n")
  }
  return(results)
}

# ── 7. Execute & Visualize ───────────────────────────────────────────────────
final_res <- run_benchmark(n_sims = 5, n_sp = 200)

cat("\n--- FINAL SCOREBOARD ---\n")
print(aggregate(cbind(RMSE, Cor) ~ Method, data = final_res, mean))

# Plots
p1 <- ggplot(final_res, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot(alpha=0.6) + geom_jitter(width=0.1) +
  theme_minimal() + 
  labs(title = "Imputation Error (RMSE)", subtitle = "Lower is Better") +
  theme(legend.position = "none")

p2 <- ggplot(final_res, aes(x = Method, y = Cor, fill = Method)) +
  geom_boxplot(alpha=0.6) + geom_jitter(width=0.1) +
  theme_minimal() + 
  labs(title = "Correlation (R)", subtitle = "Higher is Better") +
  theme(legend.position = "none")

# Display side-by-side
tryCatch({ p1 + p2 }, error = function(e) { print(p1); print(p2) })
