# ==============================================================================
# PHYLO-DAE (IMPROVED): Phylogenetic Denoising Graph Auto-Encoder
# ==============================================================================
# Improvements vs your original:
#   1) Spectral features: use SMALLEST non-zero Laplacian eigenvectors (correct)
#   2) Masking: learned mask token + explicit mask indicator (no ambiguity vs 0)
#   3) Corruption: only corrupt truly observed entries (never NA placeholders)
#   4) Adjacency: symmetric normalized kernel with self-loops (stable propagation)
#   5) Training stability: lower dropout + gradient clipping + weight decay
#   6) Inference: iterative refinement for missing entries
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(patchwork)

set.seed(1)
torch_manual_seed(1)

# ── 1. Device Configuration ───────────────────────────────────────────────────
device <- if (cuda_is_available()) {
  torch_device("cuda")
} else if (backends_mps_is_available()) {
  torch_device("mps")
} else {
  torch_device("cpu")
}
cat("Running on:", as.character(device), "\n")

# ── 2. Spectral Features (FIXED) ──────────────────────────────────────────────
# Correct spectral embedding uses the SMALLEST non-zero eigenvectors of the Laplacian.
get_spectral_features <- function(tree, k = 8) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 0
  
  L <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric = TRUE)
  
  # eigenvalues are sorted decreasing by default in R::eigen for symmetric=TRUE
  # so we sort ascending ourselves:
  ord <- order(eig$values, decreasing = FALSE)
  vecs <- eig$vectors[, ord, drop = FALSE]
  
  # First eigenvector corresponds to eigenvalue ~0 (constant); skip it.
  # Take next k vectors.
  n <- nrow(A)
  if (k + 1 > n) stop("k too large for number of tips.")
  coords <- vecs[, 2:(k + 1), drop = FALSE]
  
  torch_tensor(coords, dtype = torch_float(), device = device)
}

# ── 3. Stable Adjacency (Symmetric norm + self-loops) ─────────────────────────
get_adj_symnorm <- function(tree) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  
  # add self-loops for stability
  diag(A) <- 1
  
  rs <- rowSums(A) + 1e-8
  Dinv_sqrt <- diag(1 / sqrt(rs))
  A_norm <- Dinv_sqrt %*% A %*% Dinv_sqrt
  
  torch_tensor(A_norm, dtype = torch_float(), device = device)
}

# ── 4. Rphylopars helper: mean + SE for all tips ──────────────────────────────
rphylopars_impute_with_se <- function(tree, obs_trait, env_vec, model = "BM") {
  n_sp <- length(tree$tip.label)
  
  df <- data.frame(
    species = tree$tip.label,
    trait   = obs_trait,
    env     = env_vec
  )
  
  p <- suppressWarnings(phylopars(df, tree, model = model, pheno_error = FALSE))
  
  trait_col <- which(colnames(p$anc_recon) %in% c("trait", "Trait"))[1]
  if (is.na(trait_col)) trait_col <- 1
  
  mu_all  <- as.numeric(p$anc_recon[1:n_sp, trait_col])
  var_all <- as.numeric(p$anc_var[1:n_sp, trait_col])
  se_all  <- sqrt(pmax(var_all, 0))
  
  list(mu = mu_all, se = se_all)
}

# ── 5. Model Architecture (mask-aware denoising) ──────────────────────────────


PhyloDAE <- nn_module(
  "PhyloDAE",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim) {
    total_input <- input_dim + coord_dim + cov_dim
    
    self$enc1 <- nn_linear(total_input, hidden_dim)
    self$enc2 <- nn_linear(hidden_dim, hidden_dim)
    
    self$graph_mix <- nn_linear(hidden_dim, hidden_dim)
    
    self$dec1 <- nn_linear(hidden_dim, hidden_dim)
    self$dec2 <- nn_linear(hidden_dim, input_dim)
    
    self$act  <- nn_relu()
    self$drop <- nn_dropout(0.15)
    
    # learned mask token for x when masked
    self$mask_token <- nn_parameter(torch_zeros(1, 1, device = device))
  },
  
  forward = function(x, coords, covs, adj) {
    combined <- torch_cat(list(x, coords, covs), dim = 2)
    
    h <- self$enc1(combined); h <- self$act(h); h <- self$drop(h)
    h <- self$enc2(h);       h <- self$act(h)
    
    # message passing
    neigh <- torch_matmul(adj, h)
    h <- h + self$graph_mix(neigh)
    
    h <- self$dec1(h); h <- self$act(h)
    out <- self$dec2(h)
    out
  }
)

# ── 6. Imputation Function (correct denoising signal) ─────────────────────────
impute_phylo_dae <- function(sim_data,
                             epochs = 2000,
                             lr = 0.003,
                             k_eigen = 8,
                             corruption_rate = 0.4,
                             refine_steps = 5) {
  
  trait_vec <- sim_data$obs_trait
  tree <- sim_data$tree
  env_vec <- sim_data$env
  
  is_missing <- is.na(trait_vec)
  
  # Fill missing with mean (placeholder), then scale
  mu_global <- mean(trait_vec, na.rm = TRUE)
  X_filled <- trait_vec
  X_filled[is_missing] <- mu_global
  
  x_mean <- mean(X_filled)
  x_sd <- sd(X_filled); if (x_sd == 0) x_sd <- 1
  X_scaled <- (X_filled - x_mean) / x_sd
  
  # Covariates: include a simple nonlinear expansion of env (helps small N)
  env1 <- env_vec
  env2 <- env_vec^2
  env_mat <- cbind(env1, env2)
  
  # tensors
  X_full <- torch_tensor(matrix(X_scaled, ncol = 1), dtype = torch_float(), device = device)
  Env_t  <- torch_tensor(env_mat, dtype = torch_float(), device = device)  # [N,2]
  
  # observed mask (TRUE where we have real observed trait)
  M_obs <- torch_tensor(matrix(!is_missing, ncol = 1), dtype = torch_bool(), device = device)
  
  coords_t <- get_spectral_features(tree, k = k_eigen)
  adj_t    <- get_adj_symnorm(tree)
  
  # covariates include: env1, env2, mask_indicator (changes each epoch)
  cov_dim <- 2 + 1
  
  model <- PhyloDAE(input_dim = 1, hidden_dim = 64, coord_dim = k_eigen, cov_dim = cov_dim)
  model$to(device = device)
  
  opt <- optim_adam(model$parameters, lr = lr, weight_decay = 1e-4)
  
  for (ep in 1:epochs) {
    model$train()
    opt$zero_grad()
    
    # corruption ONLY among observed entries
    u <- torch_rand_like(X_full)
    masked_bool <- (u < corruption_rate) & M_obs
    
    # mask indicator feature (float)
    Mask_t <- masked_bool$to(dtype = torch_float())  # [N,1]
    covs_t <- torch_cat(list(Env_t, Mask_t), dim = 2)  # [N,3]
    
    # apply mask token to x input where masked_bool
    X_in <- X_full$clone()
    if (masked_bool$sum()$item() > 0) {
      X_in <- torch_where(masked_bool,
                          model$mask_token$expand_as(X_full),
                          X_in)
    }
    
    pred <- model(X_in, coords_t, covs_t, adj_t)
    
    # loss ONLY on masked observed entries (true denoising objective)
    if (masked_bool$sum()$item() > 0) {
      loss <- nnf_mse_loss(pred[masked_bool], X_full[masked_bool])
      loss$backward()
      
      # stability
      nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
      
      opt$step()
    }
    
    if (ep %% 500 == 0) cat(".")
  }
  cat("\n")
  
  # ── inference: iterative refinement on TRUE missing entries ──
  model$eval()
  with_no_grad({
    # mask indicator is 0 at inference
    Mask0 <- torch_zeros_like(X_full)
    covs0 <- torch_cat(list(Env_t, Mask0), dim = 2)
    
    obs_f  <- M_obs$to(dtype = torch_float())
    miss_f <- (1 - obs_f)
    
    X_iter <- X_full$clone()
    for (t in 1:refine_steps) {
      pred_t <- model(X_iter, coords_t, covs0, adj_t)
      X_iter <- X_full * obs_f + pred_t * miss_f
    }
    
    final_pred_scaled <- as.numeric(X_iter$cpu())
  })
  
  final_pred <- (final_pred_scaled * x_sd) + x_mean
  final_pred
}

# ── 7. Benchmark ─────────────────────────────────────────────────────────────
run_benchmark <- function(n_sims = 5, n_sp = 200) {
  
  results <- data.frame()
  cat(sprintf("Running Phylo-DAE Benchmark (N=%d)...\n", n_sp))
  
  for (i in 1:n_sims) {
    cat(sprintf("[Sim %d/%d] ", i, n_sims))
    
    tree <- rtree(n_sp)
    tree$edge.length <- tree$edge.length + 1e-5
    
    env <- rTraitCont(tree, model = "BM", sigma = 0.5)
    env_std <- as.numeric(scale(env))
    
    phylo_noise <- rTraitCont(tree, model = "BM", sigma = 0.5)
    true_trait <- 10 / (1 + exp(-1.2 * env_std)) + phylo_noise
    
    obs_trait <- true_trait
    missing_idx <- sample(1:n_sp, size = 0.3 * n_sp)
    obs_trait[missing_idx] <- NA
    
    dat <- list(tree = tree, env = env_std, obs_trait = obs_trait, missing_idx = missing_idx)
    
    # Rphylopars baseline
    cat("Rphylopars... ")
    phy_pred <- tryCatch({
      df <- data.frame(species = tree$tip.label, trait = obs_trait, env = env_std)
      suppressWarnings({
        p <- phylopars(df, tree, model = "BM", pheno_error = FALSE)
      })
      p$anc_recon[1:n_sp, 1][missing_idx]
    }, error = function(e) rep(NA, length(missing_idx)))
    
    # Phylo-DAE
    cat("Phylo-DAE... ")
    dae_full <- impute_phylo_dae(dat, epochs = 2000, lr = 0.003, corruption_rate = 0.45)
    dae_pred <- dae_full[missing_idx]
    
    truth <- true_trait[missing_idx]
    
    rmse_phy <- sqrt(mean((truth - phy_pred)^2, na.rm = TRUE))
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
  
  results
}

# ── 8. Execute & Visualize ───────────────────────────────────────────────────
final_res <- run_benchmark(n_sims = 5, n_sp = 200)

cat("\n--- FINAL SCOREBOARD ---\n")
print(aggregate(cbind(RMSE, Cor) ~ Method, data = final_res, mean))

p1 <- ggplot(final_res, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot(alpha = 0.6) + geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "Imputation Error (RMSE)", subtitle = "Lower is Better") +
  theme(legend.position = "none")

p2 <- ggplot(final_res, aes(x = Method, y = Cor, fill = Method)) +
  geom_boxplot(alpha = 0.6) + geom_jitter(width = 0.1) +
  theme_minimal() +
  labs(title = "Correlation (R)", subtitle = "Higher is Better") +
  theme(legend.position = "none")

tryCatch({ p1 + p2 }, error = function(e) { print(p1); print(p2) })