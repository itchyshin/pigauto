# ==============================================================================
# RESIDUAL PHYLO-DAE (FULL SCRIPT): Phylogenetic Denoising Graph Auto-Encoder
# ==============================================================================
# Key idea: NEVER try to beat Rphylopars from scratch.
#          Instead, use Rphylopars as a strong baseline (mu_phy) and train a NN
#          to learn ONLY the residual correction delta = trait - mu_phy.
#
# Improvements included:
#   1) Correct spectral coords: smallest non-zero Laplacian eigenvectors
#   2) Stable adjacency: symmetric normalization + self-loops
#   3) Mask token + explicit mask indicator feature
#   4) Corrupt only truly observed entries (never NA placeholders)
#   5) Residual learning on top of Rphylopars
#   6) Shrink residuals toward 0 (won't drift away from baseline)
#   7) Iterative refinement for missing entries
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(patchwork)

set.seed(1)
torch_manual_seed(1)

# ── 1. Device ────────────────────────────────────────────────────────────────
device <- if (cuda_is_available()) {
  torch_device("cuda")
} else if (backends_mps_is_available()) {
  torch_device("mps")
} else {
  torch_device("cpu")
}
cat("Running on:", as.character(device), "\n")

# ── 2. Spectral Features (correct: smallest non-zero eigenvectors) ────────────
get_spectral_features <- function(tree, k = 8) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 0
  
  L <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric = TRUE)
  
  ord <- order(eig$values, decreasing = FALSE)
  vecs <- eig$vectors[, ord, drop = FALSE]
  
  n <- nrow(A)
  if (k + 1 > n) stop("k too large for number of tips.")
  coords <- vecs[, 2:(k + 1), drop = FALSE]
  
  torch_tensor(coords, dtype = torch_float(), device = device)
}

# ── 3. Stable adjacency (symmetric norm + self-loops) ─────────────────────────
get_adj_symnorm <- function(tree) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  
  diag(A) <- 1  # self-loops
  
  rs <- rowSums(A) + 1e-8
  Dinv_sqrt <- diag(1 / sqrt(rs))
  A_norm <- Dinv_sqrt %*% A %*% Dinv_sqrt
  
  torch_tensor(A_norm, dtype = torch_float(), device = device)
}

# ── 4. Baseline helper (Rphylopars mean for ALL tips) ─────────────────────────
rphylopars_mu_all <- function(tree, trait_vec, env_vec, model = "BM") {
  n_sp <- length(tree$tip.label)
  df <- data.frame(species = tree$tip.label, trait = trait_vec, env = env_vec)
  p <- suppressWarnings(phylopars(df, tree, model = model, pheno_error = FALSE))
  mu <- as.numeric(p$anc_recon[1:n_sp, 1])
  mu
}

# ── 5. Residual Phylo-DAE model ───────────────────────────────────────────────
ResidualPhyloDAE <- nn_module(
  "ResidualPhyloDAE",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim) {
    total_input <- input_dim + coord_dim + cov_dim
    
    self$enc1 <- nn_linear(total_input, hidden_dim)
    self$enc2 <- nn_linear(hidden_dim, hidden_dim)
    
    # very mild message passing (safe)
    self$msg <- nn_linear(hidden_dim, hidden_dim)
    self$alpha <- nn_parameter(torch_tensor(0.05, device = device))
    
    self$dec1 <- nn_linear(hidden_dim, hidden_dim)
    self$dec2 <- nn_linear(hidden_dim, 1)  # delta
    
    self$act  <- nn_relu()
    self$drop <- nn_dropout(0.10)
    
    self$mask_token <- nn_parameter(torch_zeros(1, 1, device = device))
  },
  
  forward = function(x, coords, covs, adj) {
    combined <- torch_cat(list(x, coords, covs), dim = 2)
    
    h <- self$enc1(combined); h <- self$act(h); h <- self$drop(h)
    h <- self$enc2(h);       h <- self$act(h)
    
    m <- torch_matmul(adj, h)
    h <- h + self$alpha * self$msg(m)
    
    h <- self$dec1(h); h <- self$act(h)
    delta <- self$dec2(h)
    delta
  }
)

# ── 6. Residual imputer ───────────────────────────────────────────────────────
impute_phylo_dae_residual <- function(sim_data,
                                      epochs = 3000,
                                      lr = 0.003,
                                      k_eigen = 8,
                                      corruption_rate = 0.55,
                                      refine_steps = 8,
                                      lambda_shrink = 0.03,
                                      hidden_dim = 64,
                                      rphylopars_model = "BM") {
  
  trait_vec <- sim_data$obs_trait
  tree <- sim_data$tree
  env_vec <- sim_data$env
  n_sp <- length(tree$tip.label)
  
  is_missing <- is.na(trait_vec)
  
  # ---- baseline mean mu_phy for all tips ----
  mu_phy <- tryCatch({
    rphylopars_mu_all(tree, trait_vec, env_vec, model = rphylopars_model)
  }, error = function(e) {
    rep(mean(trait_vec, na.rm = TRUE), n_sp)
  })
  
  # fill missing with baseline
  X_filled <- trait_vec
  X_filled[is_missing] <- mu_phy[is_missing]
  
  # scaling (shared)
  x_mean <- mean(X_filled)
  x_sd <- sd(X_filled); if (x_sd == 0) x_sd <- 1
  
  X_scaled  <- (X_filled - x_mean) / x_sd
  MU_scaled <- (mu_phy  - x_mean) / x_sd
  DELTA_true <- X_scaled - MU_scaled
  
  # covariates: env, env^2, MU_scaled
  env1 <- env_vec
  env2 <- env_vec^2
  cov_mat <- cbind(env1, env2, MU_scaled)
  
  # tensors
  X_full    <- torch_tensor(matrix(X_scaled, ncol = 1), dtype = torch_float(), device = device)
  MU_full   <- torch_tensor(matrix(MU_scaled, ncol = 1), dtype = torch_float(), device = device)
  DELTA_tgt <- torch_tensor(matrix(DELTA_true, ncol = 1), dtype = torch_float(), device = device)
  Cov_t     <- torch_tensor(cov_mat, dtype = torch_float(), device = device)  # [N,3]
  
  M_obs <- torch_tensor(matrix(!is_missing, ncol = 1), dtype = torch_bool(), device = device)
  
  coords_t <- get_spectral_features(tree, k = k_eigen)
  adj_t    <- get_adj_symnorm(tree)
  
  # cov_dim = (env1, env2, MU_scaled, mask_indicator) = 4
  model <- ResidualPhyloDAE(input_dim = 1, hidden_dim = hidden_dim, coord_dim = k_eigen, cov_dim = 4)
  model$to(device = device)
  opt <- optim_adam(model$parameters, lr = lr, weight_decay = 1e-4)
  
  # ---- training ----
  for (ep in 1:epochs) {
    model$train()
    opt$zero_grad()
    
    # corrupt only observed entries
    u <- torch_rand_like(X_full)
    masked_bool <- (u < corruption_rate) & M_obs
    
    if (masked_bool$sum()$item() > 0) {
      Mask_t <- masked_bool$to(dtype = torch_float())               # [N,1]
      covs_t <- torch_cat(list(Cov_t, Mask_t), dim = 2)             # [N,4]
      
      X_in <- X_full$clone()
      X_in <- torch_where(masked_bool, model$mask_token$expand_as(X_full), X_in)
      
      delta_hat <- model(X_in, coords_t, covs_t, adj_t)
      
      loss_denoise <- nnf_mse_loss(delta_hat[masked_bool], DELTA_tgt[masked_bool])
      loss_shrink  <- nnf_mse_loss(delta_hat[masked_bool], torch_zeros_like(delta_hat[masked_bool]))
      loss <- loss_denoise + lambda_shrink * loss_shrink
      
      loss$backward()
      nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
      opt$step()
    }
    
    if (ep %% 500 == 0) cat(".")
  }
  cat("\n")
  
  # ---- inference (iterative refinement) ----
  model$eval()
  with_no_grad({
    Mask0 <- torch_zeros_like(X_full)
    covs0 <- torch_cat(list(Cov_t, Mask0), dim = 2)
    
    obs_f  <- M_obs$to(dtype = torch_float())
    miss_f <- (1 - obs_f)
    
    X_iter <- X_full$clone()
    for (t in 1:refine_steps) {
      delta_hat <- model(X_iter, coords_t, covs0, adj_t)
      pred_scaled <- MU_full + delta_hat
      X_iter <- X_full * obs_f + pred_scaled * miss_f
    }
    
    final_pred_scaled <- as.numeric(X_iter$cpu())
  })
  
  final_pred <- (final_pred_scaled * x_sd) + x_mean
  final_pred
}

# ── 7. Benchmark ─────────────────────────────────────────────────────────────
run_benchmark <- function(n_sims = 5, n_sp = 200,
                          sigmoid_beta = 1.2,
                          phylo_sigma = 0.5,
                          missing_rate = 0.3) {
  
  results <- data.frame()
  cat(sprintf("Running Residual Phylo-DAE Benchmark (N=%d)...\n", n_sp))
  
  for (i in 1:n_sims) {
    cat(sprintf("[Sim %d/%d] ", i, n_sims))
    
    tree <- rtree(n_sp)
    tree$edge.length <- tree$edge.length + 1e-5
    
    env <- rTraitCont(tree, model = "BM", sigma = 0.2)
    env_std <- as.numeric(scale(env))
    
    phylo_noise <- rTraitCont(tree, model = "BM", sigma = phylo_sigma)
    true_trait <- 10 / (1 + exp(-sigmoid_beta * env_std)) + phylo_noise
    
    obs_trait <- true_trait
    missing_idx <- sample(1:n_sp, size = missing_rate * n_sp)
    obs_trait[missing_idx] <- NA
    
    dat <- list(tree = tree, env = env_std, obs_trait = obs_trait, missing_idx = missing_idx)
    
    # ---- Rphylopars baseline (mean) ----
    cat("Rphylopars... ")
    phy_pred <- tryCatch({
      df <- data.frame(species = tree$tip.label, trait = obs_trait, env = env_std)
      suppressWarnings({
        p <- phylopars(df, tree, model = "BM", pheno_error = FALSE)
      })
      p$anc_recon[1:n_sp, 1][missing_idx]
    }, error = function(e) rep(NA, length(missing_idx)))
    
    # ---- Residual Phylo-DAE ----
    cat("Residual Phylo-DAE... ")
    dae_full <- impute_phylo_dae_residual(
      dat,
      epochs = 3000,
      lr = 0.003,
      k_eigen = 8,
      corruption_rate = 0.55,
      refine_steps = 8,
      lambda_shrink = 0.03,
      hidden_dim = 64,
      rphylopars_model = "BM"
    )
    dae_pred <- dae_full[missing_idx]
    
    truth <- true_trait[missing_idx]
    
    rmse_phy <- sqrt(mean((truth - phy_pred)^2, na.rm = TRUE))
    rmse_dae <- sqrt(mean((truth - dae_pred)^2))
    
    cor_phy <- cor(truth, phy_pred, use = "complete.obs")
    cor_dae <- cor(truth, dae_pred, use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Sim = i,
      Method = c("Rphylopars", "Residual Phylo-DAE"),
      RMSE = c(rmse_phy, rmse_dae),
      Cor  = c(cor_phy, cor_dae)
    ))
    cat("Done.\n")
  }
  
  results
}

# ── 8. Run two scenarios (like your Set A and Set B) ──────────────────────────
# Set A: strong phylo noise, milder nonlinearity (harder to beat Rphylopars)
res_A <- run_benchmark(n_sims = 5, n_sp = 200, sigmoid_beta = 1.2, phylo_sigma = 0.5)

cat("\n--- SCOREBOARD: Set A (beta=1.2, phylo_sigma=0.5) ---\n")
print(aggregate(cbind(RMSE, Cor) ~ Method, data = res_A, mean))

pA1 <- ggplot(res_A, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot(alpha = 0.6) + geom_jitter(width = 0.1) +
  theme_minimal() + labs(title = "Set A RMSE") + theme(legend.position = "none")

pA2 <- ggplot(res_A, aes(x = Method, y = Cor, fill = Method)) +
  geom_boxplot(alpha = 0.6) + geom_jitter(width = 0.1) +
  theme_minimal() + labs(title = "Set A Correlation") + theme(legend.position = "none")

print(pA1 + pA2)

# Set B: weaker phylo noise, stronger nonlinearity (good for NN correction)
res_B <- run_benchmark(n_sims = 5, n_sp = 200, sigmoid_beta = 3.0, phylo_sigma = 0.2)

cat("\n--- SCOREBOARD: Set B (beta=3.0, phylo_sigma=0.2) ---\n")
print(aggregate(cbind(RMSE, Cor) ~ Method, data = res_B, mean))

pB1 <- ggplot(res_B, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot(alpha = 0.6) + geom_jitter(width = 0.1) +
  theme_minimal() + labs(title = "Set B RMSE") + theme(legend.position = "none")

pB2 <- ggplot(res_B, aes(x = Method, y = Cor, fill = Method)) +
  geom_boxplot(alpha = 0.6) + geom_jitter(width = 0.1) +
  theme_minimal() + labs(title = "Set B Correlation") + theme(legend.position = "none")

print(pB1 + pB2)