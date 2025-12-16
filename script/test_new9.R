# ==============================================================================
# PURE PHYLO-DAE (FIXED): denoising GAE with categorical + continuous covariates
# ==============================================================================
# Key fixes vs your original:
#  1) mask token (NOT 0) + mask indicator covariate
#  2) corruption only among truly observed tips (never NA placeholders)
#  3) correct spectral embedding: smallest non-zero Laplacian eigenvectors
#  4) stable adjacency: symmetric normalized kernel + self loops
#  5) optional categorical covariates (one-hot) + env^2 expansion
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(patchwork)

set.seed(1)
torch_manual_seed(1)

# ── Device ───────────────────────────────────────────────────────────────────
device <- if (cuda_is_available()) {
  torch_device("cuda")
} else if (backends_mps_is_available()) {
  torch_device("mps")
} else {
  torch_device("cpu")
}
cat("Running on:", as.character(device), "\n")

# ── Helpers ──────────────────────────────────────────────────────────────────
one_hot <- function(f) {
  f <- factor(f)
  mm <- model.matrix(~ f - 1)
  colnames(mm) <- gsub("^f", "", colnames(mm))
  mm
}

# Correct spectral features: smallest non-zero eigenvectors
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
  coords <- vecs[, 2:(k + 1), drop = FALSE]  # skip constant eigenvector
  torch_tensor(coords, dtype = torch_float(), device = device)
}

# Stable adjacency: sym-normalized + self-loops
get_adj_symnorm <- function(tree) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 1
  
  rs <- rowSums(A) + 1e-8
  Dinv_sqrt <- diag(1 / sqrt(rs))
  A_norm <- Dinv_sqrt %*% A %*% Dinv_sqrt
  torch_tensor(A_norm, dtype = torch_float(), device = device)
}

# ── Pure Phylo-DAE model (mask-aware) ─────────────────────────────────────────
PhyloDAE <- nn_module(
  "PhyloDAE",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim) {
    total_input <- input_dim + coord_dim + cov_dim
    
    self$enc1 <- nn_linear(total_input, hidden_dim)
    self$enc2 <- nn_linear(hidden_dim, hidden_dim)
    
    self$msg  <- nn_linear(hidden_dim, hidden_dim)
    self$alpha <- nn_parameter(torch_tensor(0.10, device = device)) # message strength
    
    self$dec1 <- nn_linear(hidden_dim, hidden_dim)
    self$dec2 <- nn_linear(hidden_dim, input_dim)
    
    self$act  <- nn_relu()
    self$drop <- nn_dropout(0.15)
    
    # learned mask token for x
    self$mask_token <- nn_parameter(torch_zeros(1, 1, device = device))
  },
  
  forward = function(x, coords, covs, adj) {
    combined <- torch_cat(list(x, coords, covs), dim = 2)
    
    h <- self$enc1(combined); h <- self$act(h); h <- self$drop(h)
    h <- self$enc2(h);       h <- self$act(h)
    
    m <- self$msg(torch_matmul(adj, h))
    h <- h + self$alpha * m
    
    h <- self$dec1(h); h <- self$act(h)
    self$dec2(h)
  }
)

# ── Imputation engine ─────────────────────────────────────────────────────────
impute_phylo_dae <- function(sim_data,
                             epochs = 2500,
                             lr = 0.003,
                             k_eigen = 8,
                             corruption_rate = 0.50,
                             refine_steps = 8,
                             hidden_dim = 96) {
  
  trait_vec <- sim_data$obs_trait
  tree      <- sim_data$tree
  env_cont  <- sim_data$env_cont
  env_cat   <- sim_data$env_cat  # can be NULL
  
  is_missing <- is.na(trait_vec)
  
  # Fill missing with mean only as placeholder (training never uses NA targets)
  mu_global <- mean(trait_vec, na.rm = TRUE)
  X_filled <- trait_vec
  X_filled[is_missing] <- mu_global
  
  # scale
  x_mean <- mean(X_filled)
  x_sd <- sd(X_filled); if (x_sd == 0) x_sd <- 1
  X_scaled <- (X_filled - x_mean) / x_sd
  
  # Covariates: env, env^2, (optional) one-hot(cat)
  env1 <- env_cont
  env2 <- env_cont^2
  cov_mat <- cbind(env1, env2)
  
  if (!is.null(env_cat)) {
    cov_mat <- cbind(cov_mat, one_hot(env_cat))
  }
  
  # tensors
  X_full <- torch_tensor(matrix(X_scaled, ncol = 1), dtype = torch_float(), device = device)
  Cov_t  <- torch_tensor(cov_mat, dtype = torch_float(), device = device)
  M_obs  <- torch_tensor(matrix(!is_missing, ncol = 1), dtype = torch_bool(), device = device)
  
  coords_t <- get_spectral_features(tree, k = k_eigen)
  adj_t    <- get_adj_symnorm(tree)
  
  # add mask-indicator feature
  cov_dim <- ncol(cov_mat) + 1
  
  model <- PhyloDAE(input_dim = 1, hidden_dim = hidden_dim, coord_dim = k_eigen, cov_dim = cov_dim)
  model$to(device = device)
  
  opt <- optim_adam(model$parameters, lr = lr, weight_decay = 2e-4)
  
  for (ep in 1:epochs) {
    model$train()
    opt$zero_grad()
    
    # corrupt ONLY among observed entries
    u <- torch_rand_like(X_full)
    masked_bool <- (u < corruption_rate) & M_obs
    
    if (masked_bool$sum()$item() > 0) {
      Mask_t <- masked_bool$to(dtype = torch_float())
      covs_t <- torch_cat(list(Cov_t, Mask_t), dim = 2)
      
      # apply mask token to x where masked
      X_in <- torch_where(masked_bool, model$mask_token$expand_as(X_full), X_full)
      
      pred <- model(X_in, coords_t, covs_t, adj_t)
      
      # denoising loss only on masked observed entries
      loss <- nnf_mse_loss(pred[masked_bool], X_full[masked_bool])
      loss$backward()
      
      nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
      opt$step()
    }
    
    if (ep %% 500 == 0) cat(".")
  }
  cat("\n")
  
  # inference: iterative refinement on true missing entries
  model$eval()
  with_no_grad({
    Mask0 <- torch_zeros_like(X_full)
    covs0 <- torch_cat(list(Cov_t, Mask0), dim = 2)
    
    obs_f  <- M_obs$to(dtype = torch_float())
    miss_f <- (1 - obs_f)
    
    X_iter <- X_full$clone()
    for (t in 1:refine_steps) {
      pred_t <- model(X_iter, coords_t, covs0, adj_t)
      X_iter <- X_full * obs_f + pred_t * miss_f
    }
    
    final_scaled <- as.numeric(X_iter$cpu())
  })
  
  (final_scaled * x_sd) + x_mean
}

# ── Benchmark with categorical covariate in simulation ────────────────────────
run_benchmark_cat <- function(n_sims = 5, n_sp = 200) {
  
  results <- data.frame()
  cat(sprintf("Running pure Phylo-DAE benchmark with categorical covariate (N=%d)\n", n_sp))
  
  for (i in 1:n_sims) {
    cat(sprintf("[Sim %d/%d] ", i, n_sims))
    
    tree <- rtree(n_sp) ; tree$edge.length <- tree$edge.length + 1e-5
    tree <- chronos(tree, lambda = 0)
    
    env_cont <- as.numeric(scale(rTraitCont(tree, model = "BM", sigma = 0.6)))
    env_cat  <- factor(sample(c("A","B","C"), n_sp, replace = TRUE))
    
    # strong category + nonlinearity + phylo effect
    cat_shift <- c(A = 6, B = -6, C = 0)
    mu <- cat_shift[env_cat] + 4 * tanh(1.4 * env_cont) + 2 * sin(2 * env_cont) * (env_cat == "B")
    
    phy_eff <- rTraitCont(tree, model = "BM", sigma = 0.8)
    true_trait <- mu + phy_eff + rnorm(n_sp, 0, 1.5)
    
    # missingness depends on category (harder, and makes cat matter)
    p_miss <- ifelse(env_cat=="A", 0.15, ifelse(env_cat=="B", 0.55, 0.25))
    missing_idx <- which(runif(n_sp) < p_miss)
    
    obs_trait <- true_trait
    obs_trait[missing_idx] <- NA
    
    dat <- list(tree = tree, env_cont = env_cont, env_cat = env_cat, obs_trait = obs_trait)
    
    # Rphylopars baseline: env_cont only (misspecified on purpose)
    cat("Rphylopars... ")
    phy_pred <- tryCatch({
      df <- data.frame(species = tree$tip.label, trait = obs_trait, env = env_cont)
      p <- suppressWarnings(phylopars(df, tree, model = "BM", pheno_error = FALSE))
      as.numeric(p$anc_recon[1:n_sp, 1])[missing_idx]
    }, error = function(e) rep(NA, length(missing_idx)))
    
    # Pure Phylo-DAE with env_cont + env_cat
    cat("Phylo-DAE... ")
    dae_full <- impute_phylo_dae(dat, epochs = 2500, lr = 0.003,
                                 corruption_rate = 0.60, refine_steps = 10,
                                 hidden_dim = 96)
    dae_pred <- dae_full[missing_idx]
    
    truth <- true_trait[missing_idx]
    
    rmse_phy <- sqrt(mean((truth - phy_pred)^2, na.rm = TRUE))
    rmse_dae <- sqrt(mean((truth - dae_pred)^2))
    
    cor_phy <- cor(truth, phy_pred, use = "complete.obs")
    cor_dae <- cor(truth, dae_pred)
    
    # diagnostic: if ~0, your NN is collapsing to the baseline-like solution
    mean_abs_diff <- mean(abs(dae_pred - phy_pred), na.rm = TRUE)
    
    results <- rbind(results, data.frame(
      Sim = i,
      Method = c("Rphylopars (env_cont only)", "Phylo-DAE (env+cat)"),
      RMSE = c(rmse_phy, rmse_dae),
      Cor  = c(cor_phy, cor_dae),
      MeanAbsDiffPred = c(mean_abs_diff, mean_abs_diff)
    ))
    
    cat(sprintf("Done. mean|DAE-Rphy|=%.3f\n", mean_abs_diff))
  }
  
  results
}

# ── Run + plots ───────────────────────────────────────────────────────────────
final_res <- run_benchmark_cat(n_sims = 5, n_sp = 200)

cat("\n--- FINAL SCOREBOARD ---\n")
print(aggregate(cbind(RMSE, Cor) ~ Method, data = final_res, mean))

cat("\nDiagnostic: mean |DAE - Rphylopars| on missing entries (should NOT be ~0)\n")
print(aggregate(MeanAbsDiffPred ~ Sim, data = final_res, mean))

p1 <- ggplot(final_res, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot(alpha=0.6) + geom_jitter(width=0.1) +
  theme_minimal() + theme(legend.position="none") +
  labs(title="RMSE (lower is better)")

p2 <- ggplot(final_res, aes(x = Method, y = Cor, fill = Method)) +
  geom_boxplot(alpha=0.6) + geom_jitter(width=0.1) +
  theme_minimal() + theme(legend.position="none") +
  labs(title="Correlation (higher is better)")

print(p1 + p2)