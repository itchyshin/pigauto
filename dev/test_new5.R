# ==============================================================================
# CALIBRATED PHYLO-DAE: Aleatoric Uncertainty & Non-Linear Imputation
# ==============================================================================
# 1. Predicts both Mean (Value) and Variance (Uncertainty) per species.
# 2. Uses Heteroscedastic Loss (NLL) to calibrate confidence intervals.
# 3. Beats Linear Baselines (Rphylopars) on non-linear biological data.
# 4. Adds Rphylopars SE (sqrt(anc_var)) and compares uncertainty ribbons visually.
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(patchwork)

# ── 1. Device Setup ──────────────────────────────────────────────────────────
device <- if (cuda_is_available()) {
  torch_device("cuda")
} else if (backends_mps_is_available()) {
  torch_device("mps") # Mac Metal
} else {
  torch_device("cpu")
}
cat("Using device:", as.character(device), "\n")

# ── 2. Feature Engineering ───────────────────────────────────────────────────
get_spectral_features <- function(tree, k = 8) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2)); diag(A) <- 0
  L <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric = TRUE)
  n <- nrow(A)
  coords <- eig$vectors[, (n - k + 1):n]
  torch_tensor(coords, dtype = torch_float(), device = device)
}

get_stable_adj <- function(tree) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2)); diag(A) <- 0
  row_sums <- rowSums(A); row_sums[row_sums == 0] <- 1
  A_norm <- diag(1/row_sums) %*% A
  torch_tensor(A_norm, dtype = torch_float(), device = device)
}

# ── 2b. Rphylopars helper: mean + SE for all tips ────────────────────────────
rphylopars_impute_with_se <- function(tree, obs_trait, env_vec, model = "BM") {
  n_sp <- length(tree$tip.label)
  
  df <- data.frame(
    species = tree$tip.label,
    trait   = obs_trait,
    env     = env_vec
  )
  
  p <- suppressWarnings(phylopars(df, tree, model = model, pheno_error = FALSE))
  
  # Find trait column robustly (fallback to first col)
  trait_col <- which(colnames(p$anc_recon) %in% c("trait", "Trait"))[1]
  if (is.na(trait_col)) trait_col <- 1
  
  mu_all  <- as.numeric(p$anc_recon[1:n_sp, trait_col])
  var_all <- as.numeric(p$anc_var[1:n_sp, trait_col])
  se_all  <- sqrt(pmax(var_all, 0))  # guard tiny negatives due to numerics
  
  list(mu = mu_all, se = se_all)
}

# ── 3. The Calibrated ResNet Model (Improved Architecture) ───────────────────
PhyloDAE_Calibrated <- nn_module(
  "PhyloDAE_Calibrated",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim) {
    total_input <- input_dim + coord_dim + cov_dim
    
    # ── THE DEEP PART (Non-Linear) ──
    self$enc1 <- nn_linear(total_input, hidden_dim)
    self$enc2 <- nn_linear(hidden_dim, hidden_dim)
    self$graph_mix <- nn_linear(hidden_dim, hidden_dim)
    self$dec1 <- nn_linear(hidden_dim, hidden_dim)
    
    # Heads for the Deep Path
    self$head_mu_deep <- nn_linear(hidden_dim, input_dim)
    self$head_var <- nn_linear(hidden_dim, input_dim) # Uncertainty comes from deep features
    
    # ── THE SHORTCUT (Linear Baseline) ──
    # This layer mimics Rphylopars. It connects Input -> Output directly.
    self$linear_shortcut <- nn_linear(total_input, input_dim)
    
    self$act  <- nn_relu()
    self$drop <- nn_dropout(0.2) # Reduced dropout (0.5 was too high for N=200)
  },
  
  forward = function(x, coords, covs, adj) {
    combined <- torch_cat(list(x, coords, covs), dim = 2)
    
    # Path A: The Linear Shortcut (Simple & Robust)
    linear_pred <- self$linear_shortcut(combined)
    
    # Path B: The Deep Network (Complex & Flexible)
    h <- self$enc1(combined); h <- self$act(h); h <- self$drop(h)
    h <- self$enc2(h); h <- self$act(h)
    
    # Graph Mixing
    h_mix <- self$graph_mix(torch_matmul(adj, h))
    h <- h + h_mix
    
    # Decode
    h <- self$dec1(h); h <- self$act(h)
    deep_resid <- self$head_mu_deep(h) # This is the "correction"
    
    # ── FINAL PREDICTION ──
    # Result = Linear Guess + Non-Linear Correction
    mu <- linear_pred + deep_resid
    
    # Uncertainty
    logvar <- self$head_var(h)
    
    return(list(mu = mu, logvar = logvar))
  }
)
# ── 4. The Heteroscedastic Loss Function (NLL) ───────────────────────────────
gaussian_nll_loss <- function(pred_mu, pred_logvar, target, mask) {
  # NLL = 0.5 * [ exp(-logvar)*(y-y_hat)^2 + logvar ]
  precision <- torch_exp(-pred_logvar)
  squared_error <- (pred_mu - target)^2
  loss <- 0.5 * (precision * squared_error + pred_logvar)
  
  mask_f <- mask$to(dtype = torch_float())
  torch_sum(loss * mask_f) / torch_sum(mask_f + 1e-6)
}

# ── 5. Main Imputation Function ──────────────────────────────────────────────
impute_calibrated <- function(sim_data, epochs = 2000, lr = 0.005, k_eigen = 8,
                              corruption_rate = 0.3) {
  
  trait_vec <- sim_data$obs_trait
  tree <- sim_data$tree
  env_vec <- sim_data$env
  
  # A) Data Prep
  is_missing <- is.na(trait_vec)
  mu_global <- mean(trait_vec, na.rm = TRUE)
  X_filled <- trait_vec
  X_filled[is_missing] <- mu_global
  
  x_mean <- mean(X_filled)
  x_sd <- sd(X_filled)
  if (x_sd == 0) x_sd <- 1
  X_scaled <- (X_filled - x_mean) / x_sd
  
  # Tensors (Forced [N, 1])
  X_full <- torch_tensor(matrix(X_scaled, ncol = 1), dtype = torch_float(), device = device)
  Env_t  <- torch_tensor(matrix(env_vec, ncol = 1), dtype = torch_float(), device = device)
  M_true_obs <- torch_tensor(matrix(!is_missing, ncol = 1), dtype = torch_bool(), device = device)
  
  adj_t <- get_stable_adj(tree)
  coords_t <- get_spectral_features(tree, k = k_eigen)
  
  # B) Model
  model <- PhyloDAE_Calibrated(input_dim = 1, hidden_dim = 128, coord_dim = k_eigen, cov_dim = 1)
  model$to(device = device)
  opt <- optim_adam(model$parameters, lr = lr)
  
  # C) Training
  for (i in 1:epochs) {
    model$train()
    opt$zero_grad()
    
    # Dynamic Denoising
    noise_mask <- torch_rand_like(X_full) > corruption_rate
    X_corrupted <- X_full * noise_mask
    
    # Forward Pass
    out <- model(X_corrupted, coords_t, Env_t, adj_t)
    
    # NLL Loss on hidden data (observed but masked by corruption)
    target_mask <- M_true_obs & (!noise_mask)
    if (target_mask$sum()$item() > 0) {
      loss <- gaussian_nll_loss(out$mu, out$logvar, X_full, target_mask)
      loss$backward()
      opt$step()
    }
    
    if (i %% 500 == 0) cat(".")
  }
  cat("\n")
  
  # D) Inference
  model$eval()
  with_no_grad({
    final_out <- model(X_full, coords_t, Env_t, adj_t)
    
    # Mean
    pred_mu_scaled <- as.numeric(final_out$mu$cpu())
    final_imputed <- (pred_mu_scaled * x_sd) + x_mean
    
    # Uncertainty (SD)
    pred_logvar <- as.numeric(final_out$logvar$cpu())
    pred_std_scaled <- sqrt(exp(pred_logvar))
    final_se <- pred_std_scaled * x_sd
  })
  
  list(imputed = final_imputed, se = final_se)
}

# ── 6. Benchmark Loop (Performance) ──────────────────────────────────────────
run_performance_benchmark <- function(n_sims = 5, n_sp = 200) {
  
  results <- data.frame()
  cat(sprintf("Running Performance Benchmark (N=%d)...\n", n_sp))
  
  for (i in 1:n_sims) {
    cat(sprintf("[Sim %d/%d] ", i, n_sims))
    
    # 1) Sim Data (Sigmoid)
    tree <- rtree(n_sp); tree$edge.length <- tree$edge.length + 1e-5
    env <- rTraitCont(tree, model = "BM", sigma = 0.5)
    env_std <- as.numeric(scale(env))
    phylo_noise <- rTraitCont(tree, model = "BM", sigma = 0.8)
    true_trait <- 10 / (1 + exp(-1 * env_std)) + phylo_noise
    
    obs_trait <- true_trait
    missing_idx <- sample(1:n_sp, size = 0.3 * n_sp)
    obs_trait[missing_idx] <- NA
    dat <- list(tree = tree, env = env_std, obs_trait = obs_trait)
    
    # 2) Rphylopars (mean + SE)
    cat("Rphylopars... ")
    phy <- tryCatch({
      rphylopars_impute_with_se(tree, obs_trait, env_std, model = "BM")
    }, error = function(e) list(mu = rep(NA, n_sp), se = rep(NA, n_sp)))
    
    phy_pred <- phy$mu[missing_idx]
    phy_se   <- phy$se[missing_idx]   # <- available if you want calibration later
    
    # 3) Phylo-DAE Calibrated
    cat("Phylo-DAE... ")
    res <- impute_calibrated(dat, epochs = 1500, lr = 0.005, corruption_rate = 0.3)
    dae_pred <- res$imputed[missing_idx]
    dae_se   <- res$se[missing_idx]   # <- also available
    
    # 4) Metrics
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

# ── 7. Calibration Test (Visualizing Uncertainty) ────────────────────────────
run_calibration_check <- function(n_sp = 300) {
  cat("\nRunning Calibration Stress Test (Heteroscedastic Noise)...\n")
  set.seed(999)
  
  tree <- rtree(n_sp); tree$edge.length <- tree$edge.length + 1e-5
  env_raw <- rTraitCont(tree, model = "BM", sigma = 0.5)
  env <- (env_raw - min(env_raw)) / (max(env_raw) - min(env_raw))
  
  # Noise grows with environment
  noise_level <- 0.1 + (1.5 * env)
  true_trait <- (5 * env) + rnorm(n_sp, mean = 0, sd = noise_level)
  
  obs_trait <- true_trait
  missing_idx <- sample(1:n_sp, size = 0.3 * n_sp)
  obs_trait[missing_idx] <- NA
  dat <- list(tree = tree, env = env, obs_trait = obs_trait)
  
  # ---- Phylo-DAE calibrated ----
  dae <- impute_calibrated(dat, epochs = 2500, lr = 0.005, corruption_rate = 0.3)
  dae_pred <- dae$imputed[missing_idx]
  dae_se   <- dae$se[missing_idx]
  
  # ---- Rphylopars mean + SE ----
  phy <- tryCatch({
    rphylopars_impute_with_se(tree, obs_trait, env, model = "BM")
  }, error = function(e) list(mu = rep(NA, n_sp), se = rep(NA, n_sp)))
  phy_pred <- phy$mu[missing_idx]
  phy_se   <- phy$se[missing_idx]
  
  truth <- true_trait[missing_idx]
  env_miss <- env[missing_idx]
  
  df_plot <- rbind(
    data.frame(Method = "Phylo-DAE",   Env = env_miss, Truth = truth, Pred = dae_pred, SE = dae_se),
    data.frame(Method = "Rphylopars",  Env = env_miss, Truth = truth, Pred = phy_pred, SE = phy_se)
  )
  
  df_plot$Lower  <- df_plot$Pred - 1.96 * df_plot$SE
  df_plot$Upper  <- df_plot$Pred + 1.96 * df_plot$SE
  df_plot$Inside <- with(df_plot, Truth >= Lower & Truth <= Upper)
  
  # Coverage per method
  cov_tab <- aggregate(Inside ~ Method, df_plot, mean)
  cov_tab$Coverage_pct <- 100 * cov_tab$Inside
  cat("\nCoverage by method (nominal 95%):\n")
  print(cov_tab[, c("Method", "Coverage_pct")], row.names = FALSE)
  
  p_cal <- ggplot(df_plot, aes(x = Env)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.20) +
    geom_point(aes(y = Truth, color = Inside), size = 1.8) +
    geom_line(aes(y = Pred), linewidth = 0.8) +
    facet_wrap(~Method, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(
      title = "Uncertainty ribbons (95% CI) vs truth",
      subtitle = "Points colored by whether truth falls inside the 95% interval"
    )
  
  p_cal
}

# ── 8. MASTER EXECUTION ──────────────────────────────────────────────────────
# A) Run Performance Benchmark (Boxplots)
perf_res <- run_performance_benchmark(n_sims = 5, n_sp = 200)

cat("\n--- PERFORMANCE SUMMARY ---\n")
print(aggregate(cbind(RMSE, Cor) ~ Method, data = perf_res, mean))

# Plot RMSE & Cor
p1 <- ggplot(perf_res, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot(alpha = 0.6) + geom_jitter(width = 0.1) +
  theme_minimal() + theme(legend.position = "none") +
  labs(title = "RMSE (Lower is Better)")

p2 <- ggplot(perf_res, aes(x = Method, y = Cor, fill = Method)) +
  geom_boxplot(alpha = 0.6) + geom_jitter(width = 0.1) +
  theme_minimal() + theme(legend.position = "none") +
  labs(title = "Correlation (Higher is Better)")

# B) Run Calibration Check (Ribbons: Phylo-DAE vs Rphylopars)
p3 <- run_calibration_check(n_sp = 300)

# C) Display All
(p1 + p2) / p3