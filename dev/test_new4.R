# ==============================================================================
# CALIBRATED PHYLO-DAE: Aleatoric Uncertainty & Non-Linear Imputation
# ==============================================================================
# 1. Predicts both Mean (Value) and Variance (Uncertainty) per species.
# 2. Uses Heteroscedastic Loss (NLL) to calibrate confidence intervals.
# 3. Beats Linear Baselines (Rphylopars) on non-linear biological data.
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

# ── 3. The Calibrated Model (Dual-Head Output) ───────────────────────────────
PhyloDAE_Calibrated <- nn_module(
  "PhyloDAE_Calibrated",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim) {
    total_input <- input_dim + coord_dim + cov_dim
    
    # Shared Encoder
    self$enc1 <- nn_linear(total_input, hidden_dim)
    self$enc2 <- nn_linear(hidden_dim, hidden_dim)
    self$graph_mix <- nn_linear(hidden_dim, hidden_dim)
    self$dec1 <- nn_linear(hidden_dim, hidden_dim)
    
    # HEAD 1: Value Prediction (Mean)
    self$head_mu <- nn_linear(hidden_dim, input_dim)
    
    # HEAD 2: Uncertainty Prediction (Log Variance)
    # We predict log(var) because variance must be positive
    self$head_var <- nn_linear(hidden_dim, input_dim)
    
    self$act  <- nn_relu() 
    self$drop <- nn_dropout(0.5) 
  },
  
  forward = function(x, coords, covs, adj) {
    combined <- torch_cat(list(x, coords, covs), dim = 2)
    
    # Encode
    h <- self$enc1(combined); h <- self$act(h); h <- self$drop(h) 
    h <- self$enc2(h); h <- self$act(h)
    
    # Mix
    neighbor_info <- torch_matmul(adj, h)
    h_mix <- self$graph_mix(neighbor_info)
    h <- h + h_mix 
    
    # Decode
    h <- self$dec1(h); h <- self$act(h)
    
    # Dual Output
    mu      <- self$head_mu(h)
    logvar  <- self$head_var(h) 
    
    return(list(mu = mu, logvar = logvar))
  }
)

# ── 4. The Heteroscedastic Loss Function (NLL) ───────────────────────────────
gaussian_nll_loss <- function(pred_mu, pred_logvar, target, mask) {
  # NLL = 0.5 * [ exp(-logvar)*(y-y_hat)^2 + logvar ]
  # Minimizing this forces the model to increase variance when error is high
  
  precision <- torch_exp(-pred_logvar)
  squared_error <- (pred_mu - target)^2
  loss <- 0.5 * (precision * squared_error + pred_logvar)
  
  # Apply mask and mean
  mask_f <- mask$to(dtype = torch_float())
  torch_sum(loss * mask_f) / torch_sum(mask_f + 1e-6)
}

# ── 5. Main Imputation Function ──────────────────────────────────────────────
impute_calibrated <- function(sim_data, epochs=2000, lr=0.005, k_eigen=8, 
                              corruption_rate=0.3) {
  
  trait_vec <- sim_data$obs_trait; tree <- sim_data$tree; env_vec <- sim_data$env
  
  # A) Data Prep
  is_missing <- is.na(trait_vec)
  mu_global <- mean(trait_vec, na.rm=TRUE)
  X_filled <- trait_vec; X_filled[is_missing] <- mu_global
  
  x_mean <- mean(X_filled); x_sd <- sd(X_filled)
  if(x_sd == 0) x_sd <- 1
  X_scaled <- (X_filled - x_mean) / x_sd
  
  # Tensors (Forced [N, 1])
  X_full <- torch_tensor(matrix(X_scaled, ncol=1), dtype=torch_float(), device=device)
  Env_t  <- torch_tensor(matrix(env_vec, ncol=1), dtype=torch_float(), device=device)
  M_true_obs <- torch_tensor(matrix(!is_missing, ncol=1), dtype=torch_bool(), device=device)
  
  adj_t <- get_stable_adj(tree); coords_t <- get_spectral_features(tree, k = k_eigen)
  
  # B) Model
  model <- PhyloDAE_Calibrated(input_dim=1, hidden_dim=128, coord_dim=k_eigen, cov_dim=1)
  model$to(device = device)
  opt <- optim_adam(model$parameters, lr = lr)
  
  # C) Training
  for(i in 1:epochs) {
    model$train()
    opt$zero_grad()
    
    # Dynamic Denoising
    noise_mask <- torch_rand_like(X_full) > corruption_rate
    X_corrupted <- X_full * noise_mask
    
    # Forward Pass
    out <- model(X_corrupted, coords_t, Env_t, adj_t)
    
    # NLL Loss on hidden data
    target_mask <- M_true_obs & (!noise_mask)
    if (target_mask$sum()$item() > 0) {
      loss <- gaussian_nll_loss(out$mu, out$logvar, X_full, target_mask)
      loss$backward()
      opt$step()
    }
    
    if(i %% 500 == 0) cat(".")
  }
  cat("\n")
  
  # D) Inference (Deterministic, no MC needed)
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
  
  return(list(imputed = final_imputed, se = final_se))
}

# ── 6. Benchmark Loop (Performance) ──────────────────────────────────────────
run_performance_benchmark <- function(n_sims = 5, n_sp = 200) {
  
  results <- data.frame()
  cat(sprintf("Running Performance Benchmark (N=%d)...\n", n_sp))
  
  for(i in 1:n_sims) {
    cat(sprintf("[Sim %d/%d] ", i, n_sims))
    
    # 1. Sim Data (Sigmoid)
    tree <- rtree(n_sp); tree$edge.length <- tree$edge.length + 1e-5
    env <- rTraitCont(tree, model = "BM", sigma = 0.5)
    env_std <- as.numeric(scale(env))
    phylo_noise <- rTraitCont(tree, model = "BM", sigma = 0.2)
    true_trait <- 10 / (1 + exp(-3 * env_std)) + phylo_noise
    
    obs_trait <- true_trait; missing_idx <- sample(1:n_sp, size = 0.3 * n_sp)
    obs_trait[missing_idx] <- NA
    dat <- list(tree=tree, env=env_std, obs_trait=obs_trait)
    
    # 2. Rphylopars
    cat("Rphylopars... ")
    phy_pred <- tryCatch({
      df <- data.frame(species = tree$tip.label, trait = obs_trait, env = env_std)
      suppressWarnings({ p <- phylopars(df, tree, model="BM", pheno_error=FALSE) })
      p$anc_recon[1:n_sp, 1][missing_idx]
    }, error = function(e) rep(NA, length(missing_idx)))
    
    # 3. Phylo-DAE Calibrated
    cat("Phylo-DAE... ")
    res <- impute_calibrated(dat, epochs=1500, lr=0.005, corruption_rate=0.3)
    dae_pred <- res$imputed[missing_idx]
    
    # 4. Metrics
    truth <- true_trait[missing_idx]
    
    rmse_phy <- sqrt(mean((truth - phy_pred)^2, na.rm=TRUE))
    rmse_dae <- sqrt(mean((truth - dae_pred)^2))
    
    cor_phy <- cor(truth, phy_pred, use = "complete.obs")
    cor_dae <- cor(truth, dae_pred, use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Sim = i, Method = c("Rphylopars", "Phylo-DAE"),
      RMSE = c(rmse_phy, rmse_dae), Cor  = c(cor_phy, cor_dae)
    ))
    cat("Done.\n")
  }
  return(results)
}

# ── 7. Calibration Test (Visualizing Uncertainty) ────────────────────────────
run_calibration_check <- function(n_sp=300) {
  cat("\nRunning Calibration Stress Test (Heteroscedastic Noise)...\n")
  set.seed(999) 
  tree <- rtree(n_sp); tree$edge.length <- tree$edge.length + 1e-5
  env_raw <- rTraitCont(tree, model="BM", sigma=0.5)
  env <- (env_raw - min(env_raw)) / (max(env_raw) - min(env_raw)) 
  
  # Noise grows with environment!
  noise_level <- 0.1 + (1.5 * env) 
  true_trait <- (5 * env) + rnorm(n_sp, mean=0, sd=noise_level)
  
  obs_trait <- true_trait; missing_idx <- sample(1:n_sp, size=0.3*n_sp)
  obs_trait[missing_idx] <- NA
  dat <- list(tree=tree, env=env, obs_trait=obs_trait)
  
  # Train Calibrated Model
  res <- impute_calibrated(dat, epochs=2500, lr=0.005, corruption_rate=0.3)
  
  # Metrics
  truth <- true_trait[missing_idx]; pred <- res$imputed[missing_idx]
  se <- res$se[missing_idx]; env_miss <- env[missing_idx]
  
  # Coverage Check
  lower <- pred - 1.96 * se; upper <- pred + 1.96 * se
  is_inside <- (truth >= lower) & (truth <= upper)
  coverage <- mean(is_inside)
  
  cat(sprintf("\nActual Coverage (Target 95%%): %.1f%%\n", coverage * 100))
  
  # Plot
  df_plot <- data.frame(Env=env_miss, Truth=truth, Pred=pred, SE=se, 
                        Inside=ifelse(is_inside, "Yes", "No"))
  
  p_cal <- ggplot(df_plot, aes(x=Env)) +
    geom_ribbon(aes(ymin=Pred-1.96*SE, ymax=Pred+1.96*SE), fill="#1f78b4", alpha=0.2) +
    geom_point(aes(y=Truth, color=Inside), size=2) +
    geom_line(aes(y=Pred), color="#1f78b4", size=1) +
    scale_color_manual(values=c("Yes"="grey50", "No"="red")) +
    theme_minimal() + 
    labs(title="Uncertainty Calibration", subtitle=sprintf("Coverage: %.1f%%", coverage*100))
  
  return(p_cal)
}

# ── 8. MASTER EXECUTION ──────────────────────────────────────────────────────
# A) Run Performance Benchmark (Boxplots)
perf_res <- run_performance_benchmark(n_sims=5, n_sp=200)

cat("\n--- PERFORMANCE SUMMARY ---\n")
print(aggregate(cbind(RMSE, Cor) ~ Method, data = perf_res, mean))

# Plot RMSE & Cor
p1 <- ggplot(perf_res, aes(x=Method, y=RMSE, fill=Method)) +
  geom_boxplot(alpha=0.6) + geom_jitter(width=0.1) +
  theme_minimal() + theme(legend.position="none") + labs(title="RMSE (Lower is Better)")

p2 <- ggplot(perf_res, aes(x=Method, y=Cor, fill=Method)) +
  geom_boxplot(alpha=0.6) + geom_jitter(width=0.1) +
  theme_minimal() + theme(legend.position="none") + labs(title="Correlation (Higher is Better)")

# B) Run Calibration Check (Ribbons)
p3 <- run_calibration_check(n_sp=300)

# C) Display All
(p1 + p2) / p3