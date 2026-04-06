###############################################################################
## Phylogenetic VGAE imputer w/ separate phylo & env encoders + MC‐dropout
###############################################################################

library(torch)
library(ape)
library(MASS)    # for mvrnorm
library(Matrix)  # for nearPD
library(Metrics) # for rmse()
library(ggplot2)

#— detect device ------------------------------------------------------------
if (cuda_is_available()) {
  device <- torch_device("cuda")
} else if (backends_mps_is_available()) {
  device <- torch_device("mps")
} else {
  device <- torch_device("cpu")
}
#cat("Using device:", device, "\n")

#— 1. Validation ------------------------------------------------------------
validate_inputs <- function(trait_data, phylo_tree, env_data, species_id) {
  if (! (is.matrix(trait_data) || is.data.frame(trait_data)))
    stop("`trait_data` must be a matrix or data.frame.")
  if (! (inherits(phylo_tree, "phylo") || is.matrix(phylo_tree)))
    stop("`phylo_tree` must be ape::phylo or an adjacency matrix.")
  if (nrow(env_data) != length(species_id) ||
      nrow(env_data) != nrow(as.matrix(trait_data)))
    stop("Rows of `env_data`, `species_id` & `trait_data` must match.")
}

#— 2. Phylogeny → adjacency -------------------------------------------------
tree_to_graph <- function(phylo_tree) {
  if (inherits(phylo_tree, "phylo")) {
    D <- cophenetic(phylo_tree)
    A <- 1/(1 + D); diag(A) <- 0
    dimnames(A) <- list(phylo_tree$tip.label, phylo_tree$tip.label)
    A
  } else phylo_tree
}
normalize_adj <- function(A) {
  nm   <- dimnames(A)
  Dinv <- diag(1/sqrt(rowSums(A))); Dinv[!is.finite(Dinv)] <- 0
  B    <- Dinv %*% A %*% Dinv
  if (!is.null(nm)) dimnames(B) <- nm
  B
}

#— 3. Prepare tensors -------------------------------------------------------
prepare_tensors <- function(trait_data, env_data, species_id, A_hat) {
  X     <- as.matrix(trait_data)
  M_obs <- !is.na(X); X[!M_obs] <- 0
  X_t   <- torch_tensor(X, dtype = torch_float(), device = device)
  M_t   <- torch_tensor(M_obs, dtype = torch_bool(), device = device)
  
  E_t   <- torch_tensor(as.matrix(env_data), dtype = torch_float(), device = device)
  
  sp_f   <- factor(species_id, levels = unique(species_id))
  # **1-based indexing** for torch
  sp_idx <- torch_tensor(as.integer(sp_f), dtype = torch_long(), device = device)
  sp_lev <- levels(sp_f)
  
  Xsp <- t(sapply(sp_lev, function(s){
    rows <- which(species_id == s)
    colMeans(X[rows,, drop = FALSE], na.rm = TRUE)
  }))
  M_sp   <- !is.na(Xsp); Xsp[!M_sp] <- 0
  X_sp_t <- torch_tensor(Xsp, dtype = torch_float(), device = device)
  M_sp_t <- torch_tensor(M_sp, dtype = torch_bool(), device = device)
  
  idx  <- match(sp_lev, rownames(A_hat))
  if (any(is.na(idx))) idx <- seq_along(sp_lev)
  A_sub <- A_hat[idx, idx, drop = FALSE]
  A_t   <- torch_tensor(A_sub, dtype = torch_float(), device = device)
  
  list(
    X      = X_t,    # n_obs × p
    M      = M_t,    # mask
    env    = E_t,    # n_obs × envdim
    sp_idx = sp_idx, # 1-based species index
    X_sp   = X_sp_t, # n_sp × p
    M_sp   = M_sp_t, # mask for X_sp
    A      = A_t     # n_sp × n_sp
  )
}

#— 4. Masked MSE -------------------------------------------------------------
masked_mse <- function(pred, target, mask) {
  mask_f <- mask$to(dtype = torch_float())
  diff2  <- (pred - target)$pow(2)
  torch_mean(diff2 * mask_f)
}

#— 5. Encoder & Decoder -----------------------------------------------------
VGAE_Encoder <- nn_module(
  initialize = function(inp, hid, lat) {
    self$gc1 <- nn_linear(inp, hid)
    self$mu  <- nn_linear(hid, lat)
    self$lv  <- nn_linear(hid, lat)
  },
  forward = function(X_sp, A) {
    H <- torch_relu(self$gc1(torch_matmul(A, X_sp)))
    list(mu = self$mu(H), logvar = self$lv(H))
  }
)

EnvEncoder <- nn_module(
  initialize = function(envdim, hid, lat) {
    self$f1 <- nn_linear(envdim, hid)
    self$f2 <- nn_linear(hid, lat)
  },
  forward = function(e) {
    h <- torch_relu(self$f1(e))
    self$f2(h)
  }
)

Decoder <- nn_module(
  initialize = function(lat, out) {
    self$fc <- nn_linear(lat, out)
  },
  forward = function(z) self$fc(z)
)

#— 6. Training + MC‐dropout inference ---------------------------------------
train_vgae_mc <- function(prep, lat, epochs, nsamp, dropout_rate = 0.3) {
  p <- prep$X$size(2); h <- p * 2
  
  enc_phy <- VGAE_Encoder(p, h, lat)$to(device = device)
  enc_env <- EnvEncoder(prep$env$size(2), h, lat)$to(device = device)
  dec     <- Decoder(lat, p)$to(device = device)
  
  opt <- optim_adam(
    c(enc_phy$parameters, enc_env$parameters, dec$parameters),
    lr = 1e-2
  )
  
  # — training loop —
  for (ep in seq_len(epochs)) {
    enc_phy$train(); enc_env$train(); dec$train()
    
    Rphy <- enc_phy(prep$X_sp, prep$A)
    Renv <- enc_env(prep$env)
    
    zsp  <- Rphy$mu +
      torch_exp(0.5 * Rphy$logvar) * torch_randn_like(Rphy$mu)
    zobs <- zsp$index_select(1, prep$sp_idx) + Renv
    
    out      <- dec(zobs)
    loss_rec <- masked_mse(out, prep$X, prep$M)
    loss_kl  <- -0.5 * torch_mean(1 + Rphy$logvar -
                                    Rphy$mu$pow(2) -
                                    torch_exp(Rphy$logvar))
    loss     <- loss_rec + loss_kl
    
    opt$zero_grad(); loss$backward(); opt$step()
    
    if (ep %% 100 == 0) {
      cat(sprintf("Epoch %4d | recon=%.4f | kl=%.4f\n",
                  ep, loss_rec$item(), loss_kl$item()))
    }
  }
  
  # — MC-dropout for uncertainty —
  enc_phy$eval(); enc_env$eval(); dec$eval()
  n_obs <- prep$X$size(1)
  pmat  <- array(NA_real_, c(nsamp, n_obs, p))
  
  for (i in seq_len(nsamp)) {
    Rphy <- enc_phy(prep$X_sp, prep$A)
    Renv <- enc_env(prep$env)
    zsp  <- Rphy$mu +
      torch_exp(0.5 * Rphy$logvar) * torch_randn_like(Rphy$mu)
    zobs <- zsp$index_select(1, prep$sp_idx) + Renv
    
    # manual MC-dropout on zobs:
    # mask  <- torch_rand(zobs$size(), device = device) > dropout_rate
    # mask_f <- mask$to(dtype = torch_float())
    # zdrop <- (mask_f * zobs) / (1 - dropout_rate)
    # 
    # pmat[i,,] <- as_array(dec(zdrop))
    
    pmat[i,,] <- as_array(dec(zobs))
  }
  
  list(
    mean = apply(pmat, 2:3, mean),
    sd   = apply(pmat, 2:3, sd)
  )
}

#— 7. High-level wrapper ----------------------------------------------------
#' Impute missing with phylo+env VGAE + MC-dropout

impute_phylo <- function(trait_data,
                         phylo_tree,
                         env_data,
                         species_id,
                         latent_dim           = 16,
                         epochs               = 300,
                         n_samples            = 50,
                         dropout_rate         = 0.3,
                         mask_obs_uncertainty = TRUE,
                         restore_observed     = TRUE) {
  ## 0) Sanity checks
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  
  ## 1) Build phylo adjacency and prepare tensors
  A_hat <- normalize_adj(tree_to_graph(phylo_tree))
  prep  <- prepare_tensors(trait_data, env_data, species_id, A_hat)
  
  ## 2) Train & MC‐dropout inference
  fit    <- train_vgae_mc(prep, latent_dim, epochs, n_samples, dropout_rate)
  mean_m <- fit$mean   # n_obs × p
  sd_m   <- fit$sd     # n_obs × p
  
  ## 3) Mask or restore observed entries
  orig_mat <- as.matrix(trait_data)
  obs_mask <- !is.na(orig_mat)
  if (mask_obs_uncertainty) sd_m  [obs_mask] <- 0
  if (restore_observed)    mean_m[obs_mask] <- orig_mat[obs_mask]
  
  ## 4) Build outputs
  original_data   <- as.data.frame(orig_mat)
  completed_data  <- as.data.frame(mean_m)
  uncertainty     <- as.data.frame(sd_m)
  colnames(completed_data) <- colnames(trait_data)
  colnames(uncertainty)    <- colnames(trait_data)
  
  ## 5) Long‐format only‐missing table
  miss_idx <- which(is.na(orig_mat), arr.ind = TRUE)
  imputed_missing <- data.frame(
    obs_id      = miss_idx[,1],
    species_id  = species_id[ miss_idx[,1] ],
    trait       = colnames(original_data)[ miss_idx[,2] ],
    imputed     = completed_data[ miss_idx ],
    uncertainty =     uncertainty[ miss_idx ],
    stringsAsFactors = FALSE
  )
  
  ## 6) Return everything, including inputs for reference
  list(
    original_data   = original_data,
    completed_data  = completed_data,
    uncertainty     = uncertainty,
    imputed_missing = imputed_missing,
    env_data        = env_data,
    species_id      = species_id,
    phylo_tree      = phylo_tree
  )
}
###############################################################################
## 8. Simulation & Evaluation
###############################################################################
set.seed(42)

# A) phylogeny + phylo-correlation
n_sp <- 500
tree <- rtree(n_sp)
V_phy <- cov2cor(vcv(tree))

# B) correlated environmental covariates
p_env    <- 10
R_env    <- matrix(runif(p_env^2, .05, .8), p_env, p_env)
R_env    <- (R_env + t(R_env))/2; diag(R_env) <- 1
Sigma_env <- as.matrix(nearPD(R_env)$mat)
env_data  <- mvrnorm(n_sp, rep(0,p_env), Sigma_env)
colnames(env_data) <- paste0("env",1:p_env)

# C) phylogenetic species-level traits
p_tr <- 10
Z    <- sapply(1:p_tr, function(i) mvrnorm(1, rep(0,n_sp), V_phy))
R_tr <- matrix(runif(p_tr^2, .05, .8), p_tr, p_tr)
R_tr <- (R_tr + t(R_tr))/2; diag(R_tr) <- 1
Sigma_tr <- as.matrix(nearPD(R_tr)$mat)
Xsp      <- Z %*% chol(Sigma_tr)

# D) random env→trait slopes
B <- matrix(runif(p_tr * p_env, -1,1), p_tr, p_env)

# E) build obs-level data + noise
obs_traits <- Xsp + env_data %*% t(B) +
  matrix(rnorm(n_sp*p_tr,0, .1), n_sp, p_tr)
colnames(obs_traits) <- paste0("trait",1:p_tr)

# F) mask 30% at random
mask_mat        <- matrix(runif(n_sp*p_tr), n_sp, p_tr) < .30
obs_traits_mask <- obs_traits
obs_traits_mask[mask_mat] <- NA

# G) impute & evaluate
res <- impute_phylo(
  trait_data   = as.data.frame(obs_traits_mask),
  phylo_tree   = tree,
  env_data     = as.data.frame(env_data),
  species_id   = tree$tip.label,
  latent_dim   = 32,
  epochs       = 2000,
  n_samples    = 100,
  dropout_rate = 0
)

# Full completed matrix
head(res$completed_data)

# Just the imputed entries
head(res$imputed_missing)

completed <- as.matrix(res$completed_data)
truth     <- obs_traits
rmse_val  <- rmse(truth[mask_mat], completed[mask_mat])
cat("RMSE on masked entries:", round(rmse_val,3), "\n")

# scatterplot
df_mask <- data.frame(
  truth   = truth[mask_mat],
  imputed = completed[mask_mat]
)
ggplot(df_mask, aes(x=truth, y=imputed)) +
  geom_point(alpha=0.4) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  coord_equal() + theme_minimal() +
  labs(title="Imputed vs True (masked entries)",
       x="True", y="Imputed")