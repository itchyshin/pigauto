# compare to cont1 - I am trying to implement more complex DL models here 

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

## ── 0. small utilities ────────────────────────────────────────────────────
validate_inputs <- function(trait_data, phylo_tree, env_data, species_id) {
  if (! (is.matrix(trait_data) || is.data.frame(trait_data)))
    stop("`trait_data` must be a matrix or data.frame.")
  if (! (inherits(phylo_tree, "phylo") || is.matrix(phylo_tree)))
    stop("`phylo_tree` must be ape::phylo or an adjacency matrix.")
  if (nrow(env_data) != length(species_id) ||
      nrow(env_data) != nrow(as.matrix(trait_data)))
    stop("Rows of `env_data`, `species_id`, `trait_data` must match.")
}

tree_to_graph <- function(phylo_tree) {
  if (inherits(phylo_tree, "phylo")) {
    D <- cophenetic(phylo_tree)
    A <- 1/(1 + D); diag(A) <- 0
    dimnames(A) <- list(phylo_tree$tip.label, phylo_tree$tip.label)
    A
  } else phylo_tree
}
normalize_adj <- function(A) {
  nm <- dimnames(A)
  Dinv <- diag(1/sqrt(rowSums(A))); Dinv[!is.finite(Dinv)] <- 0
  B <- Dinv %*% A %*% Dinv
  if (!is.null(nm)) dimnames(B) <- nm
  B
}

make_mlp <- function(in_dim, hidden_dim, out_dim,
                     n_layers, dropout = 0.0) {
  layers <- list()
  if (n_layers == 1) {
    layers[[1]] <- nn_linear(in_dim, out_dim)
  } else {
    layers <- append(layers, list(nn_linear(in_dim, hidden_dim), nn_relu()))
    if (dropout > 0) layers <- append(layers, list(nn_dropout(p = dropout)))
    if (n_layers > 2) {
      for (i in 2:(n_layers-1)) {
        layers <- append(layers, list(nn_linear(hidden_dim, hidden_dim), nn_relu()))
        if (dropout > 0) layers <- append(layers, list(nn_dropout(p = dropout)))
      }
    }
    layers <- append(layers, list(nn_linear(hidden_dim, out_dim)))
  }
  do.call(nn_sequential, layers)
}

## ── 1. tensors ────────────────────────────────────────────────────────────
prepare_tensors <- function(trait_data, env_data, species_id, A_hat) {
  # 1) raw numeric matrix, impute NAs with 0
  X     <- as.matrix(trait_data)
  M_obs <- !is.na(X)
  X[!M_obs] <- 0
  
  # 2) turn *everything* into FLOAT tensors
  X_t   <- torch_tensor(X,    dtype = torch_float(), device = device)
  M_t   <- torch_tensor(M_obs * 1, dtype = torch_float(), device = device)
  E_t   <- torch_tensor(as.matrix(env_data),
                        dtype = torch_float(), device = device)
  
  # 3) species indexing
  sp_f   <- factor(species_id, levels = unique(species_id))
  sp_idx <- torch_tensor(as.integer(sp_f),
                         dtype = torch_long(), device = device)
  
  # 4) species‐level averages
  sp_lev <- levels(sp_f)
  Xsp    <- t(sapply(sp_lev, function(s) {
    rows <- which(species_id == s)
    colMeans(X[rows, , drop = FALSE], na.rm = TRUE)
  }))
  M_sp   <- !is.na(Xsp)
  Xsp[!M_sp] <- 0
  X_sp_t <- torch_tensor(Xsp, dtype = torch_float(), device = device)
  M_sp_t <- torch_tensor(M_sp * 1, dtype = torch_float(), device = device)
  
  # 5) adjacency
  idx   <- match(sp_lev, rownames(A_hat))
  if (any(is.na(idx))) idx <- seq_along(sp_lev)
  A_sub <- A_hat[idx, idx, drop = FALSE]
  A_t   <- torch_tensor(A_sub, dtype = torch_float(), device = device)
  
  list(
    X      = X_t,
    M      = M_t,
    env    = E_t,
    sp_idx = sp_idx,
    X_sp   = X_sp_t,
    M_sp   = M_sp_t,
    A      = A_t
  )
}
#— 4. Masked MSE (stable) ----------------------------------------------------
masked_mse <- function(pred, target, mask) {
  diffsq <- (pred - target)$pow(2)
  torch_sum(diffsq * mask) / torch_sum(mask)
}

## ── 2. modules ────────────────────────────────────────────────────────────
VGAE_Encoder <- nn_module(
  initialize = function(inp, hid, lat, n_layers, dropout) {
    self$gcn <- make_mlp(inp, hid, hid, n_layers, dropout)
    self$mu  <- nn_linear(hid, lat)
    self$lv  <- nn_linear(hid, lat)
  },
  forward = function(X_sp, A) {
    H <- self$gcn(torch_matmul(A, X_sp))
    list(mu = self$mu(H), logvar = self$lv(H))
  }
)

EnvEncoder <- nn_module(
  initialize = function(envdim, hid, lat, n_layers, dropout) {
    self$mlp <- make_mlp(envdim, hid, lat, n_layers, dropout)
  },
  forward = function(e) self$mlp(e)
)

Decoder <- nn_module(
  initialize = function(lat, hid, out, n_layers, dropout) {
    self$mlp <- make_mlp(lat, hid, out, n_layers, dropout)
  },
  forward = function(z) self$mlp(z)
)

## ── 3. training + MC‑dropout ─────────────────────────────────────────────
train_vgae_mc <- function(prep, latent_dim,
                          hidden_dim    = NULL,
                          gcn_layers    = 2,
                          env_layers    = 2,
                          dec_layers    = 2,
                          enc_dropout   = 0.1,
                          dec_dropout   = 0.1,
                          epochs        = 300,
                          nsamp         = 50,
                          mc_drop_rate  = 0.3) {
  
  p <- prep$X$size(2)
  if (is.null(hidden_dim)) hidden_dim <- 2 * p
  
  enc_phy <- VGAE_Encoder(p, hidden_dim, latent_dim,
                          gcn_layers, enc_dropout)$to(device)
  enc_env <- EnvEncoder(prep$env$size(2), hidden_dim, latent_dim,
                        env_layers, enc_dropout)$to(device)
  dec     <- Decoder(latent_dim, hidden_dim, p,
                     dec_layers, dec_dropout)$to(device)
  
  opt <- optim_adam(c(enc_phy$parameters, enc_env$parameters, dec$parameters),
                    lr = 1e-3)
  
  loss_hist <- numeric(epochs)
  
  for (ep in seq_len(epochs)) {
    enc_phy$train(); enc_env$train(); dec$train()
    
    r   <- enc_phy(prep$X_sp, prep$A)
    zsp <- r$mu + torch_exp(0.5 * r$logvar) * torch_randn_like(r$mu)
    z   <- zsp$index_select(1, prep$sp_idx) + enc_env(prep$env)
    
    out  <- dec(z)
    lrec <- masked_mse(out, prep$X, prep$M)
    lkl  <- -0.5 * torch_mean(1 + r$logvar - r$mu$pow(2) - torch_exp(r$logvar))
    loss <- lrec + lkl
    
    opt$zero_grad(); loss$backward(); opt$step()
    loss_hist[ep] <- as.numeric(loss)
    
    if (ep %% 100 == 0)
      cat(sprintf("Epoch %4d | recon=%.4f | kl=%.4f\n",
                  ep, as.numeric(lrec), as.numeric(lkl)))
  }
  
  ## ── MC‑dropout draws ─────────────────────────────────────────────
  enc_phy$eval(); enc_env$eval(); dec$eval()
  n_obs <- prep$X$size(1)
  pmat  <- array(NA_real_, c(nsamp, n_obs, p))
  
  for (i in seq_len(nsamp)) {
    r   <- enc_phy(prep$X_sp, prep$A)
    zsp <- r$mu + torch_exp(0.5 * r$logvar) * torch_randn_like(r$mu)
    z   <- zsp$index_select(1, prep$sp_idx) + enc_env(prep$env)
    
    zdrop <- nnf_dropout(z, p = mc_drop_rate, training = TRUE)
    pmat[i,,] <- as_array(dec(zdrop))
  }
  
  list(mean = apply(pmat, 2:3, mean),
       sd   = apply(pmat, 2:3,  sd),
       loss_history = loss_hist)
}

## ── 4. user wrapper ───────────────────────────────────────────────────────
impute_phylo <- function(trait_data, phylo_tree, env_data, species_id,
                         latent_dim   = 32,
                         hidden_dim   = NULL,
                         gcn_layers   = 2,
                         env_layers   = 2,
                         dec_layers   = 2,
                         enc_dropout  = 0.1,
                         dec_dropout  = 0.1,
                         epochs       = 300,
                         n_samples    = 50,
                         mc_drop_rate = 0.3,
                         mask_obs_uncertainty = TRUE,
                         restore_observed     = TRUE) {
  
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  A_hat <- normalize_adj(tree_to_graph(phylo_tree))
  prep  <- prepare_tensors(trait_data, env_data, species_id, A_hat)
  
  fit  <- train_vgae_mc(prep, latent_dim, hidden_dim,
                        gcn_layers, env_layers, dec_layers,
                        enc_dropout, dec_dropout,
                        epochs, n_samples, mc_drop_rate)
  
  mean_m <- fit$mean
  sd_m   <- fit$sd
  orig   <- as.matrix(trait_data)
  obs_m  <- !is.na(orig)
  if (mask_obs_uncertainty) sd_m[obs_m] <- 0
  if (restore_observed)     mean_m[obs_m] <- orig[obs_m]
  
  completed  <- as.data.frame(mean_m)
  uncertainty<- as.data.frame(sd_m)
  colnames(completed)   <- colnames(trait_data)
  colnames(uncertainty) <- colnames(trait_data)
  
  miss_idx <- which(is.na(orig), arr.ind = TRUE)
  imputed_missing <- data.frame(
    obs_id      = miss_idx[,1],
    species_id  = species_id[miss_idx[,1]],
    trait       = colnames(trait_data)[miss_idx[,2]],
    imputed     = completed[miss_idx],
    uncertainty = uncertainty[miss_idx],
    stringsAsFactors = FALSE
  )
  
  list(
    original_data   = as.data.frame(trait_data),
    completed_data  = completed,
    uncertainty     = uncertainty,
    imputed_missing = imputed_missing,
    loss_history    = fit$loss_history,
    env_data        = env_data,
    species_id      = species_id,
    phylo_tree      = phylo_tree
  )
}
#— 5. example usage ───────────────────────────────────────────────────────
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
  hidden_dim   = 128,   # width
  gcn_layers   = 3,     # depth of phylo encoder
  env_layers   = 3,     # depth of env encoder
  dec_layers   = 2,     # depth of decoder
  epochs       = 600,
  n_samples    = 100,
  mc_drop_rate = 0.2
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
