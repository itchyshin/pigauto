###############################################################################
## pigauto: Phylogenetic GAE imputer w/ MC‐dropout uncertainty + scaling
###############################################################################

library(torch)
library(ape)
library(MASS)    # for mvrnorm
library(Matrix)  # for nearPD
library(Metrics) # for rmse()
library(ggplot2)

#— detect device ------------------------------------------------------------
device <- if (cuda_is_available()) {
  torch_device("cuda")
} else if (backends_mps_is_available()) {
  torch_device("mps")
} else {
  torch_device("cpu")
}

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
  Dinv <- diag(1/sqrt(rowSums(A))); Dinv[!is.finite(Dinv)] <- 0
  B    <- Dinv %*% A %*% Dinv
  dimnames(B) <- dimnames(A)
  B
}

#— 3. Prepare tensors -------------------------------------------------------
prepare_tensors <- function(X_scaled, env_data, species_id, A_hat) {
  # X_scaled: standardized trait matrix (N × p)
  M_obs <- !is.na(X_scaled)
  X_scaled[!M_obs] <- 0
  
  X_t   <- torch_tensor(X_scaled, dtype = torch_float(), device = device)
  M_t   <- torch_tensor(M_obs,   dtype = torch_bool(),  device = device)
  
  E_t   <- torch_tensor(as.matrix(env_data),
                        dtype = torch_float(), device = device)
  
  sp_f   <- factor(species_id, levels = unique(species_id))
  sp_idx <- torch_tensor(as.integer(sp_f),
                         dtype = torch_long(), device = device)
  sp_lev <- levels(sp_f)
  
  Xsp <- t(sapply(sp_lev, function(s){
    rows <- which(species_id == s)
    colMeans(X_scaled[rows,, drop = FALSE], na.rm = TRUE)
  }))
  M_sp   <- !is.na(Xsp); Xsp[!M_sp] <- 0
  X_sp_t <- torch_tensor(Xsp, dtype = torch_float(), device = device)
  M_sp_t <- torch_tensor(M_sp, dtype = torch_bool(),  device = device)
  
  idx  <- match(sp_lev, rownames(A_hat))
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

#— 4. Masked MSE -------------------------------------------------------------
masked_mse <- function(pred, target, mask) {
  mask_f <- mask$to(dtype = torch_float())
  diff2  <- (pred - target)$pow(2)
  torch_mean(diff2 * mask_f)
}

#— 5. GAE Encoder & Decoder w/ dropout --------------------------------------
GAE_Encoder <- nn_module(
  initialize = function(inp, hid, lat, dropout = 0.3) {
    self$gc1  <- nn_linear(inp, hid)
    self$drop <- nn_dropout(p = dropout)
    self$fc   <- nn_linear(hid, lat)
  },
  forward = function(X_sp, A) {
    h <- torch_matmul(A, X_sp)
    h <- torch_relu(self$gc1(h))
    h <- self$drop(h)
    self$fc(h)
  }
)

Env_Encoder <- nn_module(
  initialize = function(envdim, hid, lat, dropout = 0.3) {
    self$f1   <- nn_linear(envdim, hid)
    self$drop <- nn_dropout(p = dropout)
    self$f2   <- nn_linear(hid, lat)
  },
  forward = function(e) {
    h <- torch_relu(self$f1(e))
    h <- self$drop(h)
    self$f2(h)
  }
)

GAE_Decoder <- nn_module(
  initialize = function(lat, out) {
    self$fc <- nn_linear(lat, out)
  },
  forward = function(z) self$fc(z)
)

#— 6. Training loop (no KL term) --------------------------------------------
train_gae <- function(prep, latent_dim,
                      hidden_mult  = 2,
                      epochs       = 300,
                      lr           = 1e-2,
                      dropout_rate = 0.3) {
  
  ## make sure the layer sizes are INTEGERS -------------------------------
  p <- as.integer(prep$X$size(2))           # input dimension
  h <- as.integer(round(p * hidden_mult))   # hidden dimension
  
  enc_phy <- GAE_Encoder(p, h, latent_dim, dropout_rate)$to(device = device)
  enc_env <- Env_Encoder(prep$env$size(2), h, latent_dim, dropout_rate)$to(device = device)
  dec     <- GAE_Decoder(latent_dim, p)$to(device = device)
  
  opt <- optim_adam(
    c(enc_phy$parameters, enc_env$parameters, dec$parameters),
    lr = lr,
    weight_decay = 1e-4
  )
  
  for (ep in seq_len(epochs)) {
    enc_phy$train(); enc_env$train(); dec$train()
    
    Zphy <- enc_phy(prep$X_sp, prep$A)
    Zenv <- enc_env(prep$env)
    Zobs <- Zphy$index_select(1, prep$sp_idx) + Zenv
    
    out      <- dec(Zobs)
    loss_rec <- masked_mse(out, prep$X, prep$M)
    
    opt$zero_grad(); loss_rec$backward(); opt$step()
    
    if (ep %% 100 == 0)
      cat(sprintf("Epoch %4d | recon=%.4f\n", ep, loss_rec$item()))
  }
  list(enc_phy = enc_phy, enc_env = enc_env, dec = dec)
}

#— 7. Wrapper + MC-dropout inference + back-transform -----------------------
impute_phylo <- function(trait_data,
                         phylo_tree,
                         env_data,
                         species_id,
                         latent_dim           = 16,
                         epochs               = 300,
                         n_samples            = 50,
                         hidden_mult          = 2,
                         lr                   = 1e-2,
                         dropout_rate         = 0.3,
                         mask_obs_uncertainty = TRUE,
                         restore_observed     = TRUE) {
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  
  # — Standardize each trait column —
  X_mat <- as.matrix(trait_data)
  mu    <- colMeans(X_mat, na.rm = TRUE)
  sigma <- apply(X_mat, 2, stats::sd, na.rm = TRUE)
  X_scaled <- sweep(sweep(X_mat, 2, mu, "-"), 2, sigma, "/")
  
  A_hat <- normalize_adj(tree_to_graph(phylo_tree))
  prep  <- prepare_tensors(X_scaled, env_data, species_id, A_hat)
  
  # — Train with dropout —
  models <- train_gae(prep, latent_dim, hidden_mult,
                      epochs, lr, dropout_rate)
  
  # — MC-dropout sampling (keep dropout on) —
  n_obs <- prep$X$size(1)
  p     <- prep$X$size(2)
  pmat  <- array(NA_real_, c(n_samples, n_obs, p))
  
  for (i in seq_len(n_samples)) {
    Zphy <- models$enc_phy(prep$X_sp, prep$A)
    Zenv <- models$enc_env(prep$env)
    Zobs <- Zphy$index_select(1, prep$sp_idx) + Zenv
    pmat[i,,] <- as_array(models$dec(Zobs))
  }
  
  # — back-transform to original scale —
  mean_std <- apply(pmat, 2:3, mean)
  sd_std   <- apply(pmat, 2:3, sd)
  mean_m   <- sweep(sweep(mean_std, 2, sigma, "*"), 2, mu, "+")
  sd_m     <- sweep(sd_std,   2, sigma, "*")
  
  # — mask/restore observed entries —
  orig_mat <- as.matrix(trait_data)
  obs_mask <- !is.na(orig_mat)
  if (mask_obs_uncertainty) sd_m  [obs_mask] <- 0
  if (restore_observed)    mean_m[obs_mask] <- orig_mat[obs_mask]
  
  # — assemble outputs —
  completed_data  <- as.data.frame(mean_m)
  uncertainty     <- as.data.frame(sd_m)
  colnames(completed_data) <- colnames(trait_data)
  colnames(uncertainty)    <- colnames(trait_data)
  
  miss_idx <- which(is.na(orig_mat), arr.ind = TRUE)
  imputed_missing <- data.frame(
    obs_id      = miss_idx[,1],
    species_id  = species_id[ miss_idx[,1] ],
    trait       = colnames(trait_data)[ miss_idx[,2] ],
    imputed     = completed_data[ miss_idx ],
    uncertainty =     uncertainty[ miss_idx ],
    stringsAsFactors = FALSE
  )
  
  list(
    original_data   = as.data.frame(orig_mat),
    completed_data  = completed_data,
    uncertainty     = uncertainty,
    imputed_missing = imputed_missing,
    env_data        = env_data,
    species_id      = species_id,
    phylo_tree      = phylo_tree
  )
}

#— 8. Quick test -------------------------------------------------------------
set.seed(42)
n_sp <- 500
tree <- rtree(n_sp)
V_phy <- cov2cor(vcv(tree))

p_env <- 10
R_env <- (matrix(runif(p_env^2, .05, .8), p_env, p_env) +
            t(matrix(runif(p_env^2, .05, .8), p_env, p_env))) / 2
diag(R_env) <- 1
Sigma_env <- as.matrix(nearPD(R_env)$mat)
env_data  <- as.data.frame(mvrnorm(n_sp, rep(0,p_env), Sigma_env))

p_tr <- 10
Z    <- sapply(1:p_tr, function(i) mvrnorm(1, rep(0, n_sp), V_phy))
Sigma_tr <- as.matrix(nearPD((matrix(runif(p_tr^2, .05, .8), p_tr, p_tr) +
                                t(matrix(runif(p_tr^2, .05, .8), p_tr, p_tr))) / 2)$mat)
Xsp      <- Z %*% chol(Sigma_tr)

B <- matrix(runif(p_tr * p_env, -1, 1), p_tr, p_env)
obs_traits <- Xsp + as.matrix(env_data) %*% t(B) +
  matrix(rnorm(n_sp * p_tr, 0, .1), n_sp, p_tr)
colnames(obs_traits) <- paste0("trait", 1:p_tr)

mask_mat <- matrix(runif(n_sp * p_tr), n_sp, p_tr) < 0.30
obs_traits_mask <- obs_traits
obs_traits_mask[mask_mat] <- NA

res <- impute_phylo(
  trait_data   = as.data.frame(obs_traits_mask),
  phylo_tree   = tree,
  env_data     = env_data,
  species_id   = tree$tip.label,
  latent_dim   = 32,
  epochs       = 2000,
  n_samples    = 100,
  dropout_rate = 0
)

completed <- as.matrix(res$completed_data)
truth     <- obs_traits
df_mask   <- data.frame(truth = truth[mask_mat],
                        imputed = completed[mask_mat])


rmse_val  <- rmse(truth[mask_mat], completed[mask_mat])
cat("RMSE on masked entries:", round(rmse_val,3), "\n")


ggplot(df_mask, aes(truth, imputed)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() + theme_minimal() +
  labs(title = "Imputed vs True (masked entries)")