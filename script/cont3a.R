###############################################################################
## Phylogenetic VGAE imputer w/ per‐trait standardisation + MC‐dropout
##   + Early‐stopping + Model checkpointing + History
##   + Configurable number of hidden layers
###############################################################################

library(torch)
library(ape)
library(MASS)     # mvrnorm
library(Matrix)   # nearPD
library(Metrics)  # rmse
library(ggplot2)

#— 0) Detect device ----------------------------------------------------------
device <- if (cuda_is_available()) {
  torch_device("cuda")
} else if (backends_mps_is_available()) {
  torch_device("mps")
} else {
  torch_device("cpu")
}

#— 1) Sanity checks ----------------------------------------------------------
validate_inputs <- function(trait_data, phylo_tree, env_data, species_id) {
  if (! (is.matrix(trait_data) || is.data.frame(trait_data)))
    stop("`trait_data` must be a matrix or data.frame.")
  if (! (inherits(phylo_tree, "phylo") || is.matrix(phylo_tree)))
    stop("`phylo_tree` must be ape::phylo or an adjacency matrix.")
  if (nrow(env_data) != length(species_id) ||
      nrow(env_data) != nrow(as.matrix(trait_data)))
    stop("Rows of `env_data`, `species_id` & `trait_data` must match.")
}

#— 2) Phylogeny → adjacency --------------------------------------------------
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

#— 3) Prepare tensors --------------------------------------------------------
prepare_tensors <- function(X, env_data, species_id, A_hat) {
  M_obs <- !is.na(X)
  X[!M_obs] <- 0
  X_t   <- torch_tensor(X, dtype = torch_float(), device = device)
  M_t   <- torch_tensor(M_obs, dtype = torch_bool(),  device = device)
  
  E_t   <- torch_tensor(as.matrix(env_data), dtype = torch_float(), device = device)
  
  sp_f   <- factor(species_id, levels = unique(species_id))
  sp_idx <- torch_tensor(as.integer(sp_f), dtype = torch_long(), device = device)
  sp_lev <- levels(sp_f)
  
  Xsp <- t(sapply(sp_lev, function(s){
    rows <- which(species_id == s)
    colMeans(X[rows,, drop = FALSE], na.rm = TRUE)
  }))
  X_sp_t <- torch_tensor(Xsp, dtype = torch_float(), device = device)
  
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
    A      = A_t
  )
}

#— 4) Loss: masked MSE -------------------------------------------------------
masked_mse <- function(pred, target, mask) {
  m <- mask$to(dtype = torch_float())
  torch_mean((pred - target)$pow(2) * m)
}

#— 5) VGAE building blocks with configurable hidden depth  ------------------

VGAE_Encoder <- nn_module(
  "VGAE_Encoder",
  initialize = function(inp, hid, lat, n_hidden = 1) {
    layers <- list(nn_linear(inp, hid), nn_relu())
    if (n_hidden > 1) {
      for (i in seq_len(n_hidden - 1)) {
        layers <- c(layers, list(nn_linear(hid, hid), nn_relu()))
      }
    }
    self$enc <- do.call(nn_sequential, layers)
    self$mu  <- nn_linear(hid, lat)
    self$lv  <- nn_linear(hid, lat)
  },
  forward = function(X_sp, A) {
    H0 <- torch_matmul(A, X_sp)
    H  <- self$enc(H0)
    list(mu = self$mu(H), logvar = self$lv(H))
  }
)

EnvEncoder <- nn_module(
  "EnvEncoder",
  initialize = function(envdim, hid, lat, n_hidden = 1) {
    layers <- list(nn_linear(envdim, hid), nn_relu())
    if (n_hidden > 1) {
      for (i in seq_len(n_hidden - 1)) {
        layers <- c(layers, list(nn_linear(hid, hid), nn_relu()))
      }
    }
    layers <- c(layers, list(nn_linear(hid, lat)))
    self$enc_env <- do.call(nn_sequential, layers)
  },
  forward = function(e) {
    self$enc_env(e)
  }
)

Decoder <- nn_module(
  "Decoder",
  initialize = function(lat, out) {
    self$fc <- nn_linear(lat, out)
  },
  forward = function(z) {
    self$fc(z)
  }
)

#— 6) Train + early‐stop + checkpoint + history -----------------------------
train_vgae_mc <- function(prep,
                          latent_dim,
                          epochs,
                          n_samp,
                          lr             = 1e-2,
                          patience       = 50,
                          ckpt_path      = NULL,
                          n_hidden_layers = 1) {
  
  # dims
  p       <- as.integer(prep$X$size(2))
  h       <- as.integer(p * 2)
  lat     <- as.integer(latent_dim)
  env_dim <- as.integer(prep$env$size(2))
  
  # build models with chosen depth
  enc_phy <- VGAE_Encoder(p, h, lat, n_hidden_layers)$to(device = device)
  enc_env <- EnvEncoder(env_dim, h, lat, n_hidden_layers)$to(device = device)
  dec     <- Decoder(lat, p)$to(device = device)
  
  opt <- optim_adam(
    c(enc_phy$parameters, enc_env$parameters, dec$parameters),
    lr = lr#,
    #weight_decay = 1e-4
  )
  
  if (!is.null(ckpt_path)) {
    dir.create(ckpt_path, recursive = TRUE, showWarnings = FALSE)
    best_file <- file.path(ckpt_path, "best_state.pt")
  }
  
  history       <- numeric(epochs)
  best_loss     <- Inf
  wait          <- 0
  stopped_epoch <- epochs
  
  for (ep in seq_len(epochs)) {
    enc_phy$train(); enc_env$train(); dec$train()
    
    Rphy <- enc_phy(prep$X_sp, prep$A)
    Renv <- enc_env(prep$env)
    zsp  <- Rphy$mu +
      torch_exp(0.5 * Rphy$logvar) * torch_randn_like(Rphy$mu)
    zobs <- zsp$index_select(1, prep$sp_idx) + Renv
    
    out      <- dec(zobs)
    loss_rec <- masked_mse(out, prep$X, prep$M)
    loss_kl  <- -0.5 * torch_mean(
      1 + Rphy$logvar - Rphy$mu$pow(2) - torch_exp(Rphy$logvar)
    )
    loss     <- loss_rec + loss_kl
    
    opt$zero_grad(); loss$backward(); opt$step()
    
    cur_loss     <- loss$item()
    history[ep]  <- cur_loss
    
    if (cur_loss < best_loss) {
      best_loss     <- cur_loss
      wait          <- 0
      stopped_epoch <- ep
      if (!is.null(ckpt_path)) {
        torch_save(
          list(
            enc_phy = enc_phy$state_dict(),
            enc_env = enc_env$state_dict(),
            dec     = dec$state_dict()
          ),
          best_file
        )
      }
    } else {
      wait <- wait + 1
      if (wait >= patience) {
        cat(sprintf("Early stopping at epoch %d (best=%.4f)\n",
                    ep, best_loss))
        break
      }
    }
    
    if (ep %% 200 == 0) {
      cat(sprintf("Epoch %4d | loss=%.4f\n", ep, cur_loss))
    }
  }
  
  # reload best
  if (!is.null(ckpt_path)) {
    st <- torch_load(best_file)
    enc_phy$load_state_dict(st$enc_phy)
    enc_env$load_state_dict(st$enc_env)
    dec$load_state_dict(st$dec)
  }
  
  list(
    enc_phy = enc_phy,
    enc_env = enc_env,
    dec     = dec,
    history = history[seq_len(ep)],
    stopped = stopped_epoch
  )
}

#— 7) Wrapper: standardise ↔ train ↔ back-transform --------------------------
#’ Impute missing trait values using a phylogenetic VGAE with flexible depth
#’ @param n_hidden_layers Number of hidden layers in each encoder (default 1)
impute_phylo <- function(trait_data,
                         phylo_tree,
                         env_data,
                         species_id,
                         latent_dim     = 32,
                         epochs         = 1000,
                         n_samples      = 100,
                         lr             = 1e-2,
                         patience       = 50,
                         ckpt_path      = "checkpoints",
                         n_hidden_layers = 1,
                         mask_obs_uncertainty = TRUE,
                         restore_observed     = TRUE) {
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  
  X_raw <- as.matrix(trait_data)
  mu    <- colMeans( X_raw, na.rm = TRUE)
  sigma <- apply( X_raw, 2, sd, na.rm = TRUE)
  X_std <- sweep(sweep(X_raw, 2, mu, "-"), 2, sigma, "/")
  
  A_hat <- normalize_adj(tree_to_graph(phylo_tree))
  prep  <- prepare_tensors(X_std, env_data, species_id, A_hat)
  
  fit <- train_vgae_mc(
    prep,
    latent_dim      = latent_dim,
    epochs          = epochs,
    n_samp          = n_samples,
    lr              = lr,
    patience        = patience,
    ckpt_path       = ckpt_path,
    n_hidden_layers = n_hidden_layers
  )
  
  enc_phy <- fit$enc_phy; enc_env <- fit$enc_env; dec <- fit$dec
  
  enc_phy$eval(); enc_env$eval(); dec$eval()
  n_obs <- prep$X$size(1)
  pmat  <- array(NA_real_, c(n_samples, n_obs, length(mu)))
  for (i in seq_len(n_samples)) {
    Rphy <- enc_phy(prep$X_sp, prep$A)
    Renv <- enc_env(prep$env)
    zsp  <- Rphy$mu +
      torch_exp(0.5 * Rphy$logvar) * torch_randn_like(Rphy$mu)
    zobs <- zsp$index_select(1, prep$sp_idx) + Renv
    pmat[i,,] <- as_array(dec(zobs))
  }
  
  mean_std <- apply(pmat, 2:3, mean)
  sd_std   <- apply(pmat, 2:3, sd)
  mean_m   <- sweep(sweep(mean_std, 2, sigma, "*"), 2, mu, "+")
  sd_m     <- sweep(sd_std,   2, sigma, "*")
  
  obs_mask <- !is.na(X_raw)
  if (mask_obs_uncertainty) sd_m  [obs_mask] <- 0
  if (restore_observed)    mean_m[obs_mask] <- X_raw[obs_mask]
  
  miss_idx <- which(is.na(X_raw), arr.ind = TRUE)
  imputed_missing <- data.frame(
    obs_id      = miss_idx[,1],
    species_id  = species_id[ miss_idx[,1] ],
    trait       = colnames(trait_data)[ miss_idx[,2] ],
    imputed     = mean_m[miss_idx],
    uncertainty = sd_m  [miss_idx],
    stringsAsFactors = FALSE
  )
  
  list(
    original_data   = as.data.frame(X_raw),
    completed_data  = as.data.frame(mean_m),
    uncertainty     = as.data.frame(sd_m),
    imputed_missing = imputed_missing,
    history         = fit$history,
    stopped_epoch   = fit$stopped
  )
}


#
#— 8) Quick test on simulated data ------------------------------------------
set.seed(42)
n_sp <- 500
tree <- rtree(n_sp)
V_phy <- cov2cor(vcv(tree))

# simulate env
p_env <- 10
R_env <- (matrix(runif(p_env^2, .05, .8), p_env, p_env) +
            t(matrix(runif(p_env^2, .05, .8), p_env, p_env))) / 2
diag(R_env) <- 1
env_data <- as.data.frame(mvrnorm(n_sp, rep(0,p_env), nearPD(R_env)$mat))

# simulate phylo‐traits + env effects
p_tr <- 10
Z    <- sapply(1:p_tr, function(i) mvrnorm(1, rep(0, n_sp), V_phy))
Sigma_tr <- nearPD((matrix(runif(p_tr^2, .05, .8), p_tr, p_tr) +
                      t(matrix(runif(p_tr^2, .05, .8), p_tr, p_tr))) / 2)$mat
Xsp <- Z %*% chol(Sigma_tr)
sobs_traits <- Xsp + as.matrix(env_data) %*% t(B) +
  matrix(rnorm(n_sp*p_tr,0, .1), n_sp, p_tr)
colnames(obs_traits) <- paste0("trait", 1:p_tr)

# mask 30%
mask_mat <- matrix(runif(n_sp*p_tr), n_sp, p_tr) < 0.30
obs_traits_mask <- obs_traits
obs_traits_mask[mask_mat] <- NA

# impute & evaluate
obs_traits_mask_mat <- as.matrix(obs_traits_mask)


res <- impute_phylo(
  trait_data   = as.data.frame(obs_traits_mask_mat),
  phylo_tree   = tree,
  env_data     = env_data,
  species_id   = tree$tip.label,
  latent_dim   = 128,
  epochs       = 6000,
  n_samples    = 200,
  lr           = 1e-2,
  patience     = 300,
  ckpt_path    = "checkpoints",
  n_hidden_layers = 1
)

######

# diagnostics
cat("Stopped at epoch:", res$stopped_epoch, "\n")
df_hist <- data.frame(
  epoch = seq_along(res$history),
  loss  = res$history
)
ggplot(df_hist, aes(epoch, loss)) +
  geom_line() +
  geom_vline(xintercept = res$stopped_epoch, linetype="dashed") +
  labs(title = "Training loss & early‐stopping point")

# RMSE
completed <- as.matrix(res$completed_data)
truth     <- obs_traits
cat("RMSE on masked entries:", rmse(truth[mask_mat], completed[mask_mat]), "\n")

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

