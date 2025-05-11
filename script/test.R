# this works for continous variables

library(torch)
library(ape)

#— 1. Validation ------------------------------------------------------------
validate_inputs <- function(trait_data, phylo_tree, env_data, species_id) {
  if (! (is.matrix(trait_data) || is.data.frame(trait_data)))
    stop("`trait_data` must be a matrix or data.frame.")
  if (! (inherits(phylo_tree, "phylo") || is.matrix(phylo_tree)))
    stop("`phylo_tree` must be an ape::phylo object or a numeric matrix.")
  if (nrow(env_data) != length(species_id) ||
      nrow(env_data) != nrow(as.matrix(trait_data)))
    stop("`env_data`, `species_id`, and `trait_data` must have the same number of rows.")
}

#— 2. Phylogeny → adjacency -----------------------------------------------
#— Helper: Phylo tree → adjacency via cophenetic distances ---------------
tree_to_graph <- function(phylo_tree) {
  if (inherits(phylo_tree, "phylo")) {
    D <- cophenetic(phylo_tree)
    A <- 1 / (1 + D)
    diag(A) <- 0
    # *** Ensure the adjacency keeps tip labels as names ***
    dimnames(A) <- list(phylo_tree$tip.label, phylo_tree$tip.label)
    return(A)
  } else {
    # If it’s already a matrix, leave it alone (but you can name it if you like)
    return(phylo_tree)
  }
}

#— Helper: Normalize adjacency --------------------------------------------
normalize_adj <- function(A) {
  # A should already carry dimnames from tree_to_graph
  nm <- dimnames(A)
  D_inv_sqrt <- diag(1 / sqrt(rowSums(A)))
  D_inv_sqrt[!is.finite(D_inv_sqrt)] <- 0
  B <- D_inv_sqrt %*% A %*% D_inv_sqrt
  if (!is.null(nm)) dimnames(B) <- nm
  B
}

#— 3. Data → torch tensors -----------------------------------------------
prepare_tensors <- function(trait_data, env_data, species_id, A_hat) {
  # — Observation‐level traits & mask
  X_obs    <- as.matrix(trait_data)
  mask_obs <- !is.na(X_obs)
  X_obs[!mask_obs] <- 0
  X_t      <- torch_tensor(X_obs, dtype = torch_float())
  M_t      <- torch_tensor(mask_obs, dtype = torch_float())
  
  # — Environmental covariates
  E_t <- torch_tensor(as.matrix(env_data), dtype = torch_float())
  
  # — Species‐level indexing
  sp_f    <- factor(species_id, levels = unique(species_id))
  sp_idx  <- torch_tensor(as.integer(sp_f), dtype = torch_long())
  sp_lev  <- levels(sp_f)
  
  # — Species‐level means & mask
  Xsp <- t(sapply(sp_lev, function(s) {
    rows <- which(species_id == s)
    colMeans(X_obs[rows, , drop = FALSE], na.rm = TRUE)
  }))
  mask_sp <- !is.na(Xsp)
  Xsp[!mask_sp] <- 0
  Xsp_t <- torch_tensor(Xsp, dtype = torch_float())
  Msp_t <- torch_tensor(mask_sp, dtype = torch_float())
  
  # — Subset adjacency to species‐level (silently by name if possible)
  if (!is.null(rownames(A_hat))) {
    idx <- match(sp_lev, rownames(A_hat))
    if (all(!is.na(idx))) {
      A_sub <- A_hat[idx, idx, drop = FALSE]
    } else {
      A_sub <- A_hat[seq_along(sp_lev), seq_along(sp_lev), drop = FALSE]
    }
  } else {
    A_sub <- A_hat[seq_along(sp_lev), seq_along(sp_lev), drop = FALSE]
  }
  A_t <- torch_tensor(A_sub, dtype = torch_float())
  
  # — Return exactly the slots trainers expect
  list(
    X      = X_t,
    M      = M_t,
    X_sp   = Xsp_t,
    M_sp   = Msp_t,
    env    = E_t,
    sp_idx = sp_idx,
    A      = A_t
  )
}
#— 4. Masked MSE helper ----------------------------------------------------
masked_mse <- function(pred, target, mask) {
  torch_mean(((pred - target)^2)[ mask$to(dtype=torch_bool()) ])
}

#— 5. Model definitions ----------------------------------------------------
VGAE_Encoder <- nn_module(
  initialize = function(inp, hid, lat) {
    self$gc1 <- nn_linear(inp, hid)
    self$mu  <- nn_linear(hid, lat)
    self$lv  <- nn_linear(hid, lat)
  },
  forward = function(X, A) {
    H <- torch_relu(self$gc1(torch_matmul(A, X)))
    list(mu = self$mu(H), logvar = self$lv(H))
  }
)
EnvEncoder <- nn_module(
  initialize = function(envdim, hid, lat) {
    self$f1 <- nn_linear(envdim, hid)
    self$f2 <- nn_linear(hid, lat)
  },
  forward = function(e) torch_relu(self$f1(e)) %>% self$f2()
)
Decoder <- nn_module(
  initialize = function(lat, out) self$fc <- nn_linear(lat, out),
  forward = function(z) self$fc(z)
)
build_dropout_ae <- function(input_dim, hidden_dim = 64, dropout = 0.3) {
  nn_module(
    initialize = function() {
      self$enc  <- nn_linear(input_dim, hidden_dim)
      self$drop <- nn_dropout(p = dropout)
      self$dec  <- nn_linear(hidden_dim, input_dim)
    },
    forward = function(x) {
      # respect train()/eval() state of the module
      z <- torch_relu(self$enc(x))
      z <- self$drop(z)         # no training= argument here
      self$dec(z)
    }
  )()
}

#— 6. Trainers -------------------------------------------------------------
train_vgae <- function(prep, lat, epochs, nsamp) {
  p <- prep$X$size(2); h <- p*2
  enc_phy <- VGAE_Encoder(p, h, lat)
  enc_env <- EnvEncoder(prep$env$size(2), h, lat)
  dec     <- Decoder(lat, p)
  opt     <- optim_adam(c(enc_phy$parameters, enc_env$parameters, dec$parameters), lr=1e-2)
  for (ep in 1:epochs) {
    enc_phy$train(); enc_env$train(); dec$train()
    res <- enc_phy(prep$X_sp, prep$A)
    zsp <- res$mu + torch_exp(0.5*res$logvar) * torch_randn_like(res$mu)
    z   <- zsp$index_select(1, prep$sp_idx) + enc_env(prep$env)
    out <- dec(z)
    loss_rec <- masked_mse(out, prep$X, prep$M)
    loss_kl  <- -0.5*torch_mean(1 + res$logvar - res$mu^2 - torch_exp(res$logvar))
    (loss <- loss_rec + loss_kl)
    opt$zero_grad(); loss$backward(); opt$step()
    if (ep %% 50 == 0) cat(sprintf("VGAE %d loss=%.4f\n", ep, loss$item()))
  }
  # posterior sampling
  enc_phy$eval(); enc_env$eval(); dec$eval()
  pmat <- array(NA_real_, c(nsamp, prep$X$size(1), p))
  for (i in 1:nsamp) {
    res <- enc_phy(prep$X_sp, prep$A)
    zsp <- res$mu + torch_exp(0.5*res$logvar)*torch_randn_like(res$mu)
    z   <- zsp$index_select(1, prep$sp_idx) + enc_env(prep$env)
    pmat[i,,] <- as_array(dec(z))
  }
  list(mean=apply(pmat,2:3,mean), sd=apply(pmat,2:3,sd))
  
}


#— 7. Master wrapper -------------------------------------------------------
#' Impute missing traits with phylogeny/environment‐aware autoencoders
#'
#' @param trait_data         data.frame or matrix of traits (with NA)
#' @param phylo_tree         ape::phylo or adjacency matrix
#' @param env_data           data.frame or matrix of environmental covariates
#' @param species_id         vector of species names (matching tree tips)
#' @param latent_dim         latent dimension for all models
#' @param epochs             training epochs
#' @param n_samples          samples for uncertainty (VGAE & dropout)
#' @param n_models           ensemble size
#' @param dropout            dropout rate for dropout/ensemble
#' @param mask_obs_uncertainty  logical; if TRUE, set SD=0 for original (non‐missing) data
#' @param restore_observed     logical; if TRUE, overwrite reconstructed points with original values
#' @return a list with  
#'   - `completed_data`: data.frame of imputed/completed values  
#'   - `uncertainty`:   data.frame of SDs  
#' @export
impute_phylo <- function(trait_data, phylo_tree, env_data, species_id,
                         #method = c("vgae","mc_dropout","ensemble"),
                         latent_dim = 16, epochs = 300,
                         n_samples = 50, n_models = 5, dropout = 0.3,
                         mask_obs_uncertainty = TRUE,
                         restore_observed     = TRUE) {
  
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  A_hat <- normalize_adj(tree_to_graph(phylo_tree))
  prep  <- prepare_tensors(trait_data, env_data, species_id, A_hat)
  
  # call the appropriate trainer
  out <- train_vgae(prep, latent_dim, epochs, n_samples)
              
  
  # assemble mean & SD matrices
  mean_mat <- out$mean
  sd_mat   <- out$sd
  
  # build logical mask of originally observed entries
  obs_mask <- !is.na(as.matrix(trait_data))
  
  # 1) zero‐out SD for observed if requested
  if (mask_obs_uncertainty) {
    sd_mat[obs_mask] <- 0
  }
  
  # 2) restore original values in completed_data if requested
  if (restore_observed) {
    mean_mat[obs_mask] <- as.matrix(trait_data)[obs_mask]
  }
  
  # convert to data.frames with proper names
  completed <- as.data.frame(mean_mat)
  names(completed) <- colnames(as.matrix(trait_data))
  
  uncertainty <- as.data.frame(sd_mat)
  names(uncertainty) <- names(completed)
  
  list(
    completed_data = completed,
    uncertainty    = uncertainty
  )
}
# ===================================================================
# Example quick test  (comment out if using in a package)
# ===================================================================


n_sp <- 30; n_obs <- 80
tree   <- rtree(n_sp); tree$tip.label <- paste0("sp",1:n_sp)
trait1 <- fastBM(tree); trait2 <- fastBM(tree, sig2=0.5)
species_id <- sample(tree$tip.label, n_obs, TRUE)
env        <- rnorm(n_obs,10,5)
obs_traits <- data.frame(
  trait1 = trait1[species_id] + 0.3*env + rnorm(n_obs,0,0.2),
  trait2 = trait2[species_id] - 0.5*env + rnorm(n_obs,0,0.3)
)
# introduce ~20% missing
for (col in names(obs_traits)) {
  idx <- sample(n_obs, size = 0.2*n_obs)
  obs_traits[idx, col] <- NA
}

res <- impute_phylo(
  trait_data  = obs_traits,
  phylo_tree  = tree,
  env_data    = data.frame(env = env),
  species_id  = species_id,
  latent_dim  = 8,
  epochs      = 250,
  n_samples   = 100
)


