###############################################################################
## Enhanced Phylogenetic VGAE imputer
##   – MC-dropout  · Early-stopping  · Checkpointing  · LR scheduling
###############################################################################

library(torch)
#   optional: devtools::install_github("torch/mlverse/torch-sparse")
have_sparse <- requireNamespace("torch.sparse", quietly = TRUE)

library(ape)
library(Matrix)    # nearPD, sparse ops
library(MASS)      # mvrnorm
library(Metrics)   # rmse
library(ggplot2)

# ── 0) device ────────────────────────────────────────────────────────────────
device <- if (cuda_is_available()) torch_device("cuda") else
  if (backends_mps_is_available()) torch_device("mps") else
    torch_device("cpu")

# ── 1) validation ------------------------------------------------------------
validate_inputs <- function(trait_data, phylo_tree, env_data, species_id) {
  if (!(is.matrix(trait_data) || is.data.frame(trait_data)))
    stop("`trait_data` must be matrix or data.frame")
  if (!(inherits(phylo_tree, "phylo") || is.matrix(phylo_tree)))
    stop("`phylo_tree` must be ape::phylo or adjacency matrix")
  if (nrow(env_data) != length(species_id) ||
      nrow(env_data) != nrow(as.matrix(trait_data)))
    stop("Rows of env_data, species_id & trait_data must match")
  invisible(TRUE)
}

# ── 2) phylogeny → adjacency -------------------------------------------------
normalize_adjacency <- function(A) {
  Dinv <- diag(1 / sqrt(rowSums(A)))
  Dinv[!is.finite(Dinv)] <- 0
  B <- Dinv %*% A %*% Dinv
  dimnames(B) <- dimnames(A)
  B
}

phylo_to_adjacency <- function(tree, sparse = FALSE) {
  A <- if (inherits(tree, "phylo")) {
    D <- cophenetic(tree)
    diag(D) <- 0
    1 / (1 + D)
  } else tree
  diag(A) <- 0
  A <- normalize_adjacency(A)
  if (sparse) as(A, "dgCMatrix") else A
}

# ── 3) helper to build hidden stacks ----------------------------------------
build_hidden <- function(inp, hid, n, dropout) {
  layers <- list()
  for (i in seq_len(n)) {
    layers <- append(layers, list(
      nn_linear(if (i == 1) inp else hid, hid), nn_relu(),
      nn_dropout(dropout)))
  }
  do.call(nn_sequential, layers)
}

he_init <- function(m) {
  if (inherits(m, "nn_linear"))
    nn_init_kaiming_normal_(m$weight, nonlinearity = "relu")
}

# ── 4) modules ---------------------------------------------------------------
VGAE_Encoder <- nn_module(
  initialize = function(inp, hid, lat, n_hidden, dropout) {
    self$hidden <- build_hidden(inp, hid, n_hidden, dropout)
    self$mu     <- nn_linear(hid, lat)
    self$lv     <- nn_linear(hid, lat)
    self$apply(he_init)
  },
  forward = function(X_sp, A) {
    H0 <- if (have_sparse && inherits(A, "torch_sparse_coo_tensor"))
      torch.sparse::matmul(A, X_sp)
    else torch_mm(A, X_sp)
    H   <- self$hidden(H0)
    list(mu = self$mu(H), logvar = self$lv(H))
  }
)

EnvEncoder <- nn_module(
  initialize = function(inp, hid, lat, n_hidden, dropout) {
    self$hidden <- build_hidden(inp, hid, n_hidden, dropout)
    self$out    <- nn_linear(hid, lat)
    self$apply(he_init)
  },
  forward = function(e) self$out(self$hidden(e))
)

Decoder <- nn_module(
  initialize = function(lat, out_dim, dropout) {
    self$drop <- nn_dropout(dropout)
    self$fc   <- nn_linear(lat, out_dim)
    self$apply(he_init)
  },
  forward = function(z) self$fc(self$drop(z))
)

masked_mse <- function(pred, target, mask) {
  mask <- mask$to(dtype = torch_float())
  torch_mean((pred - target)$pow(2) * mask)
}

# ── 5) tensor prep -----------------------------------------------------------
prepare_tensors <- function(X, env, sp_id, A, val_split = 0.2, sparse_adj = FALSE) {
  
  Mfull <- !is.na(X)
  X[!Mfull] <- 0
  
  ## train / val masks --------------------------------------------------------
  obs <- which(Mfull, arr.ind = TRUE)
  val_n <- floor(nrow(obs) * val_split)
  val_i <- obs[sample(nrow(obs), val_n), , drop = FALSE]
  
  Mtrain <- Mfull; Mtrain[val_i] <- FALSE
  Mval   <- matrix(FALSE, nrow(X), ncol(X)); Mval[val_i] <- TRUE
  
  ## tensors ------------------------------------------------------------------
  tns <- list(
    X        = torch_tensor(X,        dtype = torch_float(), device = device),
    M_train  = torch_tensor(Mtrain,   dtype = torch_bool(),  device = device),
    M_val    = torch_tensor(Mval,     dtype = torch_bool(),  device = device),
    env      = torch_tensor(as.matrix(env), dtype = torch_float(), device = device)
  )
  
  sp_f <- factor(sp_id, levels = unique(sp_id))
  tns$sp_idx <- torch_tensor(as.integer(sp_f), dtype = torch_long(), device = device)
  
  Xsp <- t(sapply(levels(sp_f), function(s) {
    rows <- which(sp_id == s)
    colMeans(X[rows, , drop = FALSE], na.rm = TRUE)
  }))
  tns$X_sp <- torch_tensor(Xsp, dtype = torch_float(), device = device)
  
  idx <- match(levels(sp_f), rownames(A))
  A_sub <- A[idx, idx]
  
  tns$A <- if (sparse_adj) {
    Acoo <- as(A_sub, "dgTMatrix")
    torch_sparse_coo_tensor(
      indices = torch_tensor(rbind(Acoo@i, Acoo@j), dtype = torch_long()),
      values  = torch_tensor(Acoo@x,          dtype = torch_float()),
      size    = dim(A_sub),
      device  = device
    )
  } else torch_tensor(A_sub, dtype = torch_float(), device = device)
  tns
}

# ── 6) training loop ---------------------------------------------------------
train_vgae_mc <- function(prep, latent_dim, epochs, lr, patience,
                          n_hidden, dropout, ckpt, scheduler_pat) {
  
  p       <- prep$X$size(2)
  h_dim   <- p * 2
  env_dim <- prep$env$size(2)
  
  enc_phy <- VGAE_Encoder(p, h_dim, latent_dim, n_hidden, dropout)$to(device)
  enc_env <- EnvEncoder(env_dim, h_dim, latent_dim, n_hidden, dropout)$to(device)
  dec     <- Decoder(latent_dim, p, dropout)$to(device)
  
  opt <- optim_adam(
    c(enc_phy$parameters, enc_env$parameters, dec$parameters),
    lr = lr, weight_decay = 1e-4)
  
  sched <- lr_scheduler_reduce_on_plateau(opt, factor = 0.5,
                                          patience = scheduler_pat)
  
  if (!is.null(ckpt)) {
    dir.create(ckpt, showWarnings = FALSE)
    best_file <- file.path(ckpt, "best.pt")
  }
  
  hist <- data.frame(epoch = integer(0), train = numeric(0), val = numeric(0))
  best <- Inf; wait <- 0; stop_ep <- epochs
  
  for (ep in seq_len(epochs)) {
    
    ###  train ---------------------------------------------------------------
    enc_phy$train(); enc_env$train(); dec$train()
    
    Rphy <- enc_phy(prep$X_sp, prep$A)
    Renv <- enc_env(prep$env)
    zsp  <- Rphy$mu + torch_exp(0.5 * Rphy$logvar) * torch_randn_like(Rphy$mu)
    zobs <- zsp$index_select(1, prep$sp_idx) + Renv
    pred <- dec(zobs)
    
    l_rec <- masked_mse(pred, prep$X, prep$M_train)
    l_kl  <- -0.5 * torch_mean(1 + Rphy$logvar -
                                 Rphy$mu$pow(2) - torch_exp(Rphy$logvar))
    loss  <- l_rec + l_kl
    
    opt$zero_grad(); loss$backward(); opt$step()
    
    ###  validation ----------------------------------------------------------
    enc_phy$eval(); enc_env$eval(); dec$eval()
    with_no_grad({
      Rphy_v <- enc_phy(prep$X_sp, prep$A)
      Renv_v <- enc_env(prep$env)
      zsp_v  <- Rphy_v$mu + torch_exp(0.5 * Rphy_v$logvar) *
        torch_randn_like(Rphy_v$mu)
      zobs_v <- zsp_v$index_select(1, prep$sp_idx) + Renv_v
      pred_v <- dec(zobs_v)
      val_loss <- masked_mse(pred_v, prep$X, prep$M_val)
    })
    
    ## history row -----------------------------------------------------------
    hist <- rbind(hist,
                  data.frame(epoch = ep,
                             train = loss$item(),
                             val   = val_loss$item()))
    ## scheduler & early-stop -------------------------------------------------
    sched$step(val_loss)
    if (val_loss$item() < best) {
      best <- val_loss$item(); wait <- 0; stop_ep <- ep
      if (!is.null(ckpt))
        torch_save(list(enc_phy = enc_phy$state_dict(),
                        enc_env = enc_env$state_dict(),
                        dec     = dec$state_dict()), best_file)
    } else {
      wait <- wait + 1
      if (wait >= patience) {
        message("Early stop @ epoch ", ep,
                " (best val = ", format(best, digits = 4), ")")
        break
      }
    }
    
    if (ep %% 100 == 0)
      cat(sprintf("Ep %4d | tr %.4f | val %.4f | lr %.2e\n",
                  ep, hist$train[ep], hist$val[ep], opt$param_groups[[1]]$lr))
  }
  
  if (!is.null(ckpt)) {
    st <- torch_load(best_file)
    enc_phy$load_state_dict(st$enc_phy)
    enc_env$load_state_dict(st$enc_env)
    dec$load_state_dict(st$dec)
  }
  
  list(enc_phy = enc_phy, enc_env = enc_env, dec = dec,
       history = hist, stopped = stop_ep)
}

# ── 7) main wrapper ----------------------------------------------------------
impute_phylo <- function(
    trait_data, phylo_tree, env_data, species_id,
    latent_dim = 32, epochs = 2000, n_samples = 100,
    lr = 1e-3, patience = 100, val_split = 0.2,
    ckpt_path = "checkpoints", n_hidden = 1,
    dropout_rate = 0.2, sparse_adj = FALSE,
    scheduler_patience = 10) {
  
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  
  Xraw <- as.matrix(trait_data)
  mu    <- colMeans(Xraw, na.rm = TRUE)
  sigma <- apply(Xraw, 2, sd, na.rm = TRUE)
  Xstd  <- sweep(sweep(Xraw, 2, mu, "-"), 2, sigma, "/")
  
  Adj   <- phylo_to_adjacency(phylo_tree, sparse_adj)
  prep  <- prepare_tensors(Xstd, env_data, species_id, Adj,
                           val_split, sparse_adj)
  
  fit <- train_vgae_mc(
    prep, latent_dim, epochs, lr, patience,
    n_hidden, dropout_rate, ckpt_path, scheduler_patience)
  
  fit$enc_phy$train(); fit$enc_env$train(); fit$dec$train()
  
  n <- prep$X$size(1); p <- prep$X$size(2)
  samp <- array(NA_real_, c(n_samples, n, p))
  
  for (i in seq_len(n_samples)) {
    with_no_grad({
      Rp <- fit$enc_phy(prep$X_sp, prep$A)
      Re <- fit$enc_env(prep$env)
      z  <- Rp$mu + torch_exp(0.5 * Rp$logvar) * torch_randn_like(Rp$mu)
      z  <- z$index_select(1, prep$sp_idx) + Re
      samp[i,,] <- as_array(fit$dec(z))
    })
  }
  
  mean_std <- apply(samp, 2:3, mean)
  sd_std   <- apply(samp, 2:3, sd)
  
  mean_imp <- sweep(sweep(mean_std, 2, sigma, "*"), 2, mu, "+")
  sd_imp   <- sweep(sd_std,   2, sigma, "*")
  
  obs_mask <- !is.na(Xraw)
  mean_imp[obs_mask] <- Xraw[obs_mask]
  sd_imp  [obs_mask] <- 0
  
  miss <- which(is.na(Xraw), arr.ind = TRUE)
  imputed <- data.frame(
    obs_id     = miss[,1],
    species_id = species_id[miss[,1]],
    trait      = colnames(Xraw)[miss[,2]],
    imputed    = mean_imp[miss],
    uncertainty= sd_imp[miss]
  )
  
  list(imputed = imputed,
       completed = as.data.frame(mean_imp),
       uncertainty = as.data.frame(sd_imp),
       history = fit$history,
       stopped_epoch = fit$stopped)
}

# ── 8) tiny smoke-test -------------------------------------------------------
if (interactive()) {
  set.seed(42)
  n_sp <- 200
  tree <- rtree(n_sp)
  
  ## env
  p_env <- 5
  R_env <- diag(p_env)
  env   <- mvrnorm(n_sp, rep(0, p_env), R_env)
  
  ## traits
  p_tr <- 5
  Vphy <- cov2cor(vcv(tree))
  Z    <- replicate(p_tr, mvrnorm(1, rep(0, n_sp), Vphy))
  Sigma_tr <- diag(p_tr)
  Xsp <- Z %*% chol(Sigma_tr)
  
  B <- matrix(runif(p_tr * p_env, -1, 1), p_tr, p_env)
  obs <- Xsp + env %*% t(B) + matrix(rnorm(n_sp * p_tr, 0, 0.1), n_sp, p_tr)
  colnames(obs) <- paste0("T", 1:p_tr)
  
  mask <- matrix(runif(n_sp * p_tr) < 0.3, n_sp, p_tr)
  obs[mask] <- NA
  
  res <- impute_phylo(obs, tree, env, tree$tip.label,
                      latent_dim = 64, epochs = 1000,
                      patience = 100, n_hidden = 2,
                      dropout_rate = 0.3, lr = 1e-3,
                      ckpt_path = tempfile())
  
  cat("early-stop @", res$stopped_epoch, "\n")
}