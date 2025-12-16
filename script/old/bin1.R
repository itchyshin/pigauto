# conceputal

# trait_data ──▶|         |─▶ zsp (phylo latent) ──┐
# phylo_graph──▶| enc_phy |──────────────────────┬─┴──▶ zobs ──▶ Decoder ──▶ Imputed traits
# env_data ─────▶| enc_env |─▶ Renv ──────────────┘


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



# binary 

# ---- Helper loss function for binary traits ----
masked_bce <- function(logit, target, mask) {
  # Target must be float (0/1), logit is real-valued
  m <- mask$to(dtype = torch_float())
  # nnf_binary_cross_entropy_with_logits expects float inputs
  loss <- nnf_binary_cross_entropy_with_logits(
    logit, target$to(dtype = torch_float()), reduction = "none"
  )
  torch_mean(loss * m)
}

# ---- Wrapper for binary phylogenetic imputation ----
impute_phylo_binary <- function(trait_data,
                                phylo_tree,
                                env_data,
                                species_id,
                                latent_dim     = 32,
                                epochs         = 1000,
                                n_samples      = 100,
                                lr             = 1e-2,
                                patience       = 50,
                                ckpt_path      = "checkpoints_bin",
                                n_hidden_layers = 1,
                                mask_obs_uncertainty = TRUE,
                                restore_observed     = TRUE) {
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  
  # Check all columns are binary (0/1)
  X_raw <- as.matrix(trait_data)
  if (!all(na.omit(as.vector(X_raw)) %in% 0:1)) stop("All values must be 0/1.")
  
  # NO standardisation for binary
  X_bin <- X_raw
  
  A_hat <- normalize_adj(tree_to_graph(phylo_tree))
  prep  <- prepare_tensors(X_bin, env_data, species_id, A_hat)
  
  # --- Modified train loop for binary ---
  train_vgae_binary <- function(prep,
                                latent_dim,
                                epochs,
                                n_samp,
                                lr             = 1e-2,
                                patience       = 50,
                                ckpt_path      = NULL,
                                n_hidden_layers = 1) {
    p       <- as.integer(prep$X$size(2))
    h       <- as.integer(p * 2)
    lat     <- as.integer(latent_dim)
    env_dim <- as.integer(prep$env$size(2))
    enc_phy <- VGAE_Encoder(p, h, lat, n_hidden_layers)$to(device = device)
    enc_env <- EnvEncoder(env_dim, h, lat, n_hidden_layers)$to(device = device)
    dec     <- Decoder(lat, p)$to(device = device)
    opt <- optim_adam(
      c(enc_phy$parameters, enc_env$parameters, dec$parameters),
      lr = lr
    )
    if (!is.null(ckpt_path)) {
      dir.create(ckpt_path, recursive = TRUE, showWarnings = FALSE)
      best_file <- file.path(ckpt_path, "best_state.pt")
    }
    history <- numeric(epochs)
    best_loss <- Inf
    wait <- 0
    stopped_epoch <- epochs
    for (ep in seq_len(epochs)) {
      enc_phy$train(); enc_env$train(); dec$train()
      Rphy <- enc_phy(prep$X_sp, prep$A)
      Renv <- enc_env(prep$env)
      zsp  <- Rphy$mu + torch_exp(0.5 * Rphy$logvar) * torch_randn_like(Rphy$mu)
      zobs <- zsp$index_select(1, prep$sp_idx) + Renv
      logits <- dec(zobs)
      loss_rec <- masked_bce(logits, prep$X, prep$M)
      loss_kl <- -0.5 * torch_mean(1 + Rphy$logvar - Rphy$mu$pow(2) - torch_exp(Rphy$logvar))
      loss <- loss_rec + loss_kl
      opt$zero_grad(); loss$backward(); opt$step()
      cur_loss <- loss$item()
      history[ep] <- cur_loss
      if (cur_loss < best_loss) {
        best_loss <- cur_loss
        wait <- 0
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
          cat(sprintf("Early stopping at epoch %d (best=%.4f)\n", ep, best_loss))
          break
        }
      }
      if (ep %% 200 == 0) {
        cat(sprintf("Epoch %4d | loss=%.4f\n", ep, cur_loss))
      }
    }
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
  
  fit <- train_vgae_binary(
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
  pmat  <- array(NA_real_, c(n_samples, n_obs, ncol(X_bin)))
  for (i in seq_len(n_samples)) {
    Rphy <- enc_phy(prep$X_sp, prep$A)
    Renv <- enc_env(prep$env)
    zsp  <- Rphy$mu + torch_exp(0.5 * Rphy$logvar) * torch_randn_like(Rphy$mu)
    zobs <- zsp$index_select(1, prep$sp_idx) + Renv
    logits <- dec(zobs)
    pmat[i,,] <- as_array(torch_sigmoid(logits))  # Get probabilities
  }
  # Average predictions, then threshold at 0.5 (or return probabilities)
  prob_m <- apply(pmat, 2:3, mean)
  imp_m  <- ifelse(prob_m >= 0.5, 1, 0)
  obs_mask <- !is.na(X_bin)
  if (restore_observed) imp_m[obs_mask] <- X_bin[obs_mask]
  miss_idx <- which(is.na(X_bin), arr.ind = TRUE)
  imputed_missing <- data.frame(
    obs_id      = miss_idx[,1],
    species_id  = species_id[ miss_idx[,1] ],
    trait       = colnames(trait_data)[ miss_idx[,2] ],
    imputed     = imp_m[miss_idx],
    prob        = prob_m[miss_idx],
    stringsAsFactors = FALSE
  )
  list(
    original_data   = as.data.frame(X_bin),
    completed_data  = as.data.frame(imp_m),
    prob            = as.data.frame(prob_m),
    imputed_missing = imputed_missing,
    history         = fit$history,
    stopped_epoch   = fit$stopped
  )
}

# test

# Example: trait_data is a data.frame of 0/1 (NA allowed)
res_bin <- impute_phylo_binary(
  trait_data   =  as.data.frame(traits_df_miss[,6:8]),        # Data frame of binary (0/1) traits
  phylo_tree   = tree,
  env_data     = env_data,
  species_id   = tree$tip.label,
  latent_dim   = 16,
  epochs       = 4000,
  n_samples    = 500,
  patience     = 500
)

res_bin2 <- impute_phylo_binary(
  trait_data   =  as.data.frame(traits_df_miss[,6:8]),        # Data frame of binary (0/1) traits
  phylo_tree   = tree2,
  env_data     = env_data,
  species_id   = tree2$tip.label,
  latent_dim   = 16,
  epochs       = 4000,
  n_samples    = 500,
  patience     = 500
)

completed_bin  <- as.matrix(res_bin$completed_data)
completed_bin2 <- as.matrix(res_bin2$completed_data)
truth_bin     <- traits_df[,6:8]
# RMSE only on masked entries (where miss_mat == TRUE)
cat("RMSE (tree):     ", rmse(truth_bin[miss_mat[,6:8]], completed_bin[miss_mat[,6:8]]), "\n")
cat("RMSE (star):     ", rmse(truth_bin[miss_mat[,6:8]], completed_bin2[miss_mat[,6:8]]), "\n")

# # For comparison, mean imputation (threshold colMeans at 0.5)
# mean_imp_bin <- matrix(
#   as.numeric(rep(colMeans(truth_bin, na.rm=TRUE), each=nrow(truth_bin)) >= 0.6),
#   nrow=nrow(truth_bin)
# )
# cat("RMSE (mean imp): ", rmse(truth_bin[miss_mat[,6:8]], mean_imp_bin[miss_mat[,6:8]]), "\n")


##

df_mask_bin <- data.frame(
  truth   = truth_bin[miss_mat[,6:8]],
  imputed = completed_bin[miss_mat[,6:8]]
)
ggplot(df_mask_bin, aes(x=truth, y=imputed)) +
  geom_jitter(alpha=0.3, width=0.2, height=0.2) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  coord_equal() + theme_minimal() +
  labs(title="Imputed vs True (masked, binary)", x="True", y="Imputed (0/1)")

