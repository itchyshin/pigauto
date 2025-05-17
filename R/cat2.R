library(torch)
library(ape)
library(Matrix)
library(Metrics)

# ─────────────────────────────────────────────────────────────
# Device Configuration (CRUCIAL FIX)
# ─────────────────────────────────────────────────────────────
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")

# ─────────────────────────────────────────────────────────────
# Helper Functions (Fixed Device Handling)
# ─────────────────────────────────────────────────────────────

categorical_nll <- function(logit, tgt, mask) {
  # Ensure consistent device
  tgt <- tgt$to(device = device, dtype = torch_long())
  mask <- mask$to(device = device, dtype = torch_float())
  
  l <- nnf_cross_entropy(logit, tgt, reduction = "none")
  torch_mean(l * mask)
}

prepare_tensors_cat <- function(trait_df, env_data, species_id,
                                A_hat, meta, levels_lst) {
  
  cat_cols <- meta$name
  
  trait_df[cat_cols] <- lapply(trait_df[cat_cols], function(x) {
    if (is.factor(x)) x else factor(x)
  })
  
  X_int <- do.call(cbind, lapply(cat_cols, function(nm) {
    v <- as.integer(trait_df[[nm]]) - 1L
    v[is.na(v)] <- 0L
    v
  }))
  
  mask_mat <- as.matrix(!is.na(trait_df[cat_cols]))
  stopifnot(is.logical(mask_mat))
  
  X_t   <- torch_tensor(X_int, dtype = torch_long(), device = device)
  M_t   <- torch_tensor(mask_mat, dtype = torch_bool(), device = device)
  env_t <- torch_tensor(as.matrix(env_data), dtype = torch_float(), device = device)
  
  sp_f   <- factor(species_id, levels = unique(species_id))
  sp_idx <- torch_tensor(as.integer(sp_f), dtype = torch_long(), device = device)
  
  Xsp <- t(sapply(levels(sp_f), function(s) {
    rows <- which(species_id == s)
    colMeans(X_int[rows, , drop = FALSE], na.rm = TRUE)
  }))
  Xsp_t <- torch_tensor(Xsp, dtype = torch_float(), device = device)
  
  idx <- match(levels(sp_f), rownames(A_hat))
  A_t <- torch_tensor(A_hat[idx, idx], dtype = torch_float(), device = device)
  
  list(
    X_sp       = Xsp_t,
    A          = A_t,
    env        = env_t,
    sp_idx     = sp_idx,
    mask_cat   = lapply(seq_len(M_t$size(2)), function(j) M_t[, j]$unsqueeze(2)),
    target_cat = lapply(seq_len(X_t$size(2)), function(j) X_t[, j]$unsqueeze(2))
  )
}

# ─────────────────────────────────────────────────────────────
# Neural Network Modules (Fixed Implementation)
# ─────────────────────────────────────────────────────────────

VGAE_Encoder <- nn_module(
  "VGAE_Encoder",
  initialize = function(inp, hid, lat, depth) {
    self$mlp <- nn_sequential(
      nn_linear(inp, hid),
      nn_relu(),
      nn_layer_norm(hid)
    )
    self$mu <- nn_linear(hid, lat)
    self$logvar <- nn_linear(hid, lat)
  },
  forward = function(Xsp, A) {
    H <- self$mlp(torch_matmul(A, Xsp))
    list(mu = self$mu(H), logvar = self$logvar(H))
  }
)

EnvEncoder <- nn_module(
  "EnvEncoder",
  initialize = function(idim, hid, lat, depth) {
    self$mlp <- nn_sequential(
      nn_linear(idim, hid),
      nn_relu(),
      nn_linear(hid, lat)
    )
  },
  forward = function(e) self$mlp(e)
)

BinaryHead <- nn_module(
  "BinaryHead",
  initialize = function(lat) self$lin <- nn_linear(lat, 1),
  forward = function(z) self$lin(z)
)

NominalHead <- nn_module(
  "NominalHead",
  initialize = function(lat, K) self$lin <- nn_linear(lat, K),
  forward = function(z) self$lin(z)
)

OrdinalHead <- nn_module(
  "OrdinalHead",
  initialize = function(lat, K) {
    self$beta <- nn_linear(lat, 1, bias = FALSE)
    self$cuts <- nn_parameter(torch_randn(K-1, device = device))
  },
  forward = function(z) torch_sigmoid(self$cuts - self$beta(z))
)

# ─────────────────────────────────────────────────────────────
# Loss Functions (Device-Consistent)
# ─────────────────────────────────────────────────────────────

ordinal_nll <- function(cumprob, tgt, mask, K) {
  # Ensure device consistency
  tgt <- tgt$to(device = cumprob$device)
  mask <- mask$to(device = cumprob$device)
  
  ones <- torch_ones(cumprob$size(1), 1, device = device)
  zeros <- torch_zeros_like(ones)
  Pcum <- torch_cat(list(ones, cumprob, zeros), dim = 2)
  pk <- Pcum[, 1:K] - Pcum[, 2:(K+1)]
  
  tgt_oh <- nnf_one_hot(tgt$to(dtype = torch_long()), num_classes = K)
  nll <- -torch_log(pk + 1e-7) * tgt_oh
  torch_mean(torch_sum(nll, dim = 2) * mask$to(dtype = torch_float()))
}

# ─────────────────────────────────────────────────────────────
# Core Imputation Function (Fixed)
# ─────────────────────────────────────────────────────────────

impute_phylo_cat <- function(trait_data, phylo_tree, env_data, species_id,
                             latent_dim = 32, hidden_dim = 64, depth = 1,
                             lr = 1e-2, epochs = 1000, patience = 100,
                             n_samples = 200, ckpt_path = "ckpt_cat",
                             mask_obs_uncertainty = TRUE, restore_observed = TRUE) {
  
  # Input validation
  validate_inputs <- function(td, pt, ed, si) {
    if (!inherits(pt, "phylo")) stop("Phylogeny must be ape::phylo object")
    if (nrow(ed) != length(si)) stop("Env data mismatch with species IDs")
  }
  
  infer_trait_types <- function(df) {
    data.frame(
      name = colnames(df),
      type = sapply(df, function(x) {
        if (is.numeric(x)) {
          if (all(na.omit(x) %in% 0:1)) "binary" else "cont"
        } else if (is.ordered(x)) "ordinal" else "nominal"
      }),
      n_level = sapply(df, function(x) if (is.factor(x)) nlevels(x) else NA)
    )
  }
  
  # Prepare data
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  trait_df <- as.data.frame(trait_data)
  meta <- infer_trait_types(trait_df)
  
  if (any(meta$type == "cont")) {
    stop("Continuous columns detected - use categorical-only version")
  }
  
  # Rest of the impute_phylo_cat function remains the same
  cat_cols   <- meta$name
  levels_lst <- lapply(trait_df[cat_cols], levels)
  A_hat      <- normalize_adj(tree_to_graph(phylo_tree))
  prep       <- prepare_tensors_cat(trait_df, env_data, species_id,
                                    A_hat, meta, levels_lst)
  
  # encoder + decoder heads
  enc_phy <- VGAE_Encoder(ncol(prep$X_sp), hidden_dim, latent_dim, depth)$to(device)
  enc_env <- EnvEncoder(ncol(prep$env),  hidden_dim, latent_dim, depth)$to(device)
  
  heads  <- nn_module_list()
  
  
  loss_f <- vector("list", nrow(meta))
  for (j in seq_len(nrow(meta))) {
    K  <- meta$n_level[j]
    tp <- meta$type[j]
    
    # Decoder head
    hd <- switch(tp,
                 binary  = BinaryHead(latent_dim),
                 nominal = NominalHead(latent_dim, K),
                 ordinal = OrdinalHead(latent_dim, K))
    heads$append(hd$to(device))
    
    # Main loss function
    main_loss <- switch(tp,
                        binary  = function(logit, tgt, mask) {
                          l <- nnf_binary_cross_entropy_with_logits(
                            logit,
                            tgt$to(dtype = torch_float())$unsqueeze(2),
                            reduction = "none"
                          )
                          torch_mean(l * mask$to(dtype = torch_float())$unsqueeze(2))
                        },
                        nominal = categorical_nll,
                        ordinal = function(cum, tgt, mask) ordinal_nll(cum, tgt, mask, K)
    )
    
    # L2 regularization loss (corrected)
    l2_loss <- function(out, tgt, mask) {
      torch_mean(out$pow(2))
    }
    
    loss_f[[j]] <- list(main = main_loss, reg = l2_loss)
  }
  
  
  opt  <- optim_adam(c(enc_phy$parameters, enc_env$parameters, heads$parameters),
                     lr = lr)
  best <- Inf; hist <- numeric(epochs); wait <- 0
  dir.create(ckpt_path, showWarnings = FALSE, recursive = TRUE)
  
  for (ep in seq_len(epochs)) {
    enc_phy$train(); enc_env$train()
    
    res  <- enc_phy(prep$X_sp, prep$A)
    z_sp <- res$mu + torch_exp(0.5 * res$logvar) * torch_randn_like(res$mu)
    z    <- z_sp$index_select(dim = 0, index = prep$sp_idx) +
      enc_env(prep$env)
    
    rec <- torch_tensor(0, device = device)
    for (j in seq_along(heads)) {
      out <- heads[[j]](z)
      # Sum all losses (main and reg) for this trait
      for (lf in loss_f[[j]]) {
        rec <- rec + lf(out, prep$target_cat[[j]], prep$mask_cat[[j]])
      }
    }
    kl   <- -0.5 * torch_mean(1 + res$logvar -
                                res$mu$pow(2) -
                                torch_exp(res$logvar))
    loss <- rec + kl
    
    opt$zero_grad(); loss$backward(); opt$step()
    hist[ep] <- loss$item()
    
    if (hist[ep] < best) {
      best <- hist[ep]; wait <- 0
      torch_save(list(
        enc_phy = enc_phy$state_dict(),
        enc_env = enc_env$state_dict(),
        heads   = heads$state_dict()
      ), file.path(ckpt_path, "best.pt"))
    } else if ((wait <- wait + 1) >= patience) {
      message("Early stop @ epoch ", ep)
      break
    }
    
    if (ep %% 100 == 0)
      cat(sprintf("Epoch %4d | Loss: %.4f\n", ep, hist[ep]))
  }
  
  # reload best
  ckpt <- torch_load(file.path(ckpt_path, "best.pt"))
  enc_phy$load_state_dict(ckpt$enc_phy)
  enc_env$load_state_dict(ckpt$enc_env)
  heads$load_state_dict(ckpt$heads)
  lapply(list(enc_phy, enc_env, heads), function(m) m$eval())
  
  # Monte Carlo posterior
  probs_cat <- setNames(
    lapply(cat_cols, function(nm) {
      K <- meta$n_level[meta$name == nm]
      array(0, dim = c(nrow(trait_df), K))
    }),
    cat_cols
  )
  
  for (s in seq_len(n_samples)) {
    res <- enc_phy(prep$X_sp, prep$A)
    z   <- (res$mu + torch_exp(0.5*res$logvar)*torch_randn_like(res$mu))$
      index_select(dim = 0, index = prep$sp_idx) +
      enc_env(prep$env)
    
    for (j in seq_along(heads)) {
      nm  <- cat_cols[j]
      out <- heads[[j]](z)
      pr  <- switch(meta$type[j],
                    binary  = torch_cat(list(1 - torch_sigmoid(out),
                                             torch_sigmoid(out)), dim = 2),
                    nominal = nnf_softmax(out, dim = 2),
                    { # ordinal → pk
                      cum  <- out
                      ones <- torch_ones(cum$size(1), 1, device = device)
                      zeros<- torch_zeros_like(ones)
                      Pcum <- torch_cat(list(ones, cum, zeros), dim = 2)
                      Pcum[,1:meta$n_level[j]] - Pcum[,2:(meta$n_level[j]+1)]
                    }
      )
      probs_cat[[nm]] <- probs_cat[[nm]] + as_array(pr$cpu())
    }
  }
  
  # finalize
  completed_cat <- setNames(
    lapply(cat_cols, function(nm) {
      avg <- probs_cat[[nm]] / n_samples
      cls <- apply(avg, 1, which.max)
      factor(cls,
             levels = seq_along(levels_lst[[nm]]),
             labels = levels_lst[[nm]])
    }),
    cat_cols
  )
  
  if (restore_observed) {
    for (nm in cat_cols) {
      obs <- !is.na(trait_df[[nm]])
      completed_cat[[nm]][obs] <- trait_df[[nm]][obs]
    }
  }
  
  list(
    completed_cat = as.data.frame(completed_cat),
    probs_cat     = probs_cat,
    history       = hist[seq_len(ep)],
    meta          = meta
  )
}


# ─────────────────────────────────────────────────────────────
# Execution Example
# ─────────────────────────────────────────────────────────────

# Convert to factors if needed
traits_df_miss <- as.data.frame(lapply(traits_df_miss, function(x) {
  if (!is.factor(x)) factor(x) else x
}))

# Run imputation
cat_fit <- impute_phylo_cat(
  trait_data = traits_df_miss[, !grepl("^cnt", names(traits_df_miss))],
  phylo_tree = tree,
  env_data = env_data,
  species_id = tree$tip.label,
  latent_dim = 128,
  hidden_dim = 256,
  depth = 2,
  epochs = 2000,
  patience = 300,
  n_samples = 200
)

