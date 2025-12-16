###############################################################################
## Phylogenetic VGAE imputer
##   • Per-trait z-scaling
##   • Early-stopping + checkpoint
##   • MC-dropout inference
##   • Mixed continuous / binary / nominal / ordinal support
###############################################################################

library(torch)
library(ape)
library(MASS)       # mvrnorm()
library(Matrix)     # nearPD()
library(Metrics)    # rmse()
library(ggplot2)
library(abind)      # only used for quick tensor padding

# ════════════════════════════════════════════════════════════════════════════
# 0)  Detect device -----------------------------------------------------------
# ════════════════════════════════════════════════════════════════════════════
device <- if (cuda_is_available()) {
  torch_device("cuda")
} else if (backends_mps_is_available()) {
  torch_device("mps")
} else {
  torch_device("cpu")
}

# ════════════════════════════════════════════════════════════════════════════
# 1)  Helpers + checks --------------------------------------------------------
# ════════════════════════════════════════════════════════════════════════════
validate_inputs <- function(trait_data, phylo_tree, env_data, species_id) {
  if (! (is.data.frame(trait_data) || is.matrix(trait_data)))
    stop("`trait_data` must be matrix or data.frame.")
  if (! (inherits(phylo_tree, "phylo") || is.matrix(phylo_tree)))
    stop("`phylo_tree` must be an ape::phylo or a square adjacency matrix.")
  if (nrow(env_data) != length(species_id))
    stop("`env_data` rows must equal length of `species_id`.")
  if (nrow(trait_data) != length(species_id))
    stop("`trait_data` rows must equal length of `species_id`.")
}

tree_to_graph <- function(tr) {
  if (inherits(tr, "phylo")) {
    d <- cophenetic(tr)
    A <- 1/(1 + d); diag(A) <- 0
    dimnames(A) <- list(tr$tip.label, tr$tip.label)
    A
  } else {
    tr
  }
}

normalize_adj <- function(A) {
  D <- diag(1/sqrt(rowSums(A))); D[!is.finite(D)] <- 0
  B <- D %*% A %*% D
  dimnames(B) <- dimnames(A)
  B
}

# ════════════════════════════════════════════════════════════════════════════
# 2)  Infer column types (cont / binary / nominal / ordinal) ------------------
# ════════════════════════════════════════════════════════════════════════════
infer_trait_types <- function(df) {
  out_type   <- character(ncol(df))
  out_levels <- integer  (ncol(df))
  
  for (j in seq_along(out_type)) {
    x <- df[[j]]
    if (is.numeric(x)) {
      if (all(na.omit(x) %in% 0:1)) {
        out_type[j]   <- "binary"
        out_levels[j] <- 2
      } else {
        out_type[j] <- "cont"
      }
    } else if (is.factor(x)) {
      if (is.ordered(x)) {
        out_type[j]   <- "ordinal"
        out_levels[j] <- nlevels(x)
      } else {
        out_type[j]   <- "nominal"
        out_levels[j] <- nlevels(x)
      }
    } else {
      stop("Unsupported column type in `trait_data`.")
    }
  }
  data.frame(name = colnames(df),
             type = out_type,
             n_level = ifelse(out_type == "cont", NA, out_levels),
             stringsAsFactors = FALSE)
}

# ════════════════════════════════════════════════════════════════════════════
# 3)  Prepare tensors ---------------------------------------------------------
# ════════════════════════════════════════════════════════════════════════════
prepare_tensors <- function(df, env, species_id, A_hat, meta) {
  cont_idx <- meta$type == "cont"
  
  ## continuous --------------------------------------------------------------
  Xcont <- as.matrix(df[ , cont_idx, drop = FALSE])
  Mcont <- !is.na(Xcont)
  Xcont[!Mcont] <- 0
  
  Xc_t <- torch_tensor(Xcont, dtype = torch_float(), device = device)
  Mc_t <- torch_tensor(Mcont, dtype = torch_bool(),  device = device)
  
  ## categorical -------------------------------------------------------------
  cat_idx <- !cont_idx
  if (any(cat_idx)) {
    mats  <- list()
    masks <- list()
    for (j in which(cat_idx)) {
      x   <- df[[j]]
      tp  <- meta$type[j]
      # Always code 1…K for factors; binary numeric will be 0/1 but we leave that
      if (tp == "binary") {
        v_int <- as.integer(x)             #  0 or 1
      } else {
        v_int <- as.integer(x)             #  1…K
      }
      m_cat <- !is.na(v_int)
      # fill missing with a valid class (will be masked)
      v_int[!m_cat] <- 1L
      mats[[length(mats)+1]]  <- v_int
      masks[[length(masks)+1]] <- m_cat
    }
    Xcat_int <- do.call(cbind, mats)
    storage.mode(Xcat_int) <- "integer"
    Mcat     <- do.call(cbind, masks)
    
    Xcat_t <- torch_tensor(Xcat_int, dtype = torch_long(), device = device)
    Mcat_t <- torch_tensor(Mcat,      dtype = torch_bool(),  device = device)
  } else {
    Xcat_t <- torch_tensor(matrix(0L, nrow(df), 1), dtype = torch_long(), device = device)
    Mcat_t <- torch_tensor(matrix(FALSE, nrow(df), 1), dtype = torch_bool(), device = device)
  }
  
  ## environmental -----------------------------------------------------------
  E_t <- torch_tensor(as.matrix(env), dtype = torch_float(), device = device)
  
  ## species-level -----------------------------------------------------------
  sp_f   <- factor(species_id, levels = unique(species_id))
  sp_idx <- torch_tensor(as.integer(sp_f), dtype = torch_long(), device = device)
  sp_lev <- levels(sp_f)
  
  Xsp_c <- t(sapply(sp_lev, function(s) {
    rows <- which(species_id == s)
    colMeans(Xcont[rows, , drop = FALSE], na.rm = TRUE)
  }))
  Xsp_c[is.na(Xsp_c)] <- 0
  Xsp_t <- torch_tensor(Xsp_c, dtype = torch_float(), device = device)
  
  A_sub <- A_hat[sp_lev, sp_lev, drop = FALSE]
  A_t   <- torch_tensor(A_sub, dtype = torch_float(), device = device)
  
  list(Xc     = Xc_t,
       Mc     = Mc_t,
       Xcat   = Xcat_t,
       Mcat   = Mcat_t,
       env    = E_t,
       A      = A_t,
       sp_idx = sp_idx,
       Xsp_c  = Xsp_t)
}

# ════════════════════════════════════════════════════════════════════════════
# 4)  Tiny helper to build MLP stacks ----------------------------------------
# ════════════════════════════════════════════════════════════════════════════
linear_stack <- function(in_dim, hid_dim, depth, final_dim = NULL, last_relu = FALSE) {
  L <- list(nn_linear(in_dim, hid_dim), nn_relu())
  if (depth > 1) {
    for (i in 2:depth)
      L <- c(L, list(nn_linear(hid_dim, hid_dim), nn_relu()))
  }
  if (!is.null(final_dim))
    L <- c(L, list(nn_linear(hid_dim, final_dim)))
  if (isTRUE(last_relu))
    L <- c(L, list(nn_relu()))
  do.call(nn_sequential, L)
}

# ════════════════════════════════════════════════════════════════════════════
# 5)  Encoder / decoder heads -------------------------------------------------
# ════════════════════════════════════════════════════════════════════════════
VGAE_Encoder <- nn_module(
  "VGAE_Encoder",
  initialize = function(inp, hid, lat, depth) {
    self$mlp <- linear_stack(inp, hid, depth)
    self$mu  <- nn_linear(hid, lat)
    self$lv  <- nn_linear(hid, lat)
  },
  forward = function(Xsp, A) {
    H <- self$mlp(torch_matmul(A, Xsp))
    list(mu = self$mu(H), logvar = self$lv(H))
  }
)

EnvEncoder <- nn_module(
  "EnvEncoder",
  initialize = function(idim, hid, lat, depth) {
    self$mlp <- linear_stack(idim, hid, depth, lat)
  },
  forward = function(e) self$mlp(e)
)

ContinuousHead <- nn_module(
  "ContinuousHead",
  initialize = function(lat) self$lin <- nn_linear(lat, 1),
  forward    = function(z) self$lin(z)
)

BinaryHead <- nn_module(
  "BinaryHead",
  initialize = function(lat) self$lin <- nn_linear(lat, 1),
  forward    = function(z) self$lin(z)          # logits
)

NominalHead <- nn_module(
  "NominalHead",
  initialize = function(lat, K) self$lin <- nn_linear(lat, K),
  forward    = function(z) self$lin(z)          # logits
)

OrdinalHead <- nn_module(
  "OrdinalHead",
  initialize = function(lat, K) {
    self$beta <- nn_linear(lat, 1, bias = FALSE)
    # initialize K−1 cutpoints evenly spaced between −1 and +1
    cp0 <- torch_linspace(-1, 1, steps = K-1, device = device)
    self$cut <- nn_parameter(cp0)
  },
  forward = function(z) {
    eta <- self$beta(z)        # n × 1
    torch_sigmoid(self$cut - eta)  # n × (K-1)
  }
)

# ════════════════════════════════════════════════════════════════════════════
# 5a)  Ordinal negative log-likelihood (torch-native) ------------------------
# ════════════════════════════════════════════════════════════════════════════
ordinal_nll <- function(cumprob, tgt, mask, K) {
  n <- cumprob$size(1)
  ones  <- torch_ones(n, 1, device = cumprob$device)
  zeros <- torch_zeros_like(ones)
  Pcum  <- torch_cat(list(ones, cumprob, zeros), dim = 2)  # n × (K+1)
  pk    <- Pcum[ , 1:K] - Pcum[ , 2:(K+1)]                 # n × K
  
  idx   <- tgt$unsqueeze(2)$to(dtype = torch_long()) + 1L  # 1-based
  p_obs <- pk$gather(2, idx)$squeeze(2)                    # n
  
  # clamp away from 0/1 to avoid log(0) → NaN
  p_obs <- torch_clamp(p_obs, min = 1e-7, max = 1 - 1e-7)
  
  nll   <- -(p_obs$log()) * mask$to(dtype = torch_float())
  torch_mean(nll)
}

# ════════════════════════════════════════════════════════════════════════════
# 6)  Training loop with early-stopping --------------------------------------
# ════════════════════════════════════════════════════════════════════════════
train_vgae <- function(prep, meta, depth, lat, epochs, lr, patience, ckpt) {
  p_cont <- sum(meta$type == "cont")
  J      <- nrow(meta)
  
  hid     <- max(2*lat, 16)
  enc_phy <- VGAE_Encoder(p_cont, hid, lat, depth)$to(device = device)
  enc_env <- EnvEncoder(prep$env$size(2), hid, lat, depth)$to(device = device)
  
  # ←─── this must be here
  heads   <- nn_module_list()
  loss_f  <- vector("list", J)
  cont_i  <- 0L
  
  for (j in seq_len(J)) {
    tp <- meta$type[j]
    
    if (tp == "cont") {
      cont_i <- cont_i + 1L
      heads$append(ContinuousHead(lat)$to(device = device))
      loss_f[[j]] <- function(pred, tgt, mask)
        torch_mean((pred - tgt)$pow(2) * mask$to(dtype = torch_float()))
      
    } else if (tp == "binary") {
      heads$append(BinaryHead(lat)$to(device = device))
      loss_f[[j]] <- function(logit, tgt, mask) {
        # make target 300×1 instead of 300
        tgt2  <- tgt$to(dtype = torch_float())$unsqueeze(2)
        # same for mask
        mask2 <- mask$to(dtype = torch_float())$unsqueeze(2)
        # BCE with matching shapes
        l <- nnf_binary_cross_entropy_with_logits(
          logit, tgt2, reduction = "none"
        )
        torch_mean(l * mask2)
      }
      
    } else if (tp == "nominal") {
      K <- meta$n_level[j]
      heads$append(NominalHead(lat, K)$to(device = device))
      loss_f[[j]] <- function(logit, tgt, mask) {
        l <- nnf_cross_entropy(
          logit, tgt$to(dtype = torch_long()), reduction = "none")
        torch_mean(l * mask$to(dtype = torch_float()))
      }
      
    } else {  # ordinal
      K <- meta$n_level[j]
      heads$append(OrdinalHead(lat, K)$to(device = device))
      loss_f[[j]] <- function(cprob, tgt, mask)
        ordinal_nll(cprob, tgt, mask, K)
    }
  }
  
  
  optim <- optim_adam(c(enc_phy$parameters,
                        enc_env$parameters,
                        heads$parameters), lr = lr)
  
  best <- Inf; wait <- 0; hist <- numeric(epochs); stop_ep <- epochs
  
  if (!is.null(ckpt)) {
    dir.create(ckpt, showWarnings = FALSE, recursive = TRUE)
    ckpt_file <- file.path(ckpt, "best_vgae.pt")
  }
  
  for (ep in seq_len(epochs)) {
    enc_phy$train(); enc_env$train()
    
    Rphy <- enc_phy(prep$Xsp_c, prep$A)           # species-level encoder
    Renv <- enc_env(prep$env)                     # env encoder
    eps  <- torch_randn_like(Rphy$mu)
    
    z_sp  <- Rphy$mu + torch_exp(0.5*Rphy$logvar) * eps
    z_obs <- z_sp$index_select(1, prep$sp_idx) + Renv
    
    ## reconstruction loss ---------------------------------------------------
    rec_loss <- torch_tensor(0, dtype = torch_float(), device = device)
   
    rec_loss <- rec_loss + loss_f[[j]](
      out,
      prep$Xcat[, cat_idx]$unsqueeze(2),
      prep$Mcat[, cat_idx]$unsqueeze(2)
    )
    
    ## KL divergence ---------------------------------------------------------
        # After computing `loss`
    cur <- loss$item()
    hist[ep] <- cur

    # Only update the “best” checkpoint if `cur` is a finite number and improves on previous best
    if (!is.na(cur) && is.finite(cur) && cur < best) {
      best   <- cur
      wait   <- 0
      stop_ep<- ep

      # Save model state if checkpointing is enabled
      if (!is.null(ckpt)) {
        torch_save(
          list(
            enc_phy = enc_phy$state_dict(),
            enc_env = enc_env$state_dict(),
            heads   = heads$state_dict()
          ),
          ckpt_file
        )
      }

    } else {
      # Otherwise increment patience counter
      wait <- wait + 1

      # Trigger early stop if we’ve gone too long without improvement
      if (wait >= patience) {
        message(
          "Early stop @ epoch ", ep,
          " (best loss = ", sprintf("%.4f", best), ")"
        )
        break
      }
    }

    # Optional logging every 200 epochs
    if (ep %% 200 == 0) {
      cat(sprintf("Epoch %4d  |  loss = %.4f\n", ep, cur))
    }
  } 
  ## reload best weights -----------------------------------------------------
  if (!is.null(ckpt)) {
    st <- torch_load(ckpt_file)
    enc_phy$load_state_dict(st$enc_phy)
    enc_env$load_state_dict(st$enc_env)
    heads$load_state_dict(st$heads)
  }
  
  list(enc_phy = enc_phy,
       enc_env = enc_env,
       heads   = heads,
       history = hist[1:stop_ep],
       stopped = stop_ep,
       meta    = meta)

}
# ════════════════════════════════════════════════════════════════════════════
# 7)  Main wrapper -----------------------------------------------------------
# ════════════════════════════════════════════════════════════════════════════
impute_phylo <- function(
    trait_data,
    phylo_tree,
    env_data,
    species_id,
    latent_dim    = 32,
    hidden_depth  = 1,
    epochs        = 2000,
    lr            = 1e-3,
    patience      = 100,
    n_samples     = 200,
    ckpt_path     = "ckpt",
    mask_obs_uncertainty = TRUE,
    restore_observed     = TRUE) {
  
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  meta <- infer_trait_types(as.data.frame(trait_data))
  
  ## z-scale continuous traits ----------------------------------------------
  cont_idx <- meta$type == "cont"
  Xsc      <- as.data.frame(trait_data)
  if (any(cont_idx)) {
    mu    <- colMeans(Xsc[ , cont_idx, drop = FALSE], na.rm = TRUE)
    sigma <- apply   (Xsc[ , cont_idx, drop = FALSE], 2, sd, na.rm = TRUE)
    Xsc[ , cont_idx] <- sweep(sweep(Xsc[ , cont_idx, drop = FALSE], 2, mu, "-"),
                              2, sigma, "/")
  } else {
    mu <- sigma <- NULL
  }
  
  ## tensors + training ------------------------------------------------------
  A_hat <- normalize_adj(tree_to_graph(phylo_tree))
  prep  <- prepare_tensors(Xsc, env_data, species_id, A_hat, meta)
  
  fit <- train_vgae(prep, meta,
                    depth    = hidden_depth,
                    lat      = latent_dim,
                    epochs   = epochs,
                    lr       = lr,
                    patience = patience,
                    ckpt     = ckpt_path)
  
  enc_phy <- fit$enc_phy; enc_env <- fit$enc_env; heads <- fit$heads
  enc_phy$eval(); enc_env$eval(); heads$eval()
  
  ## Monte-Carlo draws -------------------------------------------------------
  n <- nrow(trait_data); J <- ncol(trait_data)
  draws <- vector("list", n_samples)
  
  for (s in seq_len(n_samples)) {
    Rphy <- enc_phy(prep$Xsp_c, prep$A)
    Renv <- enc_env(prep$env)
    zsp  <- Rphy$mu + torch_exp(0.5*Rphy$logvar) * torch_randn_like(Rphy$mu)
    zobs <- zsp$index_select(1, prep$sp_idx) + Renv
    
    smpl <- vector("list", J)
    cont_i <- 0L
    
    for (j in seq_len(J)) {
      tp <- meta$type[j]
      hd <- heads[[j]]
      
      if (tp == "cont") {
        cont_i <- cont_i + 1L
        smpl[[j]] <- as_array(hd(zobs)[ ,1])             # n-vector
        
      } else if (tp == "binary") {
        smpl[[j]] <- as_array(torch_sigmoid(hd(zobs))$view(-1))
        
      } else if (tp == "nominal") {
        cls <- torch_argmax(hd(zobs), dim = 2)$cpu()
        smpl[[j]] <- as.integer(cls) + 1L                 # 1 … K
        
      } else {  # ordinal
        cprob <- hd(zobs)
        nlev  <- meta$n_level[j]
        ones  <- torch_ones(cprob$size(1), 1, device = device)
        zeros <- torch_zeros_like(ones)
        Pcum  <- torch_cat(list(ones, cprob, zeros), dim = 2)
        pk    <- Pcum[ , 1:nlev] - Pcum[ , 2:(nlev+1)]
        cls   <- torch_argmax(pk, dim = 2)$cpu()
        smpl[[j]] <- as.integer(cls) + 1L                 # 1 … K
      }
    }
    draws[[s]] <- smpl
  }
  
  ## aggregate ---------------------------------------------------------------
  mean_mat <- matrix(NA_real_, n, J)
  sd_mat   <- matrix(NA_real_, n, J)
  for (j in seq_len(J)) {
    stack_j <- sapply(draws, `[[`, j)                     # n × S
    mean_mat[ , j] <- rowMeans(stack_j)
    sd_mat  [ , j] <- apply(stack_j, 1, sd)
  }
  
  ## un-scale continuous back ------------------------------------------------
  if (any(cont_idx)) {
    mean_mat[ , cont_idx] <- sweep(sweep(mean_mat[ , cont_idx, drop = FALSE],
                                         2, sigma, "*"), 2, mu, "+")
    sd_mat  [ , cont_idx] <- sweep(sd_mat  [ , cont_idx, drop = FALSE],
                                   2, sigma, "*")
  }
  
  ## restore observed + zero SD if desired -----------------------------------
  obs <- !is.na(as.matrix(trait_data))
  if (mask_obs_uncertainty) sd_mat[obs] <- 0
  if (restore_observed)     mean_mat[obs] <- as.matrix(trait_data)[obs]
  
  list(completed_data = as.data.frame(mean_mat),
       uncertainty    = as.data.frame(sd_mat),
       history        = fit$history,
       stopped_epoch  = fit$stopped,
       meta           = fit$meta)
}

# ════════════════════════════════════════════════════════════════════════════
# 8)  Minimal example (uses objects from your data-generation script) ---------
# ════════════════════════════════════════════════════════════════════════════
res <- impute_phylo(
  trait_data   = as.data.frame(traits_df_miss ),
  phylo_tree   = tree,
  env_data     = env_data,
  species_id   = tree$tip.label,
  latent_dim   = 128,
  hidden_depth = 1,
  epochs       = 1000,
  lr           = 1e-2,
  patience     = 300,
  n_samples    = 200,
  ckpt_path    = "checkpoints")


