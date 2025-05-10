###############################################################################
## Phylogenetic VGAE imputer – counts ~ Gaussian
###############################################################################
library(torch)
library(ape)

# ────────────────────────────────────────────────────────────────────────────
# 0. Column-type inference ---------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
infer_col_types <- function(df) {
  sapply(names(df), function(col) {
    x <- df[[col]]
    if      (is.ordered(x))                       "ordinal"
    else if (is.factor(x))                        "nominal"
    else if (is.numeric(x)) {
      if (all(!is.na(x)) && all(x >= 0 & x <= 1)) "proportion"
      else if (all(floor(x) == x, na.rm = TRUE))  "count"
      else                                        "continuous"
    } else stop("Cannot infer type for column ", col)
  }, USE.NAMES = TRUE)
}

# ────────────────────────────────────────────────────────────────────────────
# 1. Input checks ------------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
validate_inputs <- function(trait_data, phylo_tree, env_data, species_id) {
  if (! (is.matrix(trait_data) || is.data.frame(trait_data)))
    stop("trait_data must be a matrix or data.frame.")
  if (! (inherits(phylo_tree, "phylo") || is.matrix(phylo_tree)))
    stop("phylo_tree must be an ape::phylo object or a numeric matrix.")
  if (nrow(env_data) != length(species_id) ||
      nrow(env_data) != nrow(as.matrix(trait_data)))
    stop("Row counts of trait_data, env_data, and species_id must match.")
}

# ────────────────────────────────────────────────────────────────────────────
# 2. Phylogeny → adjacency ---------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
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

# ────────────────────────────────────────────────────────────────────────────
# 3. Prepare tensors ---------------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
prepare_tensors <- function(trait_data, env_data, species_id,
                            A_hat, col_types, levels_list)
{
  ## (i) numeric traits ------------------------------------------------------
  num_types <- c("continuous","count","proportion")
  num_cols  <- names(col_types)[col_types %in% num_types]
  X_num     <- as.matrix(trait_data[, num_cols, drop = FALSE])
  mask_num  <- !is.na(X_num)
  X_num[!mask_num] <- 0
  X_t     <- torch_tensor(X_num, dtype = torch_float())
  M_num_t <- torch_tensor(mask_num, dtype = torch_bool())
  
  ## (ii) categoricals -------------------------------------------------------
  cat_cols   <- names(col_types)[col_types %in% c("nominal","ordinal")]
  target_cat <- mask_cat <- list()
  for (col in cat_cols) {
    fac   <- factor(trait_data[[col]], levels = levels_list[[col]])
    obs   <- !is.na(fac)
    codes <- as.integer(fac); codes[!obs] <- 1L
    target_cat[[col]] <- torch_tensor(codes, dtype = torch_long())
    mask_cat  [[col]] <- torch_tensor(obs,   dtype = torch_bool())
  }
  
  ## (iii) environment -------------------------------------------------------
  E_mat    <- as.matrix(env_data)
  mask_env <- !is.na(E_mat)
  E_mat[!mask_env] <- 0
  E_t     <- torch_tensor(E_mat,    dtype = torch_float())
  M_env_t <- torch_tensor(mask_env, dtype = torch_bool())
  
  ## (iv) species indexing ---------------------------------------------------
  sp_f   <- factor(species_id, levels = unique(species_id))
  sp_idx1 <- as.integer(sp_f)                    # R 1-based
  sp_idx0 <- sp_idx1 - 1L                        # torch expects 0-based
  sp_idx0 <- torch_tensor(sp_idx0, dtype = torch_long())
  
  sp_lev <- levels(sp_f)
  Xsp <- t(sapply(sp_lev, function(s){
    rows <- which(species_id == s)
    colMeans(X_num[rows,, drop = FALSE], na.rm = TRUE)
  }))
  m_sp <- !is.na(Xsp); Xsp[!m_sp] <- 0
  X_sp_t <- torch_tensor(Xsp, dtype = torch_float())
  M_sp_t <- torch_tensor(m_sp, dtype = torch_bool())
  
  ## (v) phylo adjacency -----------------------------------------------------
  idx  <- match(sp_lev, rownames(A_hat))
  if (any(is.na(idx))) idx <- seq_along(sp_lev)
  A_sub <- A_hat[idx, idx, drop = FALSE]
  A_t   <- torch_tensor(A_sub, dtype = torch_float())
  
  list(
    X_num       = X_t,  M_num       = M_num_t,
    target_cat  = target_cat, mask_cat = mask_cat,
    env         = E_t,  M_env       = M_env_t,
    sp_idx0     = sp_idx0,
    X_sp        = X_sp_t, M_sp  = M_sp_t,
    A           = A_t,
    col_types   = col_types,
    levels_list = levels_list,
    num_cols    = num_cols
  )
}

# ────────────────────────────────────────────────────────────────────────────
# 4. Decoder modules (unchanged) ---------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
ContinuousDecoder <- nn_module(
  initialize = function(latent_dim) self$fc <- nn_linear(latent_dim, 1),
  forward    = function(z)            self$fc(z)
)
ProportionDecoder <- nn_module(
  initialize = function(latent_dim) self$fc <- nn_linear(latent_dim, 1),
  forward    = function(z)            torch_sigmoid(self$fc(z))
)
NominalDecoder <- nn_module(
  initialize = function(latent_dim, n_cat) self$fc <- nn_linear(latent_dim, n_cat),
  forward    = function(z)                  self$fc(z)
)
OrdinalDecoder <- nn_module(
  initialize = function(latent_dim, n_cat) {
    self$score     <- nn_linear(latent_dim, 1)
    init_cp        <- torch_linspace(-1, 1, steps = n_cat - 1)
    self$cutpoints <- nn_parameter(init_cp)
  },
  forward = function(z) {
    eta  <- self$score(z)
    cp   <- self$cutpoints$unsqueeze(1)
    leq  <- torch_sigmoid(cp - eta)
    K1   <- leq$size(2)
    torch_cat(list(
      leq[,1,drop=FALSE],
      leq[,2:K1,drop=FALSE] - leq[,1:(K1-1),drop=FALSE],
      1 - leq[,K1,drop=FALSE]
    ), dim = 2)
  }
)

# ────────────────────────────────────────────────────────────────────────────
# 5. Loss helpers (unchanged) ------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
masked_mse <- function(pred, tgt, mask)
  torch_mean(((pred - tgt)^2)[mask])
bernoulli_nll <- function(p, tgt, mask)
  torch_mean(((-tgt*torch_log(p + 1e-6) -
                 (1 - tgt)*torch_log(1 - p + 1e-6)))[mask])
categorical_nll <- function(logits, tgt, mask) {
  lp <- nnf_log_softmax(logits, dim = 2)$gather(2, tgt$unsqueeze(2))$squeeze(2)
  torch_mean((-lp)[mask])
}
ordinal_nll <- function(probs, tgt, mask) {
  pt <- probs$gather(2, tgt$unsqueeze(2))$squeeze(2)
  torch_mean((-torch_log(pt + 1e-6))[mask])
}

# ────────────────────────────────────────────────────────────────────────────
# 6. Trainer (uses sp_idx0 everywhere) ---------------------------------------
# ────────────────────────────────────────────────────────────────────────────
train_mixed <- function(prep, latent_dim, epochs, nsamp) {
  types      <- prep$col_types
  levels_lst <- prep$levels_list
  num_cols   <- prep$num_cols
  num_types  <- c("continuous","count","proportion")
  
  ## encoders
  enc_phy <- VGAE_Encoder(ncol(prep$X_sp), ncol(prep$X_sp)*2, latent_dim)
  enc_env <- EnvEncoder(ncol(prep$env), ncol(prep$env)*2, latent_dim)
  
  ## decoders
  decoders <- setNames(vector("list", length(types)), names(types))
  for (col in names(types)) {
    tp <- types[col]
    decoders[[col]] <- switch(tp,
                              continuous = ContinuousDecoder(latent_dim),
                              count      = ContinuousDecoder(latent_dim),  # Gaussian
                              proportion = ProportionDecoder(latent_dim),
                              nominal    = NominalDecoder(latent_dim, length(levels_lst[[col]])),
                              ordinal    = OrdinalDecoder(latent_dim, length(levels_lst[[col]]))
    )
  }
  
  ## optimizer
  all_pars <- c(enc_phy$parameters, enc_env$parameters,
                unlist(lapply(decoders, function(d) d$parameters), recursive = FALSE))
  opt <- optim_adam(all_pars, lr = 1e-2)
  
  X_num    <- prep$X_num
  M_num    <- prep$M_num
  sp_idx0  <- prep$sp_idx0
  
  for (ep in seq_len(epochs)) {
    enc_phy$train(); enc_env$train(); lapply(decoders, function(d) d$train())
    
    r     <- enc_phy(prep$X_sp, prep$A)
    mu_sp <- r$mu; lv_sp <- r$logvar
    zsp   <- mu_sp + torch_exp(0.5 * lv_sp) * torch_randn_like(mu_sp)
    zobs  <- zsp$index_select(1, sp_idx0) + enc_env(prep$env)
    
    loss_rec <- torch_tensor(0, dtype = torch_float())
    cont_i <- 1L
    for (col in names(types)) {
      tp <- types[col]
      if (tp %in% num_types) {
        maskj <- M_num[,cont_i]; xj <- X_num[,cont_i,drop=FALSE]
        out   <- decoders[[col]](zobs)
        loss_rec <- loss_rec + switch(tp,
                                      continuous = masked_mse(out, xj, maskj),
                                      count      = masked_mse(out, xj, maskj),
                                      proportion = bernoulli_nll(out, xj, maskj)
        )
        cont_i <- cont_i + 1L
      } else {
        maskj  <- prep$mask_cat[[col]]
        tgtj   <- prep$target_cat[[col]]
        logits <- decoders[[col]](zobs)
        loss_rec <- loss_rec + if (tp=="nominal")
          categorical_nll(logits, tgtj, maskj) else
            ordinal_nll(logits, tgtj, maskj)
      }
    }
    
    loss_kl <- -0.5 * torch_mean(1 + lv_sp - mu_sp^2 - torch_exp(lv_sp))
    loss    <- loss_rec + loss_kl
    opt$zero_grad(); loss$backward(); opt$step()
    
    if (ep %% 50 == 0)
      cat(sprintf("Epoch %3d | recon=%.4f | KL=%.4f\n",
                  ep, loss_rec$item(), loss_kl$item()))
  }
  
  ## inference (MC sampling)
  enc_phy$eval(); enc_env$eval(); lapply(decoders, function(d) d$eval())
  
  n_obs    <- X_num$size(1)
  contpr_cols <- num_cols[types[num_cols] != "count"]
  n_contpr <- length(contpr_cols)
  nsamp    <- as.integer(nsamp)
  
  cont_prop_out <- array(NA_real_, c(nsamp, n_obs, n_contpr))
  count_out     <- matrix(NA_real_, nsamp, n_obs)
  cat_cols      <- names(types)[types %in% c("nominal","ordinal")]
  cat_draws     <- setNames(
    lapply(cat_cols, function(.) matrix(NA_integer_, nsamp, n_obs)),
    cat_cols
  )
  
  for (i in seq_len(nsamp)) {
    r    <- enc_phy(prep$X_sp, prep$A)
    zsp  <- r$mu + torch_exp(0.5*r$logvar) * torch_randn_like(r$mu)
    zobs <- zsp$index_select(1, sp_idx0) + enc_env(prep$env)
    
    # numeric
    cont_i <- 1L
    for (col in num_cols) {
      tp <- types[col]
      out <- as_array(decoders[[col]](zobs))[,1]
      if (tp == "count") {
        count_out[i,] <- pmax(rnorm(n_obs, out, sd = 1), 0)
      } else {
        if (tp == "proportion") out <- pmin(pmax(out, 0), 1)
        cont_prop_out[i,,cont_i] <- out
        cont_i <- cont_i + 1L
      }
    }
    
    # categoricals
    for (col in cat_cols) {
      logits <- decoders[[col]](zobs)
      if (types[col]=="nominal") {
        probs <- as_array(nnf_softmax(logits, dim = 2))
        cat_draws[[col]][i,] <- apply(probs, 1,
                                      function(pr) sample.int(length(pr), 1, prob = pr))
      } else {
        probs <- as_array(logits)
        cat_draws[[col]][i,] <- apply(probs, 1, function(pr){
          cs <- cumsum(pr); u <- runif(1); which(u <= cs)[1]
        })
      }
    }
  }
  
  mean_cp  <- apply(cont_prop_out, 2:3, mean)
  sd_cp    <- apply(cont_prop_out, 2:3, sd)
  mean_cnt <- apply(count_out, 2, mean)
  sd_cnt   <- apply(count_out, 2, sd)
  cat_probs   <- lapply(cat_draws, function(mat)
    t(apply(mat, 2, function(x) tabulate(x, nbins = max(x))/length(x))))
  cat_entropy <- lapply(cat_probs, function(m) -rowSums(m * log(m + 1e-8)))
  
  list(
    mean_cont_prop = mean_cp,
    sd_cont_prop   = sd_cp,
    mean_counts    = mean_cnt,
    sd_counts      = sd_cnt,
    cat_probs      = cat_probs,
    cat_entropy    = cat_entropy
  )
}

# ────────────────────────────────────────────────────────────────────────────
# 7. Wrapper + MI sampler remain unchanged (just call train_mixed()).
# ────────────────────────────────────────────────────────────────────────────
impute_phylo <- function(trait_data, phylo_tree, env_data, species_id,
                         latent_dim = 16, epochs = 300, n_samples = 50,
                         col_types  = NULL,
                         mask_obs_uncertainty = TRUE,
                         restore_observed     = TRUE) {
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  if (is.null(col_types)) col_types <- infer_col_types(trait_data)
  levels_list   <- lapply(trait_data, function(x) if (is.factor(x)) levels(x) else NULL)
  contpr_names  <- names(col_types)[col_types %in% c("continuous","proportion")]
  count_name    <- names(col_types)[col_types == "count"]
  cat_names     <- names(col_types)[col_types %in% c("nominal","ordinal")]
  
  A_hat <- normalize_adj(tree_to_graph(phylo_tree))
  prep  <- prepare_tensors(trait_data, env_data, species_id,
                           A_hat, col_types, levels_list)
  out   <- train_mixed(prep, latent_dim, epochs, n_samples)
  
  if (mask_obs_uncertainty || restore_observed) {
    orig <- as.data.frame(trait_data)
    # continuous / proportion
    for (k in seq_along(contpr_names)) {
      nm  <- contpr_names[k]
      obs <- !is.na(orig[[nm]])
      if (mask_obs_uncertainty) out$sd_cont_prop[obs,k] <- 0
      if (restore_observed)     out$mean_cont_prop[obs,k] <- orig[[nm]][obs]
    }
    # counts
    if (nzchar(count_name)) {
      obs <- !is.na(orig[[count_name]])
      if (mask_obs_uncertainty) out$sd_counts[obs]   <- 0
      if (restore_observed)     out$mean_counts[obs] <- orig[[count_name]][obs]
    }
    # categoricals
    for (col in cat_names) {
      obs <- !is.na(orig[[col]])
      if (mask_obs_uncertainty) {
        lvl <- levels_list[[col]]
        onehot <- t(sapply(orig[[col]][obs], function(v){
          z <- numeric(length(lvl)); z[match(v,lvl)] <- 1; z
        }))
        out$cat_probs[[col]][obs,]  <- onehot
        out$cat_entropy[[col]][obs] <- 0
      }
    }
  }
  
  out$original_data <- as.data.frame(trait_data)
  out$col_types     <- col_types
  out$levels_list   <- levels_list
  out$contpr_names  <- contpr_names
  out$count_name    <- count_name
  out$cat_names     <- cat_names
  class(out) <- "phylo_vgae_impute"
  out
}

# ────────────────────────────────────────────────────────────────────────────
# 8. Multiple‑imputation sampler -------------------------------------------
# ────────────────────────────────────────────────────────────────────────────
mi_create <- function(fit, M = 5) {
  orig     <- fit$original_data
  n_obs    <- nrow(orig)
  out_lst  <- vector("list", M)
  
  for (m in seq_len(M)) {
    df <- orig
    
    ## continuous / proportion
    for (k in seq_along(fit$contpr_names)) {
      nm <- fit$contpr_names[k]
      x  <- rnorm(n_obs, fit$mean_cont_prop[,k], fit$sd_cont_prop[,k])
      if (fit$col_types[nm] == "proportion") x <- pmin(pmax(x, 0), 1)
      df[[nm]] <- x
    }
    
    ## counts  (Gaussian draw → round & clamp)
    if (nzchar(fit$count_name)) {
      mu <- fit$mean_counts; sd <- fit$sd_counts
      x  <- round(rnorm(n_obs, mu, sd)); x[x < 0] <- 0
      df[[fit$count_name]] <- x
    }
    
    ## categoricals
    for (col in fit$cat_names) {
      P <- fit$cat_probs[[col]]
      codes <- apply(P, 1, function(pr) sample.int(length(pr), 1, prob = pr))
      df[[col]] <- factor(codes,
                          levels = seq_along(fit$levels_list[[col]]),
                          labels = fit$levels_list[[col]])
    }
    
    ## put observed back
    for (col in names(df)) {
      obs <- !is.na(orig[[col]])
      df[[col]][obs] <- orig[[col]][obs]
    }
    out_lst[[m]] <- df
  }
  out_lst
}

# Example usage:
set.seed(42)
tree       <- rtree(30); sp <- tree$tip.label
n          <- 150
species_id <- sample(sp, n, TRUE)
traits <- data.frame(
  cont = rnorm(n, 10, 2),
  cnt  = rpois(n, 3),
  prop = runif(n),
  nom  = factor(sample(c("A","B","C"), n, TRUE)),
  ord  = ordered(sample(c("low","med","high"), n, TRUE),
                 levels = c("low","med","high"))
)
for (j in seq_along(traits))
  traits[sample.int(n, floor(0.2*n)), j] <- NA
env <- data.frame(temp = rnorm(n), precip = runif(n))
res <- impute_phylo(traits, tree, env, species_id,
                    latent_dim = 8, epochs = 150, n_samples = 30)


