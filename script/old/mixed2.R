library(torch)
library(ape)
library(MASS)
library(Matrix)
library(Metrics)
library(abind)

# 0) device -------------------------------------------------------------
device <- if (cuda_is_available()) torch_device("cuda") else
  if (backends_mps_is_available()) torch_device("mps") else
    torch_device("cpu")

device <- torch_device("cpu") # for testing

# 1) helpers ------------------------------------------------------------
validate_inputs <- function(trait_data, phylo_tree, env_data, species_id) {
  stopifnot(is.data.frame(trait_data) || is.matrix(trait_data))
  stopifnot(inherits(phylo_tree, "phylo") || is.matrix(phylo_tree))
  stopifnot(nrow(env_data) == length(species_id))
  stopifnot(nrow(trait_data) == length(species_id))
}
tree_to_graph   <- function(tr) if (inherits(tr,"phylo")) { d<-cophenetic(tr); A<-1/(1+d); diag(A)<-0; dimnames(A)<-list(tr$tip.label,tr$tip.label); A } else tr
normalize_adj   <- function(A) { D<-diag(1/sqrt(rowSums(A))); D[!is.finite(D)]<-0; B<-D%*%A%*%D; dimnames(B)<-dimnames(A); B }

infer_trait_types <- function(df) {
  n <- ncol(df)
  out_type   <- character(n)
  out_levels <- integer(n)
  
  for (j in seq_len(n)) {
    x <- df[[j]]
    if (is.numeric(x) && all(na.omit(x) %in% 0:1)) {
      out_type[j]   <- "binary"
      out_levels[j] <- 2
    } else if (is.numeric(x)) {
      out_type[j]   <- "cont"
      out_levels[j] <- NA_integer_
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
  
  data.frame(
    name    = colnames(df),
    type    = out_type,
    n_level = out_levels,
    stringsAsFactors = FALSE
  )
}

# 2) build stacks -------------------------------------------------------
linear_stack <- function(in_dim, hid_dim, depth, final_dim=NULL) {
  layers <- list(nn_linear(in_dim, hid_dim), nn_relu())
  if (depth>1) for (i in 2:depth) layers <- c(layers,list(nn_linear(hid_dim,hid_dim),nn_relu()))
  if (!is.null(final_dim)) layers <- c(layers,list(nn_linear(hid_dim,final_dim)))
  nn_sequential(!!!layers)
}

# 3) modules ------------------------------------------------------------
VGAE_Encoder <- nn_module(
  "VGAE_Encoder",
  initialize=function(inp,hid,lat,depth){
    self$enc<-linear_stack(inp,hid,depth)
    self$mu <-nn_linear(hid,lat)
    self$lv <-nn_linear(hid,lat)
  },
  forward=function(Xsp,A){
    H<-self$enc(torch_matmul(A,Xsp))
    list(mu=self$mu(H), logvar=self$lv(H))
  }
)
EnvEncoder <- nn_module(
  "EnvEncoder",
  initialize=function(idim,hid,lat,depth) self$enc<-linear_stack(idim,hid,depth,lat),
  forward=function(e) self$enc(e)
)
ContinuousHead <- nn_module(
  "ContinuousHead",
  initialize=function(lat) self$lin<-nn_linear(lat,1),
  forward=function(z) self$lin(z)
)
BinaryHead <- nn_module(
  "BinaryHead",
  initialize=function(lat) self$lin<-nn_linear(lat,1),
  forward=function(z) self$lin(z)
)
NominalHead <- nn_module(
  "NominalHead",
  initialize=function(lat,K) self$lin<-nn_linear(lat,K),
  forward=function(z) self$lin(z)
)
OrdinalHead <- nn_module(
  "OrdinalHead",
  initialize=function(lat,K){
    self$beta<-nn_linear(lat,1,bias=FALSE)
    cp<-torch_linspace(-1,1,steps=K-1,device=device)
    self$cut<-nn_parameter(cp)
  },
  forward=function(z){
    eta<-self$beta(z)          # (n×1)
    torch_sigmoid(self$cut-eta) # (n×(K-1))
  }
)

ordinal_nll <- function(cumprob,tgt,mask,K){
  n <- cumprob$size(1)
  ones  <- torch_ones(n,1,device=device)
  zeros <- torch_zeros_like(ones)
  Pcum  <- torch_cat(list(ones,cumprob,zeros),dim=2)           # n×(K+1)
  pk    <- Pcum[,1:K] - Pcum[,2:(K+1)]                         # n×K
  idx   <- tgt$unsqueeze(2)$to(dtype=torch_long()) + 1L        # 1-based
  p_obs <- pk$gather(2,idx)$squeeze(2)
  p_obs <- torch_clamp(p_obs,1e-7,1-1e-7)
  nll <- -p_obs$log() * mask$to(dtype=torch_float())
  torch_mean(nll)
}

# -------------------------------------------------------------------------
# 4) prepare tensors  (renamed Xsp_c -> Xsp so every function sees
#                      exactly the same slot name)
# -------------------------------------------------------------------------
prepare_tensors <- function(df, env, species_id, Ahat, meta) {
  
  cont_idx <- meta$type == "cont"
  
  ## ----- continuous ------------------------------------------------------
  Xc  <- as.matrix(df[, cont_idx, drop = FALSE])
  Mc  <- !is.na(Xc); Xc[!Mc] <- 0
  Xc_t  <- torch_tensor(Xc, dtype = torch_float(), device = device)
  Mc_t  <- torch_tensor(Mc, dtype = torch_bool(),  device = device)
  
  ## ----- categorical -----------------------------------------------------
  cat_idx <- !cont_idx
  if (any(cat_idx)) {
    mats  <- masks <- list(); cnt <- 0
    for (j in which(cat_idx)) {
      cnt  <- cnt + 1
      x    <- df[[j]]
      code <- as.integer(x)          # 1…K ; NA handled below
      mask <- !is.na(code)
      code[!mask] <- 1L              # dummy level
      mats[[cnt]]  <- code
      masks[[cnt]] <- mask
    }
    Xcat  <- do.call(cbind, mats);  storage.mode(Xcat) <- "integer"
    Mcat  <- do.call(cbind, masks)
  } else {
    Xcat  <- matrix(1L, nrow(df), 1)
    Mcat  <- matrix(FALSE, nrow(df), 1)
  }
  Xcat_t <- torch_tensor(Xcat, dtype = torch_long(),  device = device)
  Mcat_t <- torch_tensor(Mcat, dtype = torch_bool(),  device = device)
  
  ## ----- environmental ---------------------------------------------------
  E_t <- torch_tensor(as.matrix(env), dtype = torch_float(), device = device)
  
  ## ----- phylogeny -------------------------------------------------------
  sp_f   <- factor(species_id, levels = unique(species_id))
  sp_idx <- torch_tensor(as.integer(sp_f), dtype = torch_long(), device = device)
  
  Xsp <- t(sapply(levels(sp_f), function(s) {
    m <- which(species_id == s)
    colMeans(Xc[m, , drop = FALSE], na.rm = TRUE)
  }))
  Xsp[is.na(Xsp)] <- 0
  Xsp_t <- torch_tensor(Xsp, dtype = torch_float(), device = device)
  
  A_t   <- torch_tensor(Ahat[levels(sp_f), levels(sp_f)],
                        dtype = torch_float(), device = device)
  
  list(
    Xc   = Xc_t,   Mc   = Mc_t,
    Xcat = Xcat_t, Mcat = Mcat_t,
    env  = E_t,
    Xsp  = Xsp_t,  A    = A_t,
    sp_idx = sp_idx
  )
}

# -------------------------------------------------------------------------
# 5) training loop  (now references prep$Xsp, not the old Xsp_c) ----------
# -------------------------------------------------------------------------
train_vgae <- function(prep, meta, depth, lat,
                       epochs, lr, patience, ckpt) {
  
  p_cont <- sum(meta$type == "cont")
  J      <- nrow(meta)
  hid    <- max(2 * lat, 16)
  
  enc_phy <- VGAE_Encoder(p_cont, hid, lat, depth)$to(device)
  enc_env <- EnvEncoder(prep$env$size(2), hid, lat, depth)$to(device)
  
  ## ---- heads & per-trait loss fns ---------------------------------------
  heads  <- nn_module_list()
  loss_f <- vector("list", J)
  
  ci <- 0
  for (j in seq_len(J)) {
    tp <- meta$type[j]
    if (tp == "cont") {
      ci <- ci + 1
      heads$append(ContinuousHead(lat)$to(device))
      loss_f[[j]] <- function(pred, tgt, mask)
        torch_mean((pred - tgt)$pow(2) *
                     mask$to(dtype = torch_float()))
    } else if (tp == "binary") {
      heads$append(BinaryHead(lat)$to(device))
      loss_f[[j]] <- function(logit, tgt, mask) {
        l <- nnf_binary_cross_entropy_with_logits(
          logit,
          tgt$to(dtype = torch_float())$unsqueeze(2),
          reduction = "none")
        torch_mean(l * mask$to(dtype = torch_float())$unsqueeze(2))
      }
    } else if (tp == "nominal") {
      K <- meta$n_level[j]
      heads$append(NominalHead(lat, K)$to(device))
      loss_f[[j]] <- function(logit, tgt, mask) {
        l <- nnf_cross_entropy(logit,
                               tgt$to(dtype = torch_long()),
                               reduction = "none")
        torch_mean(l * mask$to(dtype = torch_float()))
      }
    } else {                              # ordinal
      K <- meta$n_level[j]
      heads$append(OrdinalHead(lat, K)$to(device))
      loss_f[[j]] <- function(cum, tgt, mask)
        ordinal_nll(cum, tgt, mask, K)
    }
  }
  
  opt <- optim_adam(
    c(enc_phy$parameters,
      enc_env$parameters,
      heads$parameters),
    lr = lr)
  
  best  <- Inf; wait <- 0; hist <- numeric(epochs); stop_ep <- epochs
  if (!is.null(ckpt)) {
    dir.create(ckpt, recursive = TRUE, showWarnings = FALSE)
    ckf <- file.path(ckpt, "vgae.pt")
  }
  
  for (ep in seq_len(epochs)) {
    enc_phy$train(); enc_env$train()
    
    Rp  <- enc_phy(prep$Xsp, prep$A)
    Re  <- enc_env(prep$env)
    eps <- torch_randn_like(Rp$mu)
    zsp <- Rp$mu + torch_exp(0.5 * Rp$logvar) * eps
    z   <- zsp$index_select(1, prep$sp_idx) + Re
    
    ## reconstruction loss -------------------------------------------------
    rec <- torch_tensor(0, device = device)
    ci  <- 0; cci <- 0
    for (j in seq_len(J)) {
      tp <- meta$type[j]
      if (tp == "cont") {
        ci <- ci + 1
        pr <- heads[[j]](z)[,1]
        rec <- rec + loss_f[[j]](pr, prep$Xc[, ci], prep$Mc[, ci])
      } else {
        cci  <- cci + 1
        out  <- heads[[j]](z)
        rec  <- rec + loss_f[[j]](
          out, prep$Xcat[, cci], prep$Mcat[, cci])
      }
    }
    
    kl   <- -0.5 * torch_mean(
      1 + Rp$logvar - Rp$mu$pow(2) - torch_exp(Rp$logvar))
    loss <- rec + kl
    
    opt$zero_grad(); loss$backward(); opt$step()
    
    cur <- loss$item(); hist[ep] <- cur
    if (cur < best && is.finite(cur)) {
      best <- cur; wait <- 0; stop_ep <- ep
      if (!is.null(ckpt))
        torch_save(list(enc_phy = enc_phy$state_dict(),
                        enc_env = enc_env$state_dict(),
                        heads   = heads$state_dict()), ckf)
    } else {
      wait <- wait + 1
      if (wait >= patience) {
        message("Early stop @ epoch ", ep); break
      }
    }
    if (ep %% 200 == 0)
      cat(sprintf("Epoch %4d | loss = %.4f\n", ep, cur))
  }
  
  if (!is.null(ckpt)) {
    st <- torch_load(ckf)
    enc_phy$load_state_dict(st$enc_phy)
    enc_env$load_state_dict(st$enc_env)
    heads$load_state_dict(st$heads)
  }
  
  list(enc_phy = enc_phy,
       enc_env = enc_env,
       heads   = heads,
       history = hist[1:stop_ep],
       stopped = stop_ep)
}

# -------------------------------------------------------------------------
# >>>  NO other part of the script needs to change  <<<
#     (impute_phylo() already refers to prep$Xsp)
# -------------------------------------------------------------------------

# function

impute_phylo <- function(
    trait_data, phylo_tree, env_data, species_id,
    latent_dim   = 32, hidden_depth = 1,
    epochs       = 2000, lr = 1e-3, patience = 100,
    n_samples    = 200, ckpt_path = "ckpt",
    mask_obs_uncertainty = TRUE, restore_observed = TRUE
) {
  ## 1) prep + train (unchanged) ---------------------------------------------
  validate_inputs(trait_data, phylo_tree, env_data, species_id)
  meta <- infer_trait_types(as.data.frame(trait_data))
  
  cont_idx <- meta$type == "cont"
  # z-scale continuous
  Xsc <- as.data.frame(trait_data)
  if (any(cont_idx)) {
    mu    <- colMeans(Xsc[ , cont_idx, drop=FALSE], na.rm=TRUE)
    sigma <- apply(   Xsc[ , cont_idx, drop=FALSE], 2, sd, na.rm=TRUE)
    Xsc[ , cont_idx] <- sweep(
      sweep(Xsc[ , cont_idx, drop=FALSE], 2, mu, "-"),
      2, sigma, "/"
    )
  } else mu <- sigma <- NULL
  
  A_hat <- normalize_adj(tree_to_graph(phylo_tree))
  prep  <- prepare_tensors(Xsc, env_data, species_id, A_hat, meta)
  fit   <- train_vgae(
    prep, meta,
    depth    = hidden_depth,
    lat      = latent_dim,
    epochs   = epochs,
    lr       = lr,
    patience = patience,
    ckpt     = ckpt_path
  )
  enc_phy <- fit$enc_phy; enc_env <- fit$enc_env; heads <- fit$heads
  enc_phy$eval(); enc_env$eval(); heads$eval()
  
  ## 2) set up MC‐draw containers --------------------------------------------
  n <- nrow(trait_data); J <- ncol(trait_data)
  
  # for continuous
  n_cont <- sum(cont_idx)
  sum_cont <- matrix(0, n, n_cont)
  sumsq_cont <- matrix(0, n, n_cont)
  
  # for categorical: one list element per non‐cont trait
  cat_names <- meta$name[!cont_idx]
  draws_cat_counts <- setNames(
    lapply(cat_names, function(nm){
      K <- meta$n_level[ meta$name == nm ]
      matrix(0, n, K)
    }),
    cat_names
  )
  
  ## 3) Monte Carlo sampling -------------------------------------------------
  for (s in seq_len(n_samples)) {
    Rphy <- enc_phy(prep$Xsp_c, prep$A)
    Renv <- enc_env(prep$env)
    eps  <- torch_randn_like(Rphy$mu)
    zsp  <- Rphy$mu + torch_exp(0.5*Rphy$logvar) * eps
    zobs <- zsp$index_select(1, prep$sp_idx) + Renv
    
    # continuous part
    ci <- 0L
    for (j in seq_len(J)) if (cont_idx[j]) {
      ci <- ci + 1L
      out <- as_array(heads[[j]](zobs)[ ,1])
      sum_cont[ ,ci]   <- sum_cont[ ,ci]   + out
      sumsq_cont[ ,ci] <- sumsq_cont[ ,ci] + out^2
    }
    
    # categorical part
    cci <- 0L
    for (j in seq_len(J)) if (!cont_idx[j]) {
      cci <- cci + 1L
      nm <- cat_names[cci]
      logits <- heads[[j]](zobs)
      K <- meta$n_level[j]
      if (meta$type[j] == "binary") {
        # build two-column [P(0), P(1)] from the single logit
          p1    <- torch_sigmoid(logits)                   # n×1
          probs <- torch_cat(list(1 - p1, p1), dim = 2)    # n×2
          } else if (meta$type[j] == "nominal") {
          probs <- nnf_softmax(logits, dim = 2)            # n×K
            } else {
               #ordinal → reconstruct pk from cumulative
                cum  <- logits                                   # n×(K-1)
                ones <- torch_ones(n,1,device=device)
                zeros<- torch_zeros_like(ones)
                Pcum <- torch_cat(list(ones,cum,zeros), dim=2)   # n×(K+1)
                probs<- Pcum[,1:K] - Pcum[,2:(K+1)]              # n×K
                }
      pr_mat <- as_array(probs$cpu())
      draws_cat_counts[[nm]] <- draws_cat_counts[[nm]] + pr_mat
    }
  }
  
  ## 4) finalize continuous --------------------------------------------------
  completed_cont <- sum_cont / n_samples
  var_cont       <- (sumsq_cont / n_samples) - completed_cont^2
  uncertainty_cont <- sqrt(pmax(var_cont, 0))
  
  if (!is.null(mu)) {
    completed_cont    <- sweep(sweep(completed_cont,    2, sigma, "*"), 2, mu, "+")
    uncertainty_cont  <- sweep(uncertainty_cont, 2, sigma, "*")
  }
  colnames(completed_cont)   <- meta$name[cont_idx]
  colnames(uncertainty_cont) <- meta$name[cont_idx]
  
  ## 5) finalize categorical outputs -----------------------------------------
  # completed_cat: pick the modal (argmax) class from the **average** probs
  completed_cat <- lapply(cat_names, function(nm){
    avg_p <- draws_cat_counts[[nm]] / n_samples
    cls   <- apply(avg_p, 1, which.max)
    factor(cls,
           levels=seq_along(levels(trait_data[[nm]])),
           labels=levels(trait_data[[nm]]))
  })
  completed_cat <- as.data.frame(setNames(completed_cat, cat_names))
  
  # uncertainty_cat: percent in each level
  uncertainty_cat <- lapply(cat_names, function(nm){
    mat <- draws_cat_counts[[nm]] / n_samples * 100
    colnames(mat) <- levels(trait_data[[nm]])
    as.data.frame(mat)
  })
  names(uncertainty_cat) <- cat_names
  
  ## 6) restore / mask observed if requested ---------------------------------
  # (continuous)
  orig_cont_mat <- as.matrix(trait_data)[ , cont_idx, drop = FALSE]
  obs_m         <- !is.na(orig_cont_mat)
  
  if (mask_obs_uncertainty)
    uncertainty_cont[obs_m] <- 0
  
  if (restore_observed)
    completed_cont[obs_m]   <- orig_cont_mat[obs_m]
  
  # (categorical: if you want to force original obs → you could overwrite
  #  completed_cat and set uncertainty_cat to zero rows for those obs)
  
  ## 7) return ----------------------------------------------------------------
  list(
    completed_cont   = as.data.frame(completed_cont),
    uncertainty_cont = as.data.frame(uncertainty_cont),
    completed_cat    = completed_cat,
    uncertainty_cat  = uncertainty_cat,
    history          = fit$history,
    stopped_epoch    = fit$stopped,
    meta             = meta
  )
}

# test

# data generation 

###############################################################################
## Simulate phylogenetically-correlated continuous *and* categorical traits
## plus environmental covariates, then impose 25% MCAR missingness
###############################################################################

library(ape)        # rtree()
library(MASS)       # mvrnorm()
library(Matrix)     # nearPD()

set.seed(123)

## 1. SETTINGS
n_sp   <- 300          # number of species
p_cts  <-  5           # continuous traits
p_bin  <-  3           # binary traits
p_ord  <-  1           # ordinal traits
p_mult <-  1           # multinomial traits (3-level)
p_env  <-  6           # environmental predictors

## 2. PHYLOGENY
tree   <- rtree(n_sp)
V_phy  <- cov2cor(vcv(tree))

## 3. ENVIRONMENT
Sigma_env <- diag(p_env)
env_data  <- mvrnorm(n_sp, rep(0, p_env), Sigma_env)
colnames(env_data) <- paste0("env", 1:p_env)

## 4A. CONTINUOUS TRAITS
lambda <- runif(p_cts, 0.2, 0.8)
Sigma_list <- lapply(lambda, function(lam) lam * V_phy + (1 - lam) * diag(n_sp))
Z_phy <- sapply(Sigma_list, function(Sig) mvrnorm(n = 1, mu = rep(0, n_sp), Sigma = Sig))
B_cts <- matrix(runif(p_cts * p_env, -1, 1), p_cts, p_env)
set.seed(42)
cts_traits <- Z_phy +
  as.matrix(env_data) %*% t(B_cts) +
  matrix(rnorm(n_sp * p_cts, 0, 0.1), n_sp, p_cts)
colnames(cts_traits) <- paste0("cnt", 1:p_cts)

## 4B. BINARY TRAITS
lambda_bin <- runif(p_bin, 0.2, 0.8)
Sigma_bin_list <- lapply(lambda_bin, function(lam) lam * V_phy + (1 - lam) * diag(n_sp))
Z_phy_bin <- sapply(Sigma_bin_list, function(Sig) mvrnorm(n = 1, mu = rep(0, n_sp), Sigma = Sig))
B_bin <- matrix(runif(p_bin * p_env, -1, 1), p_bin, p_env)
linpred_bin <- Z_phy_bin + as.matrix(env_data) %*% t(B_bin)
prob_bin <- plogis(linpred_bin)
bin_traits <- matrix(rbinom(n = n_sp * p_bin, size = 1, prob = as.vector(prob_bin)), nrow = n_sp, ncol = p_bin)
colnames(bin_traits) <- paste0("bin", 1:p_bin)

## 4C. ORDINAL TRAITS
lambda_ord <- runif(p_ord, 0.2, 0.8)
Sigma_ord_list <- lapply(lambda_ord, function(lam) lam * V_phy + (1 - lam) * diag(n_sp))
Z_phy_ord <- sapply(Sigma_ord_list, function(Sig) mvrnorm(n = 1, mu = rep(0, n_sp), Sigma = Sig))
B_ord <- matrix(runif(p_ord * p_env, -1, 1), p_ord, p_env)
set.seed(42)
resid_ord <- matrix(rnorm(n_sp * p_ord, 0, 0.1), n_sp, p_ord)
linpred_ord <- Z_phy_ord + as.matrix(env_data) %*% t(B_ord) + resid_ord
ord_traits <- apply(linpred_ord, 2, function(latent) {
  thr <- quantile(latent, probs = c(1/3, 2/3))
  cut(
    latent,
    breaks = c(-Inf, thr, Inf),
    labels = c("low", "mid", "high"),
    ordered_result = TRUE
  )
})

ord_traits <- factor(ord_traits, levels = c("low", "mid", "high"), labels = c("low", "mid", "high"), ordered = TRUE)
ord_traits <- as.data.frame(ord_traits)
colnames(ord_traits) <- paste0("ord", seq_len(p_ord))

## 4D. MULTINOMIAL TRAIT (3-level)
n_mult   <- p_mult
n_levels <- 3
phylo_cor <- 0.5
lambda_mult <- c(0.5, 0.7)
V1 <- lambda_mult[1] * V_phy + (1 - lambda_mult[1]) * diag(n_sp)
V2 <- lambda_mult[2] * V_phy + (1 - lambda_mult[2]) * diag(n_sp)
Cov12 <- phylo_cor * sqrt(lambda_mult[1] * lambda_mult[2]) * V_phy
top    <- cbind(V1, Cov12)
bottom <- cbind(t(Cov12), V2)
Sigma_mult <- rbind(top, bottom)
Z_mult_latent <- matrix(mvrnorm(1, mu = rep(0, 2 * n_sp), Sigma = Sigma_mult), nrow = n_sp, ncol = 2)
B_mult <- matrix(runif((n_levels - 1) * p_env, -1, 1), n_levels - 1, p_env)
linpred_mult <- Z_mult_latent + as.matrix(env_data) %*% t(B_mult) + matrix(rnorm(n_sp * (n_levels - 1), 0, 0.1), n_sp, n_levels - 1)
softmax <- function(x) exp(x) / sum(exp(x))
logit_mat <- cbind(0, linpred_mult)
probs <- t(apply(logit_mat, 1, softmax))
mult_trait <- apply(probs, 1, function(p) sample(1:3, 1, prob = p))
mult_traits <- factor(mult_trait, levels = 1:3, labels = c("A", "B", "C"))
mult_traits <- as.data.frame(mult_traits)
colnames(mult_traits) <- paste0("mult", seq_len(n_mult))

## 5. COMBINE TRAITS INTO DATA FRAME
traits_df <- data.frame(
  cts_traits,
  bin_traits,
  ord_traits,
  mult_traits
)

## 6. IMPOSING 25% MCAR MISSINGNESS — SUPPORTS ALL TYPES
miss_mat <- matrix(runif(n_sp * ncol(traits_df)) < 0.25, n_sp, ncol(traits_df))
traits_df_miss <- traits_df # copy

for (j in seq_along(traits_df)) {
  idx <- which(miss_mat[,j])
  if (is.numeric(traits_df[[j]]) || is.integer(traits_df[[j]])) {
    traits_df_miss[idx, j] <- NA
  } else if (is.factor(traits_df[[j]])) {
    # factors, incl. ordered and multinomial
    traits_df_miss[[j]][idx] <- NA
  }
}

## 7. READY FOR VGAE PIPELINE OR EXPLORATION
env_df    <- as.data.frame(env_data)

str(traits_df_miss[1:5, ])
str(env_df[1:5, ])

# Now traits_df_miss contains 25% MCAR missing data for all trait types!

#' Classify trait columns by variable type
#'
#' @param df A data.frame containing trait columns.
#' @param binary_threshold Integer: max unique values for binary (default: 2)
#' @return A list with elements: continuous, binary, ordinal, multinomial
#' @examples
#' classify_traits(traits_df_miss)
classify_traits <- function(df, binary_threshold = 2) {
  res <- list(continuous = character(0),
              binary     = character(0),
              ordinal    = character(0),
              multinomial= character(0))
  
  for (col in names(df)) {
    x <- df[[col]]
    # Ignore completely missing columns
    if (all(is.na(x))) next
    
    # Ordinal: ordered factor
    if (is.ordered(x)) {
      res$ordinal <- c(res$ordinal, col)
    } 
    # Multinomial: unordered factor (not binary)
    else if (is.factor(x)) {
      n_lev <- nlevels(x)
      if (n_lev > binary_threshold) {
        res$multinomial <- c(res$multinomial, col)
      } else {
        res$binary <- c(res$binary, col) # fallback: factor with 2 levels = binary
      }
    }
    # Binary: integer/numeric, only 2 unique non-missing values (e.g., 0,1)
    else if (is.integer(x) || is.numeric(x)) {
      uniq <- unique(x[!is.na(x)])
      if (length(uniq) == binary_threshold &&
          all(uniq %in% c(0, 1))) {
        res$binary <- c(res$binary, col)
      } else {
        res$continuous <- c(res$continuous, col)
      }
    } 
    # Fallback: treat as continuous
    else {
      res$continuous <- c(res$continuous, col)
    }
  }
  res
}

#classify_traits(traits_df_miss)


#source("phylo_vgae_impute.R")
res <- impute_phylo(
  trait_data   = as.data.frame(traits_df_miss),
  phylo_tree   = tree,
  env_data     = env_data,
  species_id   = tree$tip.label,
  latent_dim   = 128,
  hidden_depth = 1,
  epochs       = 2000,
  lr           = 1e-2,
  patience     = 300,
  n_samples    = 200,
  ckpt_path    = "checkpoints"
)

# visualise


