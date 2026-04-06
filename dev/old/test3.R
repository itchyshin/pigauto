library(torch)
library(ape)

#— 0. Infer column types like mice -----------------------------------------
#' Infer column types for mixed‐type data (continuous, count, proportion, nominal, ordinal)
#'
#' @param df A data.frame or matrix of trait data.
#' @return A named character vector, one of
#'   \code{"continuous"}, \code{"count"}, \code{"proportion"}, \code{"nominal"}, or \code{"ordinal"}.
#' @examples
#' df <- data.frame(
#'   cont = rnorm(10),
#'   cnt  = rpois(10, 3),
#'   prop = runif(10),
#'   nom  = factor(sample(c("A","B"),10,TRUE)),
#'   ord  = ordered(sample(c("low","high"),10,TRUE), levels=c("low","high"))
#' )
#' infer_col_types(df)
infer_col_types <- function(df) {
  sapply(names(df), function(col) {
    x <- df[[col]]
    if (is.ordered(x)) {
      "ordinal"
    } else if (is.factor(x)) {
      "nominal"
    } else if (is.numeric(x)) {
      if (all(!is.na(x)) && all(x >= 0 & x <= 1)) {
        "proportion"
      } else if (all(floor(x) == x, na.rm = TRUE)) {
        "count"
      } else {
        "continuous"
      }
    } else {
      stop("Cannot infer type for column ", col)
    }
  }, USE.NAMES = TRUE)
}

#— 1. Input Validation -----------------------------------------------------
validate_inputs <- function(trait_data, phylo_tree, env_data, species_id) {
  if (! (is.matrix(trait_data) || is.data.frame(trait_data)))
    stop("trait_data must be a matrix or data.frame.")
  if (! (inherits(phylo_tree, "phylo") || is.matrix(phylo_tree)))
    stop("phylo_tree must be an ape::phylo object or a numeric matrix.")
  if (nrow(env_data) != length(species_id) ||
      nrow(env_data) != nrow(as.matrix(trait_data)))
    stop("Row counts of trait_data, env_data, and species_id must match.")
}

#— 2. Phylogeny → adjacency matrix -----------------------------------------
tree_to_graph <- function(phylo_tree) {
  if (inherits(phylo_tree, "phylo")) {
    D <- cophenetic(phylo_tree)
    A <- 1/(1 + D); diag(A) <- 0
    dimnames(A) <- list(phylo_tree$tip.label, phylo_tree$tip.label)
    A
  } else {
    phylo_tree
  }
}
normalize_adj <- function(A) {
  nm   <- dimnames(A)
  Dinv <- diag(1/sqrt(rowSums(A))); Dinv[!is.finite(Dinv)] <- 0
  B    <- Dinv %*% A %*% Dinv
  if (!is.null(nm)) dimnames(B) <- nm
  B
}

#— 3. Prepare data for torch (handles NA in env_data) ----------------------
#' @inheritParams train_mixed
#' @param trait_data   data.frame of mixed traits
#' @param env_data     numeric data.frame of covariates (may contain NAs)
#' @param species_id   character vector of species labels
#' @param A_hat        normalized adjacency matrix
#' @param col_types    named vector of column types (output of infer_col_types)
#' @param levels_list  named list of factor levels for categorical columns
#' @return A list of torch tensors and metadata for training
prepare_tensors <- function(trait_data, env_data, species_id,
                            A_hat, col_types, levels_list) {
  ## 1) Numeric trait columns
  num_types <- c("continuous","count","proportion")
  num_cols  <- names(col_types)[col_types %in% num_types]
  X_num     <- as.matrix(trait_data[, num_cols, drop = FALSE])
  mask_num  <- !is.na(X_num)
  X_num[!mask_num] <- 0
  X_t       <- torch_tensor(X_num, dtype = torch_float())
  M_num_t   <- torch_tensor(mask_num, dtype = torch_float())
  
  ## 2) Categorical trait columns
  cat_cols   <- names(col_types)[col_types %in% c("nominal","ordinal")]
  target_cat <- mask_cat <- list()
  for (col in cat_cols) {
    fac    <- factor(trait_data[[col]], levels = levels_list[[col]])
    obs    <- !is.na(fac)
    codes  <- as.integer(fac); codes[!obs] <- 1L
    target_cat[[col]] <- torch_tensor(codes, dtype = torch_long())
    mask_cat  [[col]] <- torch_tensor(obs,   dtype = torch_bool())
  }
  
  ## 3) Environmental covariates (now allowing NAs)
  E_mat    <- as.matrix(env_data)
  mask_env <- !is.na(E_mat)
  E_mat[!mask_env] <- 0
  E_t       <- torch_tensor(E_mat,    dtype = torch_float())
  M_env_t   <- torch_tensor(mask_env, dtype = torch_float())
  
  ## 4) Species‐level indexing
  sp_f   <- factor(species_id, levels = unique(species_id))
  sp_idx <- torch_tensor(as.integer(sp_f), dtype = torch_long())
  sp_lev <- levels(sp_f)
  
  ## 5) Species‐level trait means (numeric only)
  Xsp <- t(sapply(sp_lev, function(s) {
    rows <- which(species_id == s)
    colMeans(X_num[rows,, drop = FALSE], na.rm = TRUE)
  }))
  m_sp       <- !is.na(Xsp)
  Xsp[!m_sp] <- 0
  X_sp_t     <- torch_tensor(Xsp, dtype = torch_float())
  M_sp_t     <- torch_tensor(m_sp, dtype = torch_float())
  
  ## 6) Subset phylogenetic adjacency matrix to species
  idx  <- match(sp_lev, rownames(A_hat))
  if (any(is.na(idx))) idx <- seq_along(sp_lev)
  A_sub <- A_hat[idx, idx, drop = FALSE]
  A_t   <- torch_tensor(A_sub, dtype = torch_float())
  
  list(
    X_num       = X_t,
    M_num       = M_num_t,
    target_cat  = target_cat,
    mask_cat    = mask_cat,
    env         = E_t,
    M_env       = M_env_t,
    sp_idx      = sp_idx,
    X_sp        = X_sp_t,
    M_sp        = M_sp_t,
    A           = A_t,
    col_types   = col_types,
    levels_list = levels_list
  )
}

#— 4. Decoders for each type -----------------------------------------------
ContinuousDecoder <- nn_module(
  initialize = function(latent_dim) self$fc <- nn_linear(latent_dim,1),
  forward    = function(z)            self$fc(z)
)

# NB‐style CountDecoder with learned dispersion
CountDecoder <- nn_module(
  initialize = function(latent_dim) {
    self$fc_mu    <- nn_linear(latent_dim,1)
    self$fc_phi   <- nn_linear(latent_dim,1)
    self$softplus <- nn_softplus()
  },
  forward = function(z) {
    mu  <- self$softplus(self$fc_mu(z)) + 1e-6
    phi <- self$softplus(self$fc_phi(z)) + 1e-6
    list(mu = mu, phi = phi)
  }
)

ProportionDecoder <- nn_module(
  initialize = function(latent_dim) self$fc <- nn_linear(latent_dim,1),
  forward    = function(z)            torch_sigmoid(self$fc(z))
)

NominalDecoder <- nn_module(
  initialize = function(latent_dim,n_cat) self$fc <- nn_linear(latent_dim,n_cat),
  forward    = function(z)                 self$fc(z)
)

OrdinalDecoder <- nn_module(
  initialize = function(latent_dim,n_cat) {
    self$score     <- nn_linear(latent_dim,1)
    init_cp        <- torch_linspace(-1,1,steps=n_cat-1)
    self$cutpoints <- nn_parameter(init_cp)
  },
  forward = function(z) {
    eta      <- self$score(z)                   # [N,1]
    cp       <- self$cutpoints$unsqueeze(1)     # [1,K-1]
    prob_leq <- torch_sigmoid(cp - eta)         # [N, K-1]
    K1       <- prob_leq$size(2)
    p1   <- prob_leq[,1,drop=FALSE]
    pmid <- prob_leq[,2:K1,drop=FALSE] - prob_leq[,1:(K1-1),drop=FALSE]
    plast<- (1 - prob_leq[,K1,drop=FALSE])
    torch_cat(list(p1, pmid, plast), dim = 2)
  }
)

#— 5. Loss functions per type ----------------------------------------------
masked_mse      <- function(pred,tgt,mask) torch_mean(((pred-tgt)^2)[mask$to(dtype=torch_bool())])
poisson_nll     <- function(lam,tgt,mask) torch_mean(((lam - tgt*torch_log(lam + 1e-6)))[mask$to(dtype=torch_bool())])
bernoulli_nll   <- function(p, tgt,mask) torch_mean(((-tgt*torch_log(p+1e-6) - (1-tgt)*torch_log(1-p+1e-6)))[mask$to(dtype=torch_bool())])
categorical_nll <- function(logits,tgt,mask) {
  logp <- nnf_log_softmax(logits, dim = 2)
  idx  <- tgt$unsqueeze(2)
  lp   <- torch_gather(logp, 2, idx)$squeeze(2)
  torch_mean((-lp)[mask$to(dtype=torch_bool())])
}
ordinal_nll <- function(probs,tgt,mask) {
  pt <- torch_gather(probs, 2, tgt$unsqueeze(2))$squeeze(2)
  torch_mean((-torch_log(pt + 1e-6))[mask$to(dtype=torch_bool())])
}

#— 6. Mixed‐type VGAE trainer ----------------------------------------------
#' Train the mixed‐type VGAE model
#'
#' @param prep        output of \code{prepare_tensors()}
#' @param latent_dim  integer latent embedding dimension
#' @param epochs      integer number of training epochs
#' @param nsamp       integer number of MC draws for inference
#' @return A list with
#'   \item{mean_cont_prop, sd_cont_prop}{matrices for continuous & proportions}
#'   \item{mean_counts, sd_counts}{vectors for counts (NB mean & sd)}
#'   \item{cat_probs, cat_entropy}{named lists for each categorical}
train_mixed <- function(prep, latent_dim, epochs, nsamp) {
  types      <- prep$col_types
  levels_lst <- prep$levels_list
  p          <- length(types)
  
  enc_phy <- VGAE_Encoder(ncol(prep$X_sp), ncol(prep$X_sp)*2, latent_dim)
  enc_env <- EnvEncoder(ncol(prep$env),       ncol(prep$env)*2,       latent_dim)
  
  # build decoders
  decoders <- setNames(vector("list",p), names(types))
  for (col in names(types)) {
    tp <- types[col]
    decoders[[col]] <- switch(tp,
                              continuous = ContinuousDecoder(latent_dim),
                              count      = CountDecoder(latent_dim),
                              proportion = ProportionDecoder(latent_dim),
                              nominal    = NominalDecoder(latent_dim, length(levels_lst[[col]])),
                              ordinal    = OrdinalDecoder(latent_dim, length(levels_lst[[col]]))
    )
  }
  
  all_pars <- c(enc_phy$parameters, enc_env$parameters,
                unlist(lapply(decoders, function(d) d$parameters), recursive = FALSE))
  opt <- optim_adam(all_pars, lr = 1e-2)
  
  for (ep in seq_len(epochs)) {
    enc_phy$train(); enc_env$train(); lapply(decoders, function(d) d$train())
    # encode
    res   <- enc_phy(prep$X_sp, prep$A)
    mu_sp <- res$mu; lv_sp <- res$logvar
    zsp   <- mu_sp + torch_exp(0.5*lv_sp)*torch_randn_like(mu_sp)
    z_obs <- zsp$index_select(1, prep$sp_idx) + enc_env(prep$env)
    # reconstruction loss
    loss_rec <- torch_tensor(0, dtype=torch_float())
    for (j in seq_along(types)) {
      col <- names(types)[j]; tp <- types[col]
      if (tp %in% c("continuous","proportion","count")) {
        mask_j <- prep$M_num[,j]; Xj <- prep$X_num[,j,drop=FALSE]
        out    <- decoders[[col]](z_obs)
        loss_j <- switch(tp,
                         continuous = masked_mse(out, Xj, mask_j),
                         proportion = bernoulli_nll(out, Xj, mask_j),
                         count      = {
                           mu_phi <- out; mu <- mu_phi$mu; phi <- mu_phi$phi
                           # NB NLL: use torch_lgamma, etc.
                           torch_mean(
                             (torch_lgamma(phi + Xj) - torch_lgamma(phi) - torch_lgamma(Xj + 1) +
                                phi*torch_log(phi + 1e-6) + Xj*torch_log(mu + 1e-6) -
                                (phi + Xj)*torch_log(phi + mu + 1e-6)
                             )[mask_j$to(dtype=torch_bool())]
                           )
                         }
        )
      } else {
        mask_j <- prep$mask_cat[[col]]; tgt_j <- prep$target_cat[[col]]
        logits <- decoders[[col]](z_obs)
        loss_j <- if (tp=="nominal")
          categorical_nll(logits, tgt_j, mask_j)
        else
          ordinal_nll(logits, tgt_j, mask_j)
      }
      loss_rec <- loss_rec + loss_j
    }
    # KL
    loss_kl <- -0.5 * torch_mean(1 + lv_sp - mu_sp^2 - torch_exp(lv_sp))
    loss    <- loss_rec + loss_kl
    opt$zero_grad(); loss$backward(); opt$step()
    if (ep %% 50 == 0)
      cat(sprintf("Epoch %3d | rec=%.4f | kl=%.4f\n", ep, loss_rec$item(), loss_kl$item()))
  }
  
  #— 7. Inference (MC sampling) --------------------------------------------
  enc_phy$eval(); enc_env$eval(); lapply(decoders, function(d) d$eval())
  n_obs   <- prep$X_num$size(1)
  n_contpr <- sum(prep$col_types[num_idx] != "count")
  
  cont_prop_out <- array(NA_real_, c(nsamp, n_obs, n_contpr))
  count_out     <- matrix(NA_real_,    nsamp, n_obs)
  cat_cols      <- names(prep$col_types)[prep$col_types %in% c("nominal","ordinal")]
  cat_draws     <- setNames(
    lapply(cat_cols, function(.) matrix(NA_integer_, nsamp, n_obs)),
    cat_cols
  )
  
  for (i in seq_len(nsamp)) {
    r    <- enc_phy(prep$X_sp, prep$A)
    mu_s <- r$mu; lv_s <- r$logvar
    zsp  <- mu_s + torch_exp(0.5*lv_s)*torch_randn_like(mu_s)
    zobs <- zsp$index_select(1, prep$sp_idx) + enc_env(prep$env)
    # continuous & prop
    ci <- 1L
    for (j in seq_along(types)) {
      tp  <- prep$col_types[j]; dec <- decoders[[j]]
      if (tp %in% c("continuous","proportion")) {
        out <- as_array(dec(zobs))[,1]
        cont_prop_out[i,,ci] <- out; ci <- ci + 1L
      } else if (tp == "count") {
        mu_phi <- dec(zobs)
        mu      <- as_array(mu_phi$mu)
        phi     <- as_array(mu_phi$phi)
        # approximate NB variance, then Gaussian draws
        var_ct  <- mu + mu^2/phi
        sd_ct   <- sqrt(var_ct)
        draw    <- rnorm(n_obs, mean = mu, sd = sd_ct)
        draw[draw < 0] <- 0
        count_out[i,] <- draw
      }
    }
    # categoricals
    for (col in cat_cols) {
      tp     <- prep$col_types[col]
      logits <- decoders[[col]](zobs)
      if (tp == "nominal") {
        probs <- as_array(nnf_softmax(logits, dim = 2))
        cat_draws[[col]][i,] <- apply(probs, 1, function(pr) sample.int(length(pr),1,prob = pr))
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
    t(apply(mat, 2, function(x) tabulate(x, nbins = max(x))/length(x)))
  )
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

#— 8. Master wrapper -------------------------------------------------------
#' Impute mixed‐type traits via Phylogenetic VGAE
#'
#' @param trait_data            data.frame of mixed traits (continuous, factor, etc)
#' @param phylo_tree            ape::phylo or adjacency matrix
#' @param env_data              numeric data.frame of covariates (may contain NAs)
#' @param species_id            vector of species names matching rows
#' @param latent_dim            latent dimension (default 16)
#' @param epochs                training epochs (default 300)
#' @param n_samples             MC draws for uncertainty (default 50)
#' @param col_types             optional named vector of types; if NULL, inferred
#' @param mask_obs_uncertainty  logical; zero‐out SD/entropy on observed if TRUE
#' @param restore_observed      logical; restore original values on observed if TRUE
#' @return A list with:
#'   \describe{
#'     \item{mean_cont_prop, sd_cont_prop}{matrices of size n_obs × n_contpr}
#'     \item{mean_counts, sd_counts}{vectors of length n_obs for counts}
#'     \item{cat_probs, cat_entropy}{named lists for each categorical}
#'     \item{original_data, col_types, levels_list,
#'           contpr_names, count_name, cat_names}{metadata for MI}
#'   }
#' @examples
#' \dontrun{
#' set.seed(42)
#' tree       <- rtree(20); sp <- tree$tip.label
#' n          <- 100
#' species_id <- sample(sp, n, TRUE)
#' cont <- rnorm(n, 10,2); cnt <- rpois(n,3); prop <- runif(n)
#' nom  <- factor(sample(LETTERS[1:3],n,TRUE))
#' ord  <- ordered(sample(c("low","med","high"),n,TRUE), levels=c("low","med","high"))
#' traits <- data.frame(cont,cnt,prop,nom,ord)
#' env    <- data.frame(temp=rnorm(n), precip=runif(n))
#' # introduce some NAs
#' traits[sample(n,20), "cont"] <- NA
#' env[sample(n,10), "precip"] <- NA
#' res <- impute_phylo(traits, tree, env, species_id,
#'                     latent_dim=8, epochs=100, n_samples=25)
#' str(res)
#' }
#' @export
impute_phylo <- function(trait_data, phylo_tree, env_data, species_id,
                         latent_dim           = 16,
                         epochs               = 300,
                         n_samples            = 50,
                         col_types            = NULL,
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
    # continuous/proportion
    for (k in seq_along(contpr_names)) {
      nm  <- contpr_names[k]; mu <- out$mean_cont_prop[,k]; sd <- out$sd_cont_prop[,k]
      obs <- !is.na(orig[[nm]])
      if (mask_obs_uncertainty)   sd[obs] <- 0
      if (restore_observed)       mu[obs] <- orig[[nm]][obs]
      out$mean_cont_prop[,k] <- mu; out$sd_cont_prop[,k] <- sd
    }
    # counts
    if (nzchar(count_name)) {
      obs <- !is.na(orig[[count_name]])
      if (mask_obs_uncertainty)   out$sd_counts[obs] <- 0
      if (restore_observed)       out$mean_counts[obs] <- orig[[count_name]][obs]
    }
    # categoricals
    for (col in cat_names) {
      lvl <- levels_list[[col]]; obs <- !is.na(orig[[col]])
      if (mask_obs_uncertainty) {
        onehot <- t(sapply(orig[[col]][obs], function(v){
          vec <- numeric(length(lvl)); vec[match(v,lvl)] <- 1; vec
        }))
        out$cat_probs[[col]][obs,] <- onehot
        out$cat_entropy[[col]][obs] <- 0
      }
    }
  }
  
  # metadata for MI
  out$original_data <- as.data.frame(trait_data)
  out$col_types     <- col_types
  out$levels_list   <- levels_list
  out$contpr_names  <- contpr_names
  out$count_name    <- count_name
  out$cat_names     <- cat_names
  
  out
}

#— 9. Multiple‐Imputation sampler (fixed indexing) -------------------------
#' Sample M complete data‐sets from a VGAE fit + uncertainty
#'
#' @param fit A list returned by \code{impute_phylo()}
#' @param M   Number of imputed data‐sets to create
#' @return A list of \code{M} data.frames, each with no missing values.
#' @export
mi_create <- function(fit, M) {
  required <- c(
    "mean_cont_prop","sd_cont_prop",
    "mean_counts","sd_counts",
    "cat_probs","original_data",
    "col_types","levels_list",
    "contpr_names","count_name","cat_names"
  )
  missing <- setdiff(required, names(fit))
  if (length(missing)) stop("fit is missing: ", paste(missing, collapse=", "))
  
  orig       <- fit$original_data
  n_obs      <- nrow(orig)
  contpr     <- fit$contpr_names   # vector of names
  count_nm   <- fit$count_name
  cat_vars   <- fit$cat_names
  types      <- fit$col_types
  levels_lst <- fit$levels_list
  
  out <- vector("list", M)
  for (m in seq_len(M)) {
    df <- orig
    
    # 1) continuous + proportion: Normal(mu, sd)
    for (k in seq_along(contpr)) {
      nm <- contpr[k]
      mu <- fit$mean_cont_prop[, k]
      sd <- fit$sd_cont_prop[, k]
      x  <- rnorm(n_obs, mu, sd)
      if (types[nm] == "proportion") x <- pmin(pmax(x, 0), 1)
      df[[nm]] <- x
    }
    
    # 2) counts: Negative‐Binomial using mu & sd → phi
    if (nzchar(count_nm)) {
      mu  <- fit$mean_counts
      sd0 <- fit$sd_counts
      # φ = mu^2/(sd^2 - mu), clamp to >0
      phi <- pmax(mu^2 / (sd0^2 - mu + 1e-6), 1e-6)
      df[[count_nm]] <- rnbinom(n_obs, size = phi, mu = mu)
    }
    
    # 3) categoricals: sample from fitted probabilities
    for (col in cat_vars) {
      P     <- fit$cat_probs[[col]]  # n_obs × K
      codes <- apply(P, 1, function(pr) sample.int(length(pr), 1, prob = pr))
      df[[col]] <- factor(codes,
                          levels = seq_along(levels_lst[[col]]),
                          labels = levels_lst[[col]])
    }
    
    # 4) re‐insert any originally observed
    for (col in names(df)) {
      obs <- !is.na(orig[[col]])
      df[[col]][obs] <- orig[[col]][obs]
    }
    
    out[[m]] <- df
  }
  
  out
}

#—  Example usage --------------------------------------------------------
set.seed(42)
tree       <- rtree(50); sp <- tree$tip.label
n          <- 200
species_id <- sample(sp, n, TRUE)

cont <- rnorm(n, 10, 2)
cnt  <- rpois(n,      3)
prop <- runif(n)
nom  <- factor(sample(c("A","B","C"), n, TRUE))
ord  <- ordered(sample(c("low","med","high"), n, TRUE),
                levels = c("low","med","high"))

traits <- data.frame(cont, cnt, prop, nom, ord)
set.seed(99)
for (j in seq_along(traits)) {
  miss_i <- sample(n, size = floor(0.2*n))
  traits[miss_i, j] <- NA
}
env <- data.frame(temp = rnorm(n))

# after sourcing my updated impute_phylo() definition…
res <- impute_phylo(
  trait_data           = traits,
  phylo_tree           = tree,
  env_data             = env,
  species_id           = species_id,
  latent_dim           = 8,
  epochs               = 150,
  n_samples            = 30,
  mask_obs_uncertainty = TRUE,
  restore_observed     = TRUE
)

# check:
names(res)
# should include “cont_prop_draws”, “count_draws”, “cat_draws”, and “original_data”

# Tidy continuous/proportion
num_cols        <- c("cont","prop")
means_cp        <- as.data.frame(res$mean_cont_prop)
colnames(means_cp) <- num_cols
sds_cp          <- as.data.frame(res$sd_cont_prop)
colnames(sds_cp)   <- num_cols

df_num <- data.frame(
  original_cont  = traits$cont[1:6],
  imputed_cont   = means_cp$cont[1:6],
  sd_cont        = sds_cp$cont[1:6],
  original_prop  = traits$prop[1:6],
  imputed_prop   = means_cp$prop[1:6],
  sd_prop        = sds_cp$prop[1:6],
  original_count = traits$cnt[1:6],
  imputed_count  = res$mean_counts[1:6],
  sd_count       = res$sd_counts[1:6]
)
print(df_num)

# Tidy categorical
modal_from <- function(pr, levs) {
  idx <- apply(pr, 1, which.max)
  factor(idx, levels = seq_along(levs), labels = levs)
}

pred_nom <- modal_from(res$cat_probs$nom, levels(traits$nom))
pred_ord <- modal_from(res$cat_probs$ord, levels(traits$ord))

df_nom <- data.frame(
  original = as.character(traits$nom)[1:10],
  imputed  = as.character(pred_nom)[1:10],
  missing  = is.na(traits$nom)[1:10]
)
print(df_nom)

df_ord <- data.frame(
  original = as.character(traits$ord)[1:10],
  imputed  = as.character(pred_ord)[1:10],
  missing  = is.na(traits$ord)[1:10]
)
print(df_ord)

# (Optional) view entropy
data.frame(
  obs         = 1:10,
  entropy_nom = res$cat_entropy$nom[1:10],
  entropy_ord = res$cat_entropy$ord[1:10]
)
