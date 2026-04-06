# ==============================================================================
# Distilled Phylo-DAE + categorical predictors (one-hot)  [BM-RESIDUAL + FAST MODE]
# + Improvements:
#   (1) Bounded positive residual scale: res_scale = sigmoid(res_raw) * res_cap
#   (2) Supervised train-missing uses extra observed-input dropout (prevents leakage)
#   (3) Finer early stopping (eval_every=20) + shorter max_epochs
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(dplyr)
library(here)

# ----------------------------- USER SETTINGS ----------------------------------

n_sub <- 1200
cat_cols_requested <- c("Diet", "Migration")

tree_path_A <- here("avonet","Stage2_Hackett_MCC_no_neg.tre")
csv_path_A  <- here("avonet","AVONET3_BirdTree.csv")
tree_path_B <- "/mnt/data/Stage2_Hackett_MCC_no_neg.tre"
csv_path_B  <- "/mnt/data/AVONET3_BirdTree.csv"

trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

missing_frac   <- 0.25
val_frac_miss  <- 0.25
test_frac_miss <- 0.25

k_eigen    <- 12
sigma_mult <- 0.35

CFG <- list(
  hidden = 192,
  dropout = 0.10,
  lr = 2e-3,
  weight_decay = 1e-5,
  corrupt_p = 0.30,
  
  # losses
  lambda_anchor = 0.05,
  lambda_trainmiss = 0.50,   # << was 1.00; start lower to reduce overfit
  lambda_resid = 0.50,       # << stronger shrink to BM on missing cells
  lambda_distill_obs = 0.00,
  
  # NEW: supervised leakage control + bounded residual scale
  sup_observe_drop = 0.15,   # drop 15% of observed cells from input during supervised step
  res_cap = 0.30,            # max deviation multiplier on residual net
  lambda_scale = 0.00,       # optional extra penalty on res_scale^2 (try 0.01 if needed)
  
  max_epochs = 600,
  eval_every = 20,
  patience = 15,
  clip_norm = 1.0,
  seed = 1
)

# ----------------------------- helpers ----------------------------------------
rmse_vec <- function(truth, pred) sqrt(mean((truth - pred)^2))

rmse_torch <- function(pred, truth, mask_bool) {
  d <- (pred - truth)
  torch_sqrt(torch_mean(d[mask_bool]$pow(2)))
}

make_missing_splits3 <- function(X, missing_frac=0.25, val_frac=0.25, test_frac=0.25, seed=555) {
  stopifnot(val_frac + test_frac < 1)
  set.seed(seed)
  n <- nrow(X); p <- ncol(X)
  all_idx <- seq_len(n * p)
  m <- floor(missing_frac * length(all_idx))
  miss <- sample(all_idx, m)
  
  m_val  <- floor(val_frac  * length(miss))
  m_test <- floor(test_frac * length(miss))
  
  val_idx   <- miss[seq_len(m_val)]
  test_idx  <- miss[(m_val + 1):(m_val + m_test)]
  train_idx <- miss[(m_val + m_test + 1):length(miss)]
  
  list(train_idx=train_idx, val_idx=val_idx, test_idx=test_idx)
}

all_grads_finite <- function(parameters) {
  for (pp in parameters) {
    g <- pp$grad
    if (!is.null(g)) {
      ok <- torch_isfinite(g)$all()$item()
      if (!isTRUE(ok)) return(FALSE)
    }
  }
  TRUE
}

make_onehot <- function(df, cat_cols, max_levels = 40) {
  cat_cols <- intersect(cat_cols, names(df))
  if (length(cat_cols) == 0) stop("No categorical columns found in df after intersect().")
  
  dd <- df[, cat_cols, drop = FALSE]
  
  keep <- sapply(dd, function(x) {
    ux <- unique(x)
    ux <- ux[!is.na(ux) & ux != ""]
    length(ux) <= max_levels
  })
  
  if (!any(keep)) stop("All categorical columns dropped by max_levels filter (too many levels).")
  
  dropped <- names(dd)[!keep]
  if (length(dropped) > 0) {
    message("Dropping high-cardinality categorical cols: ", paste(dropped, collapse=", "))
  }
  
  dd <- dd[, keep, drop = FALSE]
  for (j in seq_along(dd)) {
    x <- as.character(dd[[j]])
    x[is.na(x) | x == ""] <- "Unknown"
    dd[[j]] <- factor(x)
  }
  
  mm <- model.matrix(~ . - 1, data = dd)
  mm <- as.matrix(mm)
  storage.mode(mm) <- "double"
  list(mat = mm, colnames = colnames(mm), kept = names(dd))
}

get_spectral_features <- function(tree, k = 12, sigma_mult = 0.35) {
  D <- cophenetic(tree)
  sigma <- median(D) * sigma_mult
  A <- exp(-(D^2) / (2 * sigma^2))
  diag(A) <- 0
  L <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric = TRUE)
  n <- nrow(A)
  idx <- (n - k):(n - 1)  # smallest non-zero eigenvectors
  eig$vectors[, idx, drop = FALSE]
}

get_stable_adj <- function(tree, sigma_mult = 0.35) {
  D <- cophenetic(tree)
  sigma <- median(D) * sigma_mult
  A <- exp(-(D^2) / (2 * sigma^2))
  diag(A) <- 0
  rs <- rowSums(A); rs[rs == 0] <- 1
  diag(1/rs) %*% A
}

# ----------------------------- device -----------------------------------------
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")
set.seed(CFG$seed); torch_manual_seed(CFG$seed)

# ----------------------------- data loading -----------------------------------
tree_path <- if (file.exists(tree_path_A)) tree_path_A else tree_path_B
csv_path  <- if (file.exists(csv_path_A))  csv_path_A  else csv_path_B

if (!file.exists(tree_path)) stop(paste0("Tree file not found. Tried:\n  ", tree_path_A, "\n  ", tree_path_B, "\ngetwd()=", getwd()))
if (!file.exists(csv_path))  stop(paste0("CSV file not found. Tried:\n  ", csv_path_A,  "\n  ", csv_path_B,  "\ngetwd()=", getwd()))

cat("Reading tree:", tree_path, "\n")
cat("Reading csv :", csv_path,  "\n")

tree <- read.tree(tree_path)
avonet <- read.csv(csv_path)
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]
stopifnot(all(df_final$Species_Key == tree_pruned$tip.label))

cat_cols_found <- intersect(cat_cols_requested, names(avonet))
if (length(cat_cols_found) == 0) {
  message("None of requested categorical columns were found. You asked for: ",
          paste(cat_cols_requested, collapse=", "))
  message("Some available column names are:")
  print(head(names(avonet), 60))
  stop("Stop: choose valid categorical columns (see names(avonet)).")
}
if (length(setdiff(cat_cols_requested, cat_cols_found)) > 0) {
  message("Dropping missing categorical columns: ",
          paste(setdiff(cat_cols_requested, cat_cols_found), collapse=", "))
}

df_cat <- avonet[match(df_final$Species_Key, avonet$Species_Key),
                 c("Species_Key", cat_cols_found), drop = FALSE]
stopifnot(all(df_cat$Species_Key == df_final$Species_Key))

if (!is.null(n_sub) && n_sub < nrow(df_final)) {
  set.seed(CFG$seed)
  keep_species <- sample(df_final$Species_Key, n_sub)
  tree_pruned <- keep.tip(tree_pruned, keep_species)
  df_final <- df_final[match(tree_pruned$tip.label, df_final$Species_Key), ]
  df_cat   <- df_cat[match(df_final$Species_Key, df_cat$Species_Key), ]
  stopifnot(all(df_final$Species_Key == tree_pruned$tip.label))
  stopifnot(all(df_cat$Species_Key == df_final$Species_Key))
  cat(sprintf("FAST MODE: using n_sub=%d tips\n", n_sub))
}

cat("Dataset:", nrow(df_final), "species x", length(trait_cols), "traits\n")

X_raw   <- as.matrix(df_final[, trait_cols])
X_log   <- log(X_raw)
X_truth <- scale(X_log)

n <- nrow(X_truth); p <- ncol(X_truth)

# ----------------------------- missingness split --------------------------------
splits <- make_missing_splits3(
  X_truth,
  missing_frac = missing_frac,
  val_frac = val_frac_miss,
  test_frac = test_frac_miss,
  seed = 555
)

train_idx <- splits$train_idx
val_idx   <- splits$val_idx
test_idx  <- splits$test_idx

mask_input <- matrix(1, nrow=n, ncol=p)
mask_input[c(train_idx, val_idx, test_idx)] <- 0

X_in <- X_truth
X_in[c(train_idx, val_idx, test_idx)] <- NA

X_fill <- X_in
X_fill[is.na(X_fill)] <- 0

t_X     <- torch_tensor(X_fill, dtype=torch_float(), device=device)
t_mask  <- torch_tensor(mask_input, dtype=torch_float(), device=device)
t_truth <- torch_tensor(X_truth, dtype=torch_float(), device=device)

train_mask_mat <- matrix(FALSE, nrow=n, ncol=p); train_mask_mat[train_idx] <- TRUE
val_mask_mat   <- matrix(FALSE, nrow=n, ncol=p); val_mask_mat[val_idx]     <- TRUE
test_mask_mat  <- matrix(FALSE, nrow=n, ncol=p); test_mask_mat[test_idx]   <- TRUE

t_train_mask <- torch_tensor(train_mask_mat, dtype=torch_bool(), device=device)
t_val_mask   <- torch_tensor(val_mask_mat,   dtype=torch_bool(), device=device)
t_test_mask  <- torch_tensor(test_mask_mat,  dtype=torch_bool(), device=device)

# ----------------------------- baseline: Rphylopars (BM) ----------------------
cat("\n--- Running Rphylopars (BM) ---\n")
df_rphylo <- data.frame(species = df_final$Species_Key, X_in)
res_rphylo <- phylopars(trait_data=df_rphylo, tree=tree_pruned, model="BM", pheno_error=FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]

rmse_train_rphylo <- rmse_vec(X_truth[train_idx], X_pred_rphylo[train_idx])
rmse_val_rphylo   <- rmse_vec(X_truth[val_idx],   X_pred_rphylo[val_idx])
rmse_test_rphylo  <- rmse_vec(X_truth[test_idx],  X_pred_rphylo[test_idx])

cat(sprintf("Rphylopars RMSE | train-miss: %.4f | val-miss: %.4f | test-miss: %.4f\n",
            rmse_train_rphylo, rmse_val_rphylo, rmse_test_rphylo))

t_bm <- torch_tensor(as.matrix(X_pred_rphylo), dtype=torch_float(), device=device)

# ----------------------------- tree features (cached) -------------------------
cache_file <- paste0("phylo_dae_cache_cat_", n, "_k", k_eigen, "_s", sigma_mult, ".rds")
cat("\n--- Building tree features (spectral + adjacency) ---\n")
t_feat0 <- Sys.time()

if (file.exists(cache_file)) {
  cache <- readRDS(cache_file)
  if (isTRUE(nrow(cache$coords) == n) && isTRUE(nrow(cache$adj) == n)) {
    coords <- cache$coords; adj <- cache$adj
  } else {
    file.remove(cache_file)
    coords <- get_spectral_features(tree_pruned, k=k_eigen, sigma_mult=sigma_mult)
    adj <- get_stable_adj(tree_pruned, sigma_mult=sigma_mult)
    saveRDS(list(coords=coords, adj=adj), cache_file)
  }
} else {
  coords <- get_spectral_features(tree_pruned, k=k_eigen, sigma_mult=sigma_mult)
  adj <- get_stable_adj(tree_pruned, sigma_mult=sigma_mult)
  saveRDS(list(coords=coords, adj=adj), cache_file)
}

cat(sprintf("Tree features ready in %.2fs\n", as.numeric(difftime(Sys.time(), t_feat0, units="secs"))))

t_coords <- torch_tensor(coords, dtype=torch_float(), device=device)
t_adj    <- torch_tensor(adj,    dtype=torch_float(), device=device)

# ----------------------------- categorical one-hot ----------------------------
oh <- make_onehot(df_cat, cat_cols_found, max_levels = 40)
C_mat <- oh$mat
cat("Using categorical cols:", paste(oh$kept, collapse=", "), "\n")
cat("One-hot dim:", ncol(C_mat), "\n")

t_cat <- torch_tensor(C_mat, dtype=torch_float(), device=device)
q_cat <- ncol(C_mat)

# ----------------------------- model: bounded residual-on-BM ------------------
DistilledPhyloDAE_cat_residBM <- nn_module(
  "DistilledPhyloDAE_cat_residBM",
  initialize = function(n_traits, k_eigen, q_cat, hidden=192, dropout=0.10, res_cap=0.30) {
    input_dim <- n_traits + n_traits + k_eigen + q_cat
    self$enc1 <- nn_linear(input_dim, hidden)
    self$enc2 <- nn_linear(hidden, hidden)
    self$mix  <- nn_linear(hidden, hidden)
    self$dec1 <- nn_linear(hidden, hidden)
    self$dec2 <- nn_linear(hidden, n_traits)
    
    self$act  <- nn_gelu()
    self$drop <- nn_dropout(dropout)
    
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits))
    
    # NEW: bounded positive scale
    self$res_raw <- nn_parameter(torch_tensor(-4, dtype=torch_float()))  # sigmoid(-4) ~ 0.018
    self$res_cap <- res_cap
  },
  
  forward = function(x, m, coords, cats, adj, bm) {
    z <- torch_cat(list(x, m, coords, cats), dim=2)
    
    h <- self$enc1(z)
    h <- self$act(h)
    h <- self$drop(h)
    
    h <- self$enc2(h)
    h <- self$act(h)
    
    neigh <- torch_matmul(adj, h)
    h <- h + self$mix(neigh)
    
    h <- self$dec1(h)
    h <- self$act(h)
    
    resid <- self$dec2(h)
    
    rs <- torch_sigmoid(self$res_raw) * self$res_cap
    bm + rs * resid
  },
  
  get_res_scale = function() {
    as.numeric((torch_sigmoid(self$res_raw) * self$res_cap)$detach()$cpu()$item())
  }
)

# ----------------------------- trainer ----------------------------------------
run_model <- function(
    hidden=CFG$hidden, dropout=CFG$dropout,
    lr=CFG$lr, weight_decay=CFG$weight_decay,
    corrupt_p=CFG$corrupt_p,
    lambda_anchor=CFG$lambda_anchor,
    lambda_trainmiss=CFG$lambda_trainmiss,
    lambda_resid=CFG$lambda_resid,
    lambda_distill_obs=CFG$lambda_distill_obs,
    sup_observe_drop=CFG$sup_observe_drop,
    res_cap=CFG$res_cap,
    lambda_scale=CFG$lambda_scale,
    max_epochs=CFG$max_epochs, eval_every=CFG$eval_every, patience=CFG$patience,
    clip_norm=CFG$clip_norm,
    seed=CFG$seed,
    verbose=TRUE
) {
  set.seed(seed); torch_manual_seed(seed)
  
  model <- DistilledPhyloDAE_cat_residBM(
    n_traits=p, k_eigen=ncol(coords), q_cat=q_cat,
    hidden=hidden, dropout=dropout, res_cap=res_cap
  )
  model$to(device=device)
  opt <- optim_adamw(model$parameters, lr=lr, weight_decay=weight_decay)
  
  best_val <- Inf
  best_state <- NULL
  patience_left <- patience
  t0 <- Sys.time()
  
  missing_bool <- (t_mask == 0)
  
  for (epoch in 1:max_epochs) {
    model$train()
    opt$zero_grad()
    
    # Corrupt only observed entries
    u <- torch_rand_like(t_X)
    corrupt <- (u < corrupt_p) & (t_mask == 1)
    if (as.numeric(corrupt$sum()$cpu()$item()) == 0) next
    corrupt_f <- corrupt$to(dtype=torch_float())
    
    mtok <- model$mask_token$to(device=device)$expand_as(t_X)
    
    # Base input: missing -> token; corrupted observed -> token
    X_in_t <- torch_where(missing_bool, mtok, t_X)
    X_in_t <- torch_where(corrupt, mtok, X_in_t)
    
    # Base dynamic mask: observed & not corrupted
    mask_dyn <- t_mask * (1 - corrupt_f)
    
    # -------- NEW: extra observed dropout used for supervised signal ----------
    # We hide additional observed cells to prevent train-missing leakage.
    if (sup_observe_drop > 0) {
      u2 <- torch_rand_like(t_X)
      extra_drop <- (u2 < sup_observe_drop) & (mask_dyn == 1)
      extra_drop_f <- extra_drop$to(dtype=torch_float())
      
      X_in_sup <- torch_where(extra_drop, mtok, X_in_t)
      mask_sup <- mask_dyn * (1 - extra_drop_f)
    } else {
      X_in_sup <- X_in_t
      mask_sup <- mask_dyn
    }
    
    # Forward pass (use the supervised-masked version)
    pred <- model(X_in_sup, mask_sup, t_coords, t_cat, t_adj, t_bm)
    
    # 1) denoise loss on corrupted observed cells (still valid)
    loss_rec <- nnf_mse_loss(pred[corrupt], t_X[corrupt])
    
    # 2) anchor on observed & not corrupted & not extra-dropped
    anchor_idx <- (mask_sup == 1)
    loss_anchor <- nnf_mse_loss(pred[anchor_idx], t_X[anchor_idx])
    
    # 3) supervised learning on TRAIN-missing cells
    loss_trainmiss <- nnf_mse_loss(pred[t_train_mask], t_truth[t_train_mask])
    
    # 4) shrink deviation from BM on missing cells
    loss_resid <- torch_mean((pred[missing_bool] - t_bm[missing_bool])$pow(2))
    
    # 5) optional distill on observed cells
    if (lambda_distill_obs > 0) {
      loss_distill_obs <- nnf_mse_loss(pred[t_mask == 1], t_bm[t_mask == 1])
    } else {
      loss_distill_obs <- torch_tensor(0, device=device)
    }
    
    # 6) optional penalty on residual scale itself
    if (lambda_scale > 0) {
      rs_t <- (torch_sigmoid(model$res_raw) * model$res_cap)
      loss_scale <- rs_t$pow(2)
    } else {
      loss_scale <- torch_tensor(0, device=device)
    }
    
    loss <- loss_rec +
      lambda_anchor      * loss_anchor +
      lambda_trainmiss   * loss_trainmiss +
      lambda_resid       * loss_resid +
      lambda_distill_obs * loss_distill_obs +
      lambda_scale       * loss_scale
    
    loss$backward()
    
    if (!all_grads_finite(model$parameters)) {
      opt$zero_grad()
      next
    }
    nn_utils_clip_grad_norm_(model$parameters, max_norm = clip_norm)
    opt$step()
    
    if (epoch %% eval_every == 0) {
      model$eval()
      with_no_grad({
        mtok_eval <- model$mask_token$to(device=device)$expand_as(t_X)
        X_eval_in <- torch_where(missing_bool, mtok_eval, t_X)
        
        pred0 <- model(X_eval_in, t_mask, t_coords, t_cat, t_adj, t_bm)
        Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
        
        train_rmse <- as.numeric(rmse_torch(Xcurr, t_truth, t_train_mask)$cpu()$item())
        val_rmse   <- as.numeric(rmse_torch(Xcurr, t_truth, t_val_mask)$cpu()$item())
        rs <- model$get_res_scale()
      })
      
      if (val_rmse + 1e-6 < best_val) {
        best_val <- val_rmse
        best_state <- model$state_dict()
        patience_left <- patience
      } else {
        patience_left <- patience_left - 1
      }
      
      if (verbose) {
        elapsed <- as.numeric(difftime(Sys.time(), t0, units="secs"))
        cat(sprintf(
          "epoch %4d | rec %.5f | anc %.5f | trainMiss %.5f | residReg %.5f | rs %.4f | trainRMSE %.4f | valRMSE %.4f | best %.4f | pat %d | t=%.1fs\n",
          epoch,
          loss_rec$item(), loss_anchor$item(), loss_trainmiss$item(), loss_resid$item(),
          rs, train_rmse, val_rmse, best_val, patience_left, elapsed
        ))
      }
      if (patience_left <= 0) break
    }
  }
  
  if (!is.null(best_state)) model$load_state_dict(best_state)
  
  model$eval()
  with_no_grad({
    mtok_eval <- model$mask_token$to(device=device)$expand_as(t_X)
    missing_bool <- (t_mask == 0)
    X_eval_in <- torch_where(missing_bool, mtok_eval, t_X)
    
    pred0 <- model(X_eval_in, t_mask, t_coords, t_cat, t_adj, t_bm)
    Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
    
    train_rmse <- as.numeric(rmse_torch(Xcurr, t_truth, t_train_mask)$cpu()$item())
    val_rmse   <- as.numeric(rmse_torch(Xcurr, t_truth, t_val_mask)$cpu()$item())
    test_rmse  <- as.numeric(rmse_torch(Xcurr, t_truth, t_test_mask)$cpu()$item())
    
    pred_mat <- as.matrix(Xcurr$cpu())
    res_scale_final <- model$get_res_scale()
  })
  
  list(pred=pred_mat,
       train_rmse=train_rmse, val_rmse=val_rmse, test_rmse=test_rmse,
       res_scale=res_scale_final)
}

# ----------------------------- run a small sweep ------------------------------
cat("\n--- Running improved Phylo-DAE + cats (bounded rs + supervised leakage control) ---\n")

cand <- list(
  list(lambda_trainmiss=0.5, lambda_resid=0.5, sup_observe_drop=0.10),
  list(lambda_trainmiss=0.5, lambda_resid=1.0, sup_observe_drop=0.10),
  list(lambda_trainmiss=0.2, lambda_resid=0.5, sup_observe_drop=0.15),
  list(lambda_trainmiss=0.5, lambda_resid=0.5, sup_observe_drop=0.20)
)

best <- NULL
best_val <- Inf

for (i in seq_along(cand)) {
  cfg <- cand[[i]]
  cat(sprintf("  [%d/%d] trainMissW=%.2f residW=%.2f supDrop=%.2f\n",
              i, length(cand), cfg$lambda_trainmiss, cfg$lambda_resid, cfg$sup_observe_drop))
  
  r <- do.call(run_model, c(cfg, list(
    hidden=CFG$hidden, dropout=CFG$dropout,
    lr=CFG$lr, weight_decay=CFG$weight_decay,
    corrupt_p=CFG$corrupt_p,
    lambda_anchor=CFG$lambda_anchor,
    lambda_distill_obs=CFG$lambda_distill_obs,
    res_cap=CFG$res_cap, lambda_scale=CFG$lambda_scale,
    max_epochs=CFG$max_epochs, eval_every=CFG$eval_every, patience=CFG$patience,
    clip_norm=CFG$clip_norm,
    seed=CFG$seed,
    verbose=TRUE
  )))
  
  cat(sprintf("      trainMiss RMSE %.4f | valMiss RMSE %.4f | testMiss RMSE %.4f | res_scale %.4f\n",
              r$train_rmse, r$val_rmse, r$test_rmse, r$res_scale))
  
  if (r$val_rmse < best_val) { best_val <- r$val_rmse; best <- r }
}

X_pred_dae_cat <- best$pred

cat("\n=========================================\n")
cat(" RESULTS (RMSE on Standardized Log Data)\n")
cat("=========================================\n")
cat(sprintf("Rphylopars (BM)                 | train-miss: %.4f | val-miss: %.4f | test-miss: %.4f\n",
            rmse_train_rphylo, rmse_val_rphylo, rmse_test_rphylo))
cat(sprintf("Improved Phylo-DAE + cats        | train-miss: %.4f | val-miss: %.4f | test-miss: %.4f | res_scale %.4f\n",
            best$train_rmse, best$val_rmse, best$test_rmse, best$res_scale))
cat("=========================================\n")

# ----------------------------- visualization ---------------------------------
df_plot <- rbind(
  data.frame(Truth = X_truth[test_idx], Prediction = X_pred_rphylo[test_idx], Method = "Rphylopars (BM)"),
  data.frame(Truth = X_truth[test_idx], Prediction = X_pred_dae_cat[test_idx], Method = "Improved Phylo-DAE + cats")
)

p_plot <- ggplot(df_plot, aes(x=Truth, y=Prediction)) +
  geom_point(alpha=0.25, size=1) +
  geom_abline(linetype="dashed") +
  facet_wrap(~Method, scales="fixed") +
  theme_minimal() +
  labs(
    title="Imputation on Held-out Cells (TEST missing cells only)",
    subtitle=sprintf("Test-miss RMSE: DAE=%.3f vs BM=%.3f (n=%d)", best$test_rmse, rmse_test_rphylo, n)
  )
print(p_plot)

diff_mat <- X_pred_dae_cat - X_pred_rphylo
cat(sprintf("\nmean(|DAE - BM|) on TEST-missing cells: %.6f\n", mean(abs(diff_mat[test_idx]))))
cat(sprintf("sd(DAE - BM) on TEST-missing cells: %.6f\n", sd(as.numeric(diff_mat[test_idx]))))