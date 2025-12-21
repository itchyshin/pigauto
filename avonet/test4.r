# ==============================================================================
# Distilled Phylo-DAE + categorical predictors (one-hot)  [BM-RESIDUAL + FAST MODE]
# - BM teacher + denoising
# - Adds one-hot categorical covariates to DAE ONLY
# - SAFE residual-on-BM: pred = BM + res_scale * residual_net(...)
#     * starts near BM, only deviates if it helps
# - Early stopping on held-out (val) cells
# - Robust to missing categorical columns + robust file paths
# - FAST mode via n_sub (subset tips to speed up dense adjacency)
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(dplyr)
library(here)

# ----------------------------- USER SETTINGS ----------------------------------

# If too slow: n_sub <- 500 or 300.  NULL uses all tips (slow with dense adj).
n_sub <- 1200

# Requested categorical predictors (script drops those that don't exist)
cat_cols_requested <- c("Diet", "Migration")  # <-- edit

# Paths: script tries ./avonet/... then /mnt/data/...
tree_path_A <- here("avonet","Stage2_Hackett_MCC_no_neg.tre")
csv_path_A  <- here("avonet","AVONET3_BirdTree.csv")

tree_path_B <- "/mnt/data/Stage2_Hackett_MCC_no_neg.tre"
csv_path_B  <- "/mnt/data/AVONET3_BirdTree.csv"

# ----------------------------- helpers ----------------------------------------
rmse_vec <- function(truth, pred) sqrt(mean((truth - pred)^2))

make_missing_splits <- function(X, missing_frac = 0.25, val_frac = 0.25, seed = 555) {
  set.seed(seed)
  n <- nrow(X); p <- ncol(X)
  all_idx <- seq_len(n * p)
  m <- floor(missing_frac * length(all_idx))
  miss <- sample(all_idx, m)
  m_val <- floor(val_frac * length(miss))
  val_idx <- miss[seq_len(m_val)]
  test_idx <- miss[(m_val + 1):length(miss)]
  list(val_idx = val_idx, test_idx = test_idx)
}

rmse_torch <- function(pred, truth, mask_bool) {
  d <- (pred - truth)
  torch_sqrt(torch_mean(d[mask_bool]$pow(2)))
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

# robust one-hot builder (drops high-cardinality + handles NA)
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
  
  mm <- model.matrix(~ . - 1, data = dd) # one-hot, no intercept
  mm <- as.matrix(mm)
  storage.mode(mm) <- "double"
  list(mat = mm, colnames = colnames(mm), kept = names(dd))
}

# tree features: spectral coords + dense normalized adjacency (cached)
get_spectral_features <- function(tree, k = 16, sigma_mult = 0.35) {
  D <- cophenetic(tree)
  sigma <- median(D) * sigma_mult
  A <- exp(-(D^2) / (2 * sigma^2))
  diag(A) <- 0
  L <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric = TRUE)
  n <- nrow(A)
  eig$vectors[, (n - k + 1):n, drop = FALSE]
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

trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]
stopifnot(all(df_final$Species_Key == tree_pruned$tip.label))

# align categorical data
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

# FAST MODE: subset tips (speeds up dense adj massively)
if (!is.null(n_sub) && n_sub < nrow(df_final)) {
  set.seed(1)
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

# ----------------------------- missingness ------------------------------------
splits <- make_missing_splits(X_truth, missing_frac = 0.25, val_frac = 0.25, seed = 555)
val_idx  <- splits$val_idx
test_idx <- splits$test_idx

mask_input <- matrix(1, nrow=n, ncol=p)
mask_input[c(val_idx, test_idx)] <- 0

X_in <- X_truth
X_in[c(val_idx, test_idx)] <- NA
X_fill <- X_in
X_fill[is.na(X_fill)] <- 0

t_X     <- torch_tensor(X_fill, dtype=torch_float(), device=device)
t_mask  <- torch_tensor(mask_input, dtype=torch_float(), device=device)
t_truth <- torch_tensor(X_truth, dtype=torch_float(), device=device)

val_mask_mat <- matrix(FALSE, nrow=n, ncol=p);  val_mask_mat[val_idx]  <- TRUE
test_mask_mat <- matrix(FALSE, nrow=n, ncol=p); test_mask_mat[test_idx] <- TRUE
t_val_mask  <- torch_tensor(val_mask_mat, dtype=torch_bool(), device=device)
t_test_mask <- torch_tensor(test_mask_mat, dtype=torch_bool(), device=device)

# ----------------------------- baseline: Rphylopars (BM) ----------------------
cat("\n--- Running Rphylopars (BM) ---\n")
df_rphylo <- data.frame(species = df_final$Species_Key, X_in)
res_rphylo <- phylopars(trait_data=df_rphylo, tree=tree_pruned, model="BM", pheno_error=FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]

rmse_val_rphylo  <- rmse_vec(X_truth[val_idx],  X_pred_rphylo[val_idx])
rmse_test_rphylo <- rmse_vec(X_truth[test_idx], X_pred_rphylo[test_idx])
cat(sprintf("Rphylopars RMSE | val: %.4f | test: %.4f\n", rmse_val_rphylo, rmse_test_rphylo))

t_bm <- torch_tensor(as.matrix(X_pred_rphylo), dtype=torch_float(), device=device)

# ----------------------------- tree features (cached) -------------------------
cache_file <- paste0("phylo_dae_cache_cat_", n, ".rds")
cat("\n--- Building tree features (spectral + adjacency) ---\n")
t_feat0 <- Sys.time()

if (file.exists(cache_file)) {
  cache <- readRDS(cache_file)
  if (isTRUE(nrow(cache$coords) == n) && isTRUE(nrow(cache$adj) == n)) {
    coords <- cache$coords; adj <- cache$adj
  } else {
    file.remove(cache_file)
    coords <- get_spectral_features(tree_pruned, k=16, sigma_mult=0.35)
    adj <- get_stable_adj(tree_pruned, sigma_mult=0.35)
    saveRDS(list(coords=coords, adj=adj), cache_file)
  }
} else {
  coords <- get_spectral_features(tree_pruned, k=16, sigma_mult=0.35)
  adj <- get_stable_adj(tree_pruned, sigma_mult=0.35)
  saveRDS(list(coords=coords, adj=adj), cache_file)
}

cat(sprintf("Tree features ready in %.2fs\n", as.numeric(difftime(Sys.time(), t_feat0, units="secs"))))

k_eigen <- ncol(coords)
t_coords <- torch_tensor(coords, dtype=torch_float(), device=device)
t_adj    <- torch_tensor(adj, dtype=torch_float(), device=device)

# ----------------------------- categorical one-hot ----------------------------
oh <- make_onehot(df_cat, cat_cols_found, max_levels = 40)
C_mat <- oh$mat
cat("Using categorical cols:", paste(oh$kept, collapse=", "), "\n")
cat("One-hot dim:", ncol(C_mat), "\n")

t_cat <- torch_tensor(C_mat, dtype=torch_float(), device=device)
q_cat <- ncol(C_mat)

# ----------------------------- model: SAFE residual-on-BM ---------------------
# Output: pred = BM + res_scale * residual
# res_scale starts small (0.01) so the model begins very close to BM.

DistilledPhyloDAE_cat_residBM <- nn_module(
  "DistilledPhyloDAE_cat_residBM",
  initialize = function(n_traits, k_eigen, q_cat, hidden=192, dropout=0.10) {
    input_dim <- n_traits + n_traits + k_eigen + q_cat
    self$enc1 <- nn_linear(input_dim, hidden)
    self$enc2 <- nn_linear(hidden, hidden)
    self$mix  <- nn_linear(hidden, hidden)
    self$dec1 <- nn_linear(hidden, hidden)
    self$dec2 <- nn_linear(hidden, n_traits)
    
    self$act  <- nn_gelu()
    self$drop <- nn_dropout(dropout)
    
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits))
    self$res_scale  <- nn_parameter(torch_tensor(0.01, dtype=torch_float()))  # << key
  },
  
  forward = function(x, m, coords, cats, adj, bm) {
    # x,m,bm: (n x p); coords: (n x k); cats: (n x q); adj: (n x n)
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
    
    # SAFE residual-on-BM
    out <- bm + self$res_scale * resid
    out
  }
)

# ----------------------------- runner -----------------------------------------
run_distilled_phylo_dae_cat <- function(
    hidden=192, dropout=0.10,
    lr=2e-3, weight_decay=1e-5,
    corrupt_p=0.30,
    lambda_distill=0.10,     # lower is ok now (residual structure is safe)
    lambda_anchor=0.05,
    lambda_resid = 0.20,     # << shrinkage: keep pred close to BM unless helpful
    max_epochs=1200, eval_every=100, patience=8,
    refine_steps=0,
    clip_norm=1.0,
    seed=1, verbose=TRUE
) {
  set.seed(seed); torch_manual_seed(seed)
  
  model <- DistilledPhyloDAE_cat_residBM(
    n_traits=p, k_eigen=k_eigen, q_cat=q_cat,
    hidden=hidden, dropout=dropout
  )
  model$to(device=device)
  opt <- optim_adamw(model$parameters, lr=lr, weight_decay=weight_decay)
  
  best_val <- Inf
  best_state <- NULL
  patience_left <- patience
  t0 <- Sys.time()
  
  for (epoch in 1:max_epochs) {
    model$train()
    opt$zero_grad()
    
    # corrupt only observed entries
    u <- torch_rand_like(t_X)
    corrupt <- (u < corrupt_p) & (t_mask == 1)
    if (as.numeric(corrupt$sum()$cpu()$item()) == 0) next
    
    corrupt_f <- corrupt$to(dtype=torch_float())
    
    # mask token
    mtok <- model$mask_token$to(device=device)$expand_as(t_X)
    X_in_t <- torch_where(corrupt, mtok, t_X)
    
    mask_dyn <- t_mask * (1 - corrupt_f)
    
    pred <- model(X_in_t, mask_dyn, t_coords, t_cat, t_adj, t_bm)
    
    # 1) denoise loss on corrupted
    loss_rec <- nnf_mse_loss(pred[corrupt], t_X[corrupt])
    
    # 2) anchor on observed & not corrupted
    anchor_idx <- (mask_dyn == 1)
    loss_anchor <- nnf_mse_loss(pred[anchor_idx], t_X[anchor_idx])
    
    # 3) distill toward BM on (true-missing OR corrupted)
    distill_idx <- (t_mask == 0) | corrupt
    loss_distill <- nnf_mse_loss(pred[distill_idx], t_bm[distill_idx])
    
    # 4) shrink residual globally (keeps you close to BM unless necessary)
    loss_resid <- torch_mean((pred - t_bm)$pow(2))
    
    loss <- loss_rec +
      lambda_anchor * loss_anchor +
      lambda_distill * loss_distill +
      lambda_resid * loss_resid
    
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
        pred0 <- model(t_X, t_mask, t_coords, t_cat, t_adj, t_bm)
        Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
        
        if (refine_steps > 0) {
          for (k in 1:refine_steps) {
            pk <- model(Xcurr, t_mask, t_coords, t_cat, t_adj, t_bm)
            Xcurr <- t_X * t_mask + pk * (1 - t_mask)
          }
        }
        
        val_rmse <- as.numeric(rmse_torch(Xcurr, t_truth, t_val_mask)$cpu()$item())
        rs <- as.numeric(model$res_scale$detach()$cpu()$item())
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
          "epoch %4d | rec %.5f | anc %.5f | dist %.5f | resid %.5f | res_scale %.4f | valRMSE %.4f | best %.4f | pat %d | t=%.1fs\n",
          epoch,
          loss_rec$item(), loss_anchor$item(), loss_distill$item(), loss_resid$item(),
          rs, val_rmse, best_val, patience_left, elapsed
        ))
      }
      if (patience_left <= 0) break
    }
  }
  
  if (!is.null(best_state)) model$load_state_dict(best_state)
  
  model$eval()
  with_no_grad({
    pred0 <- model(t_X, t_mask, t_coords, t_cat, t_adj, t_bm)
    Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
    
    val_rmse  <- as.numeric(rmse_torch(Xcurr, t_truth, t_val_mask)$cpu()$item())
    test_rmse <- as.numeric(rmse_torch(Xcurr, t_truth, t_test_mask)$cpu()$item())
    pred_mat <- as.matrix(Xcurr$cpu())
    res_scale_final <- as.numeric(model$res_scale$detach()$cpu()$item())
  })
  
  list(pred=pred_mat, val_rmse=val_rmse, test_rmse=test_rmse, res_scale=res_scale_final)
}

# ----------------------------- run a few configs ------------------------------
cat("\n--- Running Distilled Phylo-DAE + categorical predictors (ONE-HOT, residual-on-BM) ---\n")

cand <- list(
  # safer defaults; should never explode away from BM
  list(hidden=192, dropout=0.10, lr=2e-3, corrupt_p=0.30, lambda_distill=0.10, lambda_anchor=0.05, lambda_resid=0.20),
  list(hidden=256, dropout=0.10, lr=1e-3, corrupt_p=0.35, lambda_distill=0.05, lambda_anchor=0.05, lambda_resid=0.30)
)

best <- NULL
best_val <- Inf

for (i in seq_along(cand)) {
  cfg <- cand[[i]]
  cat(sprintf(
    "  [%d/%d] hidden=%d drop=%.2f lr=%.0e corrupt=%.2f dist=%.2f anc=%.2f resid=%.2f\n",
    i, length(cand), cfg$hidden, cfg$dropout, cfg$lr, cfg$corrupt_p,
    cfg$lambda_distill, cfg$lambda_anchor, cfg$lambda_resid
  ))
  
  r <- do.call(run_distilled_phylo_dae_cat, c(cfg, list(
    max_epochs=1200, eval_every=100, patience=8, refine_steps=0, seed=1, verbose=TRUE
  )))
  cat(sprintf("      val RMSE %.4f | test RMSE %.4f | res_scale %.4f\n", r$val_rmse, r$test_rmse, r$res_scale))
  
  if (r$val_rmse < best_val) { best_val <- r$val_rmse; best <- r }
}

X_pred_dae_cat <- best$pred

cat("\n=========================================\n")
cat(" RESULTS (RMSE on Standardized Log Data)\n")
cat("=========================================\n")
cat(sprintf("Rphylopars (BM)                 | val: %.4f | test: %.4f\n", rmse_val_rphylo, rmse_test_rphylo))
cat(sprintf("Phylo-DAE + cats (OH, resid-BM) | val: %.4f | test: %.4f\n", best$val_rmse, best$test_rmse))
cat("=========================================\n")

# ----------------------------- visualization ---------------------------------
df_plot <- rbind(
  data.frame(Truth = X_truth[test_idx], Prediction = X_pred_rphylo[test_idx], Method = "Rphylopars (BM)"),
  data.frame(Truth = X_truth[test_idx], Prediction = X_pred_dae_cat[test_idx], Method = "Phylo-DAE + cats (OH, resid-BM)")
)

p_plot <- ggplot(df_plot, aes(x=Truth, y=Prediction)) +
  geom_point(alpha=0.25, size=1) +
  geom_abline(linetype="dashed") +
  facet_wrap(~Method, scales="fixed") +   # fixed scales = honest comparison
  theme_minimal() +
  labs(
    title="Imputation on Held-out Cells (test split only)",
    subtitle=sprintf("Test RMSE: DAE+cats=%.3f vs BM=%.3f (n=%d)", best$test_rmse, rmse_test_rphylo, n)
  )
print(p_plot)

# Optional quick diagnostic: how far did we move from BM on test cells?
diff_mat <- X_pred_dae_cat - X_pred_rphylo
cat(sprintf("\nmean(|DAE - BM|) on TEST cells: %.6f\n", mean(abs(diff_mat[test_idx]))))
cat(sprintf("sd(DAE - BM) on TEST cells: %.6f\n", sd(as.numeric(diff_mat[test_idx]))))