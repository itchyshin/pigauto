# ==============================================================================
# Distilled Phylo-DAE + categorical predictors (one-hot)
# - Keeps BM teacher + denoising training
# - Adds one-hot categorical covariates to DAE ONLY (Rphylopars baseline unchanged)
# - Early stopping on held-out (val) cells
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(dplyr)
library(here)

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

# robust one-hot builder
make_onehot <- function(df, cat_cols, max_levels = 40) {
  dd <- df[, cat_cols, drop = FALSE]
  # coerce to factor; NA -> "Unknown"; drop absurdly high-cardinality cols
  keep <- sapply(dd, function(x) {
    ux <- unique(x)
    ux <- ux[!is.na(ux)]
    length(ux) <= max_levels
  })
  if (!any(keep)) stop("No categorical columns left after max_levels filter.")
  
  dd <- dd[, keep, drop = FALSE]
  for (j in seq_along(dd)) {
    x <- dd[[j]]
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "Unknown"
    dd[[j]] <- factor(x)
  }
  
  # one-hot without intercept; no interactions
  mm <- model.matrix(~ . - 1, data = dd)
  mm <- as.matrix(mm)
  storage.mode(mm) <- "double"
  list(mat = mm, colnames = colnames(mm), kept = names(dd))
}

# ----------------------------- device -----------------------------------------
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

# ----------------------------- data -------------------------------------------
tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

# >>> EDIT THIS: choose categorical predictors available in your AVONET file <<<
# Examples you might have (depends on your AVONET version):
# cat_cols <- c("Primary.Lifestyle", "Diet", "Migration", "Habitat", "Nest.Type")
# If unsure, run: names(avonet)
cat_cols <- c("Diet", "Migration")  # <-- CHANGE ME

df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

# bring in cats aligned to df_final order
df_cat <- avonet[match(df_final$Species_Key, avonet$Species_Key), c("Species_Key", cat_cols), drop = FALSE]
stopifnot(all(df_cat$Species_Key == df_final$Species_Key))

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

# tensors
t_X     <- torch_tensor(X_fill,  dtype=torch_float(), device=device)
t_mask  <- torch_tensor(mask_input, dtype=torch_float(), device=device)
t_truth <- torch_tensor(X_truth, dtype=torch_float(), device=device)

val_mask_mat <- matrix(FALSE, nrow=n, ncol=p);  val_mask_mat[val_idx]  <- TRUE
test_mask_mat <- matrix(FALSE, nrow=n, ncol=p); test_mask_mat[test_idx] <- TRUE
t_val_mask  <- torch_tensor(val_mask_mat,  dtype=torch_bool(), device=device)
t_test_mask <- torch_tensor(test_mask_mat, dtype=torch_bool(), device=device)

# ----------------------------- baseline: Rphylopars (BM) ----------------------
cat("\n--- Running Rphylopars (BM) ---\n")
df_rphylo <- data.frame(species = df_final$Species_Key, X_in)
res_rphylo <- phylopars(trait_data=df_rphylo, tree=tree_pruned, model="BM", pheno_error=FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]

rmse_val_rphylo  <- rmse_vec(X_truth[val_idx],  X_pred_rphylo[val_idx])
rmse_test_rphylo <- rmse_vec(X_truth[test_idx], X_pred_rphylo[test_idx])
cat(sprintf("Rphylopars RMSE | val: %.4f | test: %.4f\n", rmse_val_rphylo, rmse_test_rphylo))

# BM teacher tensor (allowed: it uses X_in where val/test are missing)
t_bm <- torch_tensor(as.matrix(X_pred_rphylo), dtype=torch_float(), device=device)

# ----------------------------- tree features (spectral + adjacency) -----------
get_spectral_features <- function(tree, k = 16) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.35
  A <- exp(-(D^2) / (2 * sigma^2))
  diag(A) <- 0
  L <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric = TRUE)
  n <- nrow(A)
  coords <- eig$vectors[, (n - k + 1):n, drop=FALSE]
  coords
}

get_stable_adj <- function(tree) {
  D <- cophenetic(tree)
  sigma <- median(D) * 0.35
  A <- exp(-(D^2) / (2 * sigma^2))
  diag(A) <- 0
  rs <- rowSums(A); rs[rs == 0] <- 1
  A_norm <- diag(1/rs) %*% A
  A_norm
}

# cache to avoid repeated eigendecomp
cache_file <- here("avonet", "phylo_dae_cache_cat.rds")
if (file.exists(cache_file)) {
  cache <- readRDS(cache_file)
  if (nrow(cache$coords) == n && nrow(cache$adj) == n) {
    coords <- cache$coords
    adj <- cache$adj
  } else {
    file.remove(cache_file)
    coords <- get_spectral_features(tree_pruned, k=16)
    adj <- get_stable_adj(tree_pruned)
    saveRDS(list(coords=coords, adj=adj), cache_file)
  }
} else {
  coords <- get_spectral_features(tree_pruned, k=16)
  adj <- get_stable_adj(tree_pruned)
  saveRDS(list(coords=coords, adj=adj), cache_file)
}

k_eigen <- ncol(coords)
t_coords <- torch_tensor(coords, dtype=torch_float(), device=device)
t_adj    <- torch_tensor(adj, dtype=torch_float(), device=device)

# ----------------------------- categorical one-hot ----------------------------
oh <- make_onehot(df_cat, cat_cols, max_levels = 40)
C_mat <- oh$mat
cat("Using categorical cols:", paste(oh$kept, collapse=", "), "\n")
cat("One-hot dim:", ncol(C_mat), "\n")

t_cat <- torch_tensor(C_mat, dtype=torch_float(), device=device)
q_cat <- ncol(C_mat)

# ----------------------------- model: Distilled Phylo-DAE + cats --------------
DistilledPhyloDAE_cat <- nn_module(
  "DistilledPhyloDAE_cat",
  initialize = function(n_traits, k_eigen, q_cat, hidden=192, dropout=0.10) {
    input_dim <- n_traits + n_traits + k_eigen + q_cat  # x + mask + coords + cats
    self$enc1 <- nn_linear(input_dim, hidden)
    self$enc2 <- nn_linear(hidden, hidden)
    
    self$mix  <- nn_linear(hidden, hidden)
    
    self$dec1 <- nn_linear(hidden, hidden)
    self$dec2 <- nn_linear(hidden, n_traits)
    
    self$act  <- nn_gelu()
    self$drop <- nn_dropout(dropout)
    
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits))
  },
  
  forward = function(x, m, coords, cats, adj) {
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
    out <- self$dec2(h)
    out
  }
)

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

# ----------------------------- runner -----------------------------------------
run_distilled_phylo_dae_cat <- function(
    hidden=192, dropout=0.10,
    lr=2e-3, weight_decay=1e-5,
    corrupt_p=0.30,
    lambda_distill=0.20,      # lower than your BM-match run, to allow deviation
    lambda_anchor=0.05,
    max_epochs=2500, eval_every=100, patience=12,
    refine_steps=6,
    clip_norm=1.0,
    seed=1, verbose=TRUE
) {
  set.seed(seed); torch_manual_seed(seed)
  
  model <- DistilledPhyloDAE_cat(n_traits=p, k_eigen=k_eigen, q_cat=q_cat, hidden=hidden, dropout=dropout)
  model$to(device=device)
  opt <- optim_adamw(model$parameters, lr=lr, weight_decay=weight_decay)
  
  best_val <- Inf
  best_state <- NULL
  patience_left <- patience
  
  for (epoch in 1:max_epochs) {
    model$train()
    opt$zero_grad()
    
    # corrupt only observed entries
    u <- torch_rand_like(t_X)
    corrupt <- (u < corrupt_p) & (t_mask == 1)
    if (as.numeric(corrupt$sum()$cpu()$item()) == 0) next
    
    corrupt_f <- corrupt$to(dtype=torch_float())
    X_in_t <- torch_where(corrupt, model$mask_token$to(device=device)$expand_as(t_X), t_X)
    mask_dyn <- t_mask * (1 - corrupt_f)
    
    pred <- model(X_in_t, mask_dyn, t_coords, t_cat, t_adj)
    
    # 1) denoise loss on corrupted
    loss_rec <- nnf_mse_loss(pred[corrupt], t_X[corrupt])
    
    # 2) anchor on observed & not corrupted
    anchor_idx <- (mask_dyn == 1)
    loss_anchor <- nnf_mse_loss(pred[anchor_idx], t_X[anchor_idx])
    
    # 3) distill toward BM on (true-missing OR corrupted) cells
    #    (lets the net deviate on missing using cats, but keeps it stable)
    distill_idx <- (t_mask == 0) | corrupt
    loss_distill <- nnf_mse_loss(pred[distill_idx], t_bm[distill_idx])
    
    loss <- loss_rec + lambda_anchor * loss_anchor + lambda_distill * loss_distill
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
        pred0 <- model(t_X, t_mask, t_coords, t_cat, t_adj)
        Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
        
        if (refine_steps > 0) {
          for (k in 1:refine_steps) {
            pk <- model(Xcurr, t_mask, t_coords, t_cat, t_adj)
            Xcurr <- t_X * t_mask + pk * (1 - t_mask)
          }
        }
        
        val_rmse <- as.numeric(rmse_torch(Xcurr, t_truth, t_val_mask)$cpu()$item())
      })
      
      if (val_rmse + 1e-6 < best_val) {
        best_val <- val_rmse
        best_state <- model$state_dict()
        patience_left <- patience
      } else {
        patience_left <- patience_left - 1
      }
      
      if (verbose) {
        cat(sprintf("epoch %4d | rec %.4f | anc %.4f | dist %.4f | valRMSE %.4f | best %.4f | pat %d\n",
                    epoch, loss_rec$item(), loss_anchor$item(), loss_distill$item(),
                    val_rmse, best_val, patience_left))
      }
      if (patience_left <= 0) break
    }
  }
  
  if (!is.null(best_state)) model$load_state_dict(best_state)
  
  model$eval()
  with_no_grad({
    pred0 <- model(t_X, t_mask, t_coords, t_cat, t_adj)
    Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
    if (refine_steps > 0) {
      for (k in 1:refine_steps) {
        pk <- model(Xcurr, t_mask, t_coords, t_cat, t_adj)
        Xcurr <- t_X * t_mask + pk * (1 - t_mask)
      }
    }
    val_rmse <- as.numeric(rmse_torch(Xcurr, t_truth, t_val_mask)$cpu()$item())
    test_rmse <- as.numeric(rmse_torch(Xcurr, t_truth, t_test_mask)$cpu()$item())
    pred_mat <- as.matrix(Xcurr$cpu())
  })
  
  list(pred=pred_mat, val_rmse=val_rmse, test_rmse=test_rmse)
}

# ----------------------------- run a few configs ------------------------------
cat("\n--- Running Distilled Phylo-DAE + categorical predictors ---\n")

cand <- list(
  list(hidden=192, dropout=0.10, lr=2e-3, corrupt_p=0.30, lambda_distill=0.20, lambda_anchor=0.05),
  list(hidden=256, dropout=0.10, lr=1e-3, corrupt_p=0.35, lambda_distill=0.10, lambda_anchor=0.05),
  list(hidden=256, dropout=0.15, lr=1e-3, corrupt_p=0.35, lambda_distill=0.05, lambda_anchor=0.05)
)

best <- NULL
best_val <- Inf

for (i in seq_along(cand)) {
  cfg <- cand[[i]]
  cat(sprintf("  [%d/%d] hidden=%d drop=%.2f lr=%.0e corrupt=%.2f dist=%.2f anc=%.2f\n",
              i, length(cand), cfg$hidden, cfg$dropout, cfg$lr, cfg$corrupt_p,
              cfg$lambda_distill, cfg$lambda_anchor))
  r <- do.call(run_distilled_phylo_dae_cat, c(cfg, list(
    max_epochs=2500, eval_every=100, patience=12, refine_steps=6, seed=1, verbose=FALSE
  )))
  cat(sprintf("      val RMSE %.4f | test RMSE %.4f\n", r$val_rmse, r$test_rmse))
  if (r$val_rmse < best_val) { best_val <- r$val_rmse; best <- r }
}

X_pred_dae_cat <- best$pred

cat("\n=========================================\n")
cat(" RESULTS (RMSE on Standardized Log Data)\n")
cat("=========================================\n")
cat(sprintf("Rphylopars (BM)                | val: %.4f | test: %.4f\n", rmse_val_rphylo, rmse_test_rphylo))
cat(sprintf("Distilled Phylo-DAE + cats      | val: %.4f | test: %.4f\n", best$val_rmse, best$test_rmse))
cat("=========================================\n")

# ----------------------------- visualization ---------------------------------
df_plot <- rbind(
  data.frame(Truth = X_truth[test_idx], Prediction = X_pred_rphylo[test_idx], Method = "Rphylopars (BM)"),
  data.frame(Truth = X_truth[test_idx], Prediction = X_pred_dae_cat[test_idx], Method = "Distilled Phylo-DAE + cats")
)

p_plot <- ggplot(df_plot, aes(x=Truth, y=Prediction)) +
  geom_point(alpha=0.25, size=1) +
  geom_abline(linetype="dashed") +
  facet_wrap(~Method, scales="free") +
  theme_minimal() +
  labs(
    title="Imputation on Held-out Cells (test split only)",
    subtitle=sprintf("Test RMSE: DAE+categ=%.3f vs Rphylopars(BM)=%.3f",
                     best$test_rmse, rmse_test_rphylo)
  )
print(p_plot)