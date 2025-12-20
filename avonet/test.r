# ==============================================================================
# AVONET BENCHMARK (CLEAN): Rphylopars vs Tuned Phylo-DAE
# - Fix spectral features (smallest non-zero Laplacian eigenvectors)
# - Proper val/test split on artificially missing entries
# - Early stopping + small hyperparam search
# - LayerNorm-based residual GNN + optional graph smoothness penalty
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

make_tree_features <- function(tree, k = 16, sigma_scale = 0.35, cache_file = NULL) {
  if (!is.null(cache_file) && file.exists(cache_file)) {
    cache <- readRDS(cache_file)
    if (all(cache$tip.label == tree$tip.label) && cache$k == k && cache$sigma_scale == sigma_scale) {
      return(cache[c("adj", "phylo")])
    }
  }
  
  D <- cophenetic(tree)
  sigma <- median(D) * sigma_scale
  
  # Gaussian kernel similarity
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 1
  
  # symmetric normalization (GCN-style)
  rs <- rowSums(A) + 1e-8
  Dinv <- diag(1 / sqrt(rs))
  adj <- Dinv %*% A %*% Dinv
  
  # Laplacian eigenmaps on similarity graph (use zero diag for Laplacian)
  A_spec <- A
  diag(A_spec) <- 0
  L <- diag(rowSums(A_spec)) - A_spec
  
  eig <- eigen(L, symmetric = TRUE)
  ord <- order(eig$values)              # smallest first
  # skip the first ~0 eigenvalue; take next k eigenvectors
  # (robustly: pick those after the smallest)
  idx <- ord[2:(k + 1)]
  phylo <- eig$vectors[, idx, drop = FALSE]
  
  out <- list(adj = adj, phylo = phylo)
  if (!is.null(cache_file)) {
    saveRDS(list(adj = adj, phylo = phylo, tip.label = tree$tip.label, k = k, sigma_scale = sigma_scale),
            cache_file)
  }
  out
}

make_missing_splits <- function(X, missing_frac = 0.25, val_frac = 0.25, seed = 555) {
  set.seed(seed)
  n <- nrow(X); p <- ncol(X)
  all_idx <- seq_len(n * p)
  m <- floor(missing_frac * length(all_idx))
  miss <- sample(all_idx, m)
  
  # split missing into val/test (both are "true missing" for evaluation)
  m_val <- floor(val_frac * length(miss))
  val_idx <- miss[seq_len(m_val)]
  test_idx <- miss[(m_val + 1):length(miss)]
  
  list(val_idx = val_idx, test_idx = test_idx)
}

idx_to_rc <- function(idx, nrow) {
  r <- ((idx - 1) %% nrow) + 1
  c <- ((idx - 1) %/% nrow) + 1
  cbind(r, c)
}

# ----------------------------- device ----------------------------------------

device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

# ----------------------------- data ------------------------------------------

tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))

avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

cat("Dataset:", nrow(df_final), "species x", length(trait_cols), "traits\n")

# log + standardize
X_raw  <- as.matrix(df_final[, trait_cols])
X_log  <- log(X_raw)
X_truth <- scale(X_log)  # ground truth on standardized log scale

# ----------------------------- tree features ---------------------------------

cache_file <- here("avonet", "phylo_cache_morpho_v2.rds")
feat <- make_tree_features(tree_pruned, k = 16, sigma_scale = 0.35, cache_file = cache_file)

t_adj   <- torch_tensor(feat$adj,   dtype = torch_float(), device = device)
t_phylo <- torch_tensor(feat$phylo, dtype = torch_float(), device = device)

# ----------------------------- missingness -----------------------------------

splits <- make_missing_splits(X_truth, missing_frac = 0.25, val_frac = 0.25, seed = 555)
val_idx  <- splits$val_idx
test_idx <- splits$test_idx

# input mask: 1 observed, 0 missing (val+test are hidden from both methods)
mask_input <- matrix(1, nrow = nrow(X_truth), ncol = ncol(X_truth))
mask_input[c(val_idx, test_idx)] <- 0

X_in <- X_truth
X_in[c(val_idx, test_idx)] <- NA

# ----------------------------- baseline: Rphylopars --------------------------

cat("\n--- Running Rphylopars (BM) ---\n")

df_rphylo <- data.frame(species = df_final$Species_Key, X_in)
start_t <- Sys.time()
res_rphylo <- phylopars(trait_data = df_rphylo, tree = tree_pruned, model = "BM", pheno_error = FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]
end_t <- Sys.time()

rmse_val_rphylo  <- rmse_vec(X_truth[val_idx],  X_pred_rphylo[val_idx])
rmse_test_rphylo <- rmse_vec(X_truth[test_idx], X_pred_rphylo[test_idx])

cat(sprintf("Rphylopars RMSE | val: %.4f | test: %.4f | time: %.2fs\n",
            rmse_val_rphylo, rmse_test_rphylo, as.numeric(difftime(end_t, start_t, units="secs"))))

# ----------------------------- Phylo-DAE model --------------------------------
# Notes:
# - LayerNorm > BatchNorm for this setting
# - Two residual "graph blocks"
# - Optional smoothness penalty encourages phylogenetic smoothness in hidden states

PhyloDAE <- nn_module(
  "PhyloDAE",
  initialize = function(n_traits, n_phylo, hidden_dim = 128, dropout = 0.1) {
    input_dim <- n_traits + n_traits + n_phylo  # x + mask + phylo coords
    
    self$in_fc   <- nn_linear(input_dim, hidden_dim)
    self$in_ln   <- nn_layer_norm(hidden_dim)
    self$act     <- nn_gelu()
    self$drop    <- nn_dropout(dropout)
    
    # graph block 1
    self$g1_msg  <- nn_linear(hidden_dim, hidden_dim)
    self$g1_gate <- nn_linear(hidden_dim, hidden_dim)
    self$g1_ln   <- nn_layer_norm(hidden_dim)
    
    # graph block 2
    self$g2_msg  <- nn_linear(hidden_dim, hidden_dim)
    self$g2_gate <- nn_linear(hidden_dim, hidden_dim)
    self$g2_ln   <- nn_layer_norm(hidden_dim)
    
    # decoder
    self$dec_fc1 <- nn_linear(hidden_dim, hidden_dim)
    self$dec_ln1 <- nn_layer_norm(hidden_dim)
    self$dec_out <- nn_linear(hidden_dim, n_traits)
    
    # trainable mask token (what the net sees for corrupted observed entries)
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits))
  },
  
  forward = function(x, m, p, adj) {
    z <- torch_cat(list(x, m, p), dim = 2)
    
    h <- self$in_fc(z)
    h <- self$in_ln(h)
    h <- self$act(h)
    h <- self$drop(h)
    
    # graph block 1 (residual gated message passing)
    msg1 <- torch_matmul(adj, h)
    gate1 <- torch_sigmoid(self$g1_gate(h))
    h <- h + gate1 * self$g1_msg(msg1)
    h <- self$g1_ln(h)
    h <- self$act(h)
    h <- self$drop(h)
    
    # graph block 2
    msg2 <- torch_matmul(adj, h)
    gate2 <- torch_sigmoid(self$g2_gate(h))
    h <- h + gate2 * self$g2_msg(msg2)
    h <- self$g2_ln(h)
    h <- self$act(h)
    h <- self$drop(h)
    
    # decode
    h2 <- self$dec_fc1(h)
    h2 <- self$dec_ln1(h2)
    h2 <- self$act(h2)
    out <- self$dec_out(h2)
    out
  }
)

# ----------------------------- training + tuning -----------------------------

run_train_eval <- function(cfg, seed = 1) {
  set.seed(seed)
  torch_manual_seed(seed)
  
  n <- nrow(X_truth); p <- ncol(X_truth)
  
  # tensors for input (val+test hidden)
  X_fill <- X_in
  X_fill[is.na(X_fill)] <- 0
  
  t_X    <- torch_tensor(X_fill,      dtype = torch_float(), device = device)
  t_mask <- torch_tensor(mask_input,  dtype = torch_float(), device = device)
  
  model <- PhyloDAE(n_traits = p, n_phylo = ncol(feat$phylo),
                    hidden_dim = cfg$hidden_dim, dropout = cfg$dropout)
  model$to(device = device)
  
  opt <- optim_adamw(model$parameters, lr = cfg$lr, weight_decay = cfg$weight_decay)
  
  best_val <- Inf
  best_state <- NULL
  patience_left <- cfg$patience
  
  # precompute val/test indices in tensor-friendly form
  val_rc  <- idx_to_rc(val_idx,  n)
  test_rc <- idx_to_rc(test_idx, n)
  
  # helper: compute RMSE on a set of linear indices
  rmse_on_idx <- function(pred_mat, idx_lin) {
    sqrt(mean((X_truth[idx_lin] - pred_mat[idx_lin])^2))
  }
  
  for (epoch in 1:cfg$max_epochs) {
    model$train()
    opt$zero_grad()
    
    # self-supervised corruption: only corrupt truly observed (mask_input == 1)
    u <- torch_rand_like(t_X)
    corrupt <- (u < cfg$corrupt_p) & (t_mask == 1)
    
    if (as.numeric(corrupt$sum()$cpu()) == 0) next
    
    corrupt_f <- corrupt$to(dtype = torch_float())
    
    # replace corrupted observed entries with mask_token
    X_in_t <- torch_where(corrupt, model$mask_token$to(device=device)$expand_as(t_X), t_X)
    mask_dyn <- t_mask * (1 - corrupt_f)  # tell model those are "missing" during this pass
    
    pred <- model(X_in_t, mask_dyn, t_phylo, t_adj)
    
    # reconstruction loss on corrupted entries only
    loss_rec <- nnf_mse_loss(pred[corrupt], t_X[corrupt])
    
    # optional smoothness penalty on hidden implied by adjacency:
    # (cheap proxy: encourage predictions to be close to graph-smoothed predictions)
    # This keeps things phylogenetically coherent without forcing strict linear BM.
    loss_smooth <- torch_tensor(0, device=device)
    if (cfg$lambda_smooth > 0) {
      pred_s <- torch_matmul(t_adj, pred)
      loss_smooth <- nnf_mse_loss(pred_s, pred)
    }
    
    loss <- loss_rec + cfg$lambda_smooth * loss_smooth
    loss$backward()
    
    # gradient clipping (stabilizes training)
    nn_utils_clip_grad_norm_(model$parameters, max_norm = cfg$clip_norm)
    
    opt$step()
    
    # --- validation every eval_every epochs
    if (epoch %% cfg$eval_every == 0) {
      model$eval()
      with_no_grad({
        # iterative refinement (few steps is enough)
        pred0 <- model(t_X, t_mask, t_phylo, t_adj)
        Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
        for (k in 1:cfg$refine_steps) {
          pk <- model(Xcurr, t_mask, t_phylo, t_adj)
          Xcurr <- t_X * t_mask + pk * (1 - t_mask)
        }
        pred_mat <- as.matrix(Xcurr$cpu())
      })
      
      val_rmse <- rmse_on_idx(pred_mat, val_idx)
      
      if (val_rmse + 1e-6 < best_val) {
        best_val <- val_rmse
        best_state <- model$state_dict()
        patience_left <- cfg$patience
      } else {
        patience_left <- patience_left - 1
      }
      
      if (cfg$verbose) {
        cat(sprintf("epoch %4d | loss %.5f | val RMSE %.4f | best %.4f | patience %d\n",
                    epoch, loss$item(), val_rmse, best_val, patience_left))
      }
      
      if (patience_left <= 0) break
    }
  }
  
  # restore best state
  if (!is.null(best_state)) model$load_state_dict(best_state)
  
  # final eval (val + test)
  model$eval()
  with_no_grad({
    pred0 <- model(t_X, t_mask, t_phylo, t_adj)
    Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
    for (k in 1:cfg$refine_steps) {
      pk <- model(Xcurr, t_mask, t_phylo, t_adj)
      Xcurr <- t_X * t_mask + pk * (1 - t_mask)
    }
    pred_mat <- as.matrix(Xcurr$cpu())
  })
  
  list(
    cfg = cfg,
    best_val = best_val,
    val_rmse = rmse_vec(X_truth[val_idx],  pred_mat[val_idx]),
    test_rmse = rmse_vec(X_truth[test_idx], pred_mat[test_idx]),
    pred = pred_mat
  )
}

cat("\n--- Tuning Phylo-DAE ---\n")

# small search (kept intentionally modest so it’s runnable)
grid <- list(
  list(hidden_dim=96,  dropout=0.10, lr=2e-3, weight_decay=1e-5, lambda_smooth=0.05),
  list(hidden_dim=128, dropout=0.10, lr=2e-3, weight_decay=1e-5, lambda_smooth=0.05),
  list(hidden_dim=128, dropout=0.15, lr=1e-3, weight_decay=2e-5, lambda_smooth=0.10),
  list(hidden_dim=192, dropout=0.10, lr=1e-3, weight_decay=1e-5, lambda_smooth=0.05)
)

base_cfg <- list(
  max_epochs   = 4000,
  patience     = 12,
  eval_every   = 100,
  corrupt_p    = 0.30,
  refine_steps = 8,
  clip_norm    = 1.0,
  verbose      = FALSE
)

results <- list()
for (g in seq_along(grid)) {
  cfg <- modifyList(base_cfg, grid[[g]])
  cat(sprintf("  [%d/%d] hidden=%d drop=%.2f lr=%.0e wd=%.0e smooth=%.2f ...\n",
              g, length(grid), cfg$hidden_dim, cfg$dropout, cfg$lr, cfg$weight_decay, cfg$lambda_smooth))
  # a couple seeds improves stability of selection
  r1 <- run_train_eval(cfg, seed=1)
  r2 <- run_train_eval(cfg, seed=2)
  results[[g]] <- list(
    cfg = cfg,
    val_rmse = mean(c(r1$val_rmse, r2$val_rmse)),
    test_rmse = mean(c(r1$test_rmse, r2$test_rmse)),
    best_run = if (r1$val_rmse <= r2$val_rmse) r1 else r2
  )
  cat(sprintf("      val RMSE ~ %.4f | test RMSE ~ %.4f\n", results[[g]]$val_rmse, results[[g]]$test_rmse))
}

best_g <- which.min(sapply(results, `[[`, "val_rmse"))
best <- results[[best_g]]$best_run
X_pred_dae <- best$pred

cat("\n=========================================\n")
cat(" RESULTS (RMSE on Standardized Log Data)\n")
cat("=========================================\n")
cat(sprintf("Rphylopars (BM) | val: %.4f | test: %.4f\n", rmse_val_rphylo, rmse_test_rphylo))
cat(sprintf("Phylo-DAE (tuned)| val: %.4f | test: %.4f\n", best$val_rmse, best$test_rmse))
cat("Best DAE config:\n")
print(best$cfg)
cat("=========================================\n")

# ----------------------------- visualization ---------------------------------

df_plot <- rbind(
  data.frame(
    Truth = X_truth[test_idx],
    Prediction = X_pred_rphylo[test_idx],
    Method = "Rphylopars (BM)"
  ),
  data.frame(
    Truth = X_truth[test_idx],
    Prediction = X_pred_dae[test_idx],
    Method = "Phylo-DAE (tuned)"
  )
)

p <- ggplot(df_plot, aes(x = Truth, y = Prediction)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_abline(linetype = "dashed") +
  facet_wrap(~Method, scales = "free") +
  theme_minimal() +
  labs(
    title = "Imputation on Held-out Cells (test split only)",
    subtitle = sprintf("Test RMSE: DAE=%.3f vs Rphylopars=%.3f",
                       best$test_rmse, rmse_test_rphylo)
  )
print(p)