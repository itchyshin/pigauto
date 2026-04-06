# ==============================================================================
# IMPROVED: Distilled Tree-DAE v2 (Knowledge Distillation Mode)
# Goal: Beat Rphylopars by using it as a "Teacher"
# ==============================================================================

library(torch)
library(ape)
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
  miss <- sample(all_idx, floor(missing_frac * length(all_idx)))
  m_val <- floor(val_frac * length(miss))
  val_idx <- miss[seq_len(m_val)]
  test_idx <- miss[(m_val + 1):length(miss)]
  list(val_idx = val_idx, test_idx = test_idx)
}

all_grads_finite <- function(parameters) {
  for (pp in parameters) {
    g <- pp$grad
    if (!is.null(g)) {
      if (!isTRUE(torch_isfinite(g)$all()$item())) return(FALSE)
    }
  }
  TRUE
}

# ----------------------------- device -----------------------------------------
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

# ----------------------------- data -------------------------------------------
tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

X_raw   <- as.matrix(df_final[, trait_cols])
X_truth <- scale(log(X_raw))
n <- nrow(X_truth); p <- ncol(X_truth)

# ----------------------------- missingness ------------------------------------
splits <- make_missing_splits(X_truth)
val_idx  <- splits$val_idx
test_idx <- splits$test_idx

mask_input <- matrix(1, nrow=n, ncol=p)
mask_input[c(val_idx, test_idx)] <- 0

X_in <- X_truth
X_in[c(val_idx, test_idx)] <- NA
X_fill <- X_in
X_fill[is.na(X_fill)] <- 0

t_X     <- torch_tensor(X_fill,  dtype=torch_float(), device=device)
t_mask  <- torch_tensor(mask_input, dtype=torch_float(), device=device)
t_truth <- torch_tensor(X_truth, dtype=torch_float(), device=device)

# masks for RMSE
t_val_mask  <- torch_tensor(matrix(FALSE, n, p) %>% {.[val_idx] <- TRUE; .}, dtype=torch_bool(), device=device)
t_test_mask <- torch_tensor(matrix(FALSE, n, p) %>% {.[test_idx] <- TRUE; .}, dtype=torch_bool(), device=device)

# ----------------------------- baseline: Rphylopars (The Teacher) -------------
cat("\n--- Running Rphylopars (Teacher BM) ---\n")
df_rphylo <- data.frame(species = df_final$Species_Key, X_in)
res_rphylo <- phylopars(trait_data=df_rphylo, tree=tree_pruned, model="BM", pheno_error=FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]

rmse_test_rphylo <- rmse_vec(X_truth[test_idx], X_pred_rphylo[test_idx])
cat(sprintf("Rphylopars RMSE | test: %.4f\n", rmse_test_rphylo))

# Convert BM prediction to tensor for distillation
t_bm <- torch_tensor(as.matrix(X_pred_rphylo), dtype=torch_float(), device=device)

# ----------------------------- model ------------------------------------------
TreeDAE_distill <- nn_module(
  "TreeDAE_distill",
  initialize = function(n_traits, N_all, hidden=192, K=4) {
    self$K <- K
    self$in_fc <- nn_linear(n_traits * 2, hidden)
    self$node_emb <- nn_embedding(N_all + 1, hidden)
    self$msg_up   <- nn_linear(hidden, hidden)
    self$msg_down <- nn_linear(hidden, hidden)
    self$gate     <- nn_linear(hidden, hidden)
    self$lns      <- nn_module_list(lapply(1:K, function(i) nn_layer_norm(hidden)))
    self$dec      <- nn_sequential(nn_linear(hidden, hidden), nn_gelu(), nn_linear(hidden, n_traits))
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits))
  },
  forward = function(x, m, parent_all, child_all, tip_idx, node_ids) {
    h_tip <- self$in_fc(torch_cat(list(x, m), 2))
    h_all <- self$node_emb(node_ids) * 0.05
    h_all$index_copy_(1, tip_idx, h_tip + h_all$index_select(1, tip_idx))
    
    for (k in 1:self$K) {
      m_up <- torch_zeros_like(h_all); m_down <- torch_zeros_like(h_all)
      m_down$index_add_(1, child_all, self$msg_down(h_all$index_select(1, parent_all)))
      m_up$index_add_(1, parent_all, self$msg_up(h_all$index_select(1, child_all)))
      gate <- torch_sigmoid(self$gate(h_all))
      h_all <- self$lns[[k]](h_all + gate * (m_up + m_down))
    }
    self$dec(h_all$index_select(1, tip_idx))
  }
)

# ----------------------------- training function ------------------------------

run_distilled_phylo_dae <- function(hidden=192, K=4, lr=5e-4, lambda_bm=5, lambda_distill=0.5, refine_steps=5, seed=1) {
  set.seed(seed); torch_manual_seed(seed)
  
  # Tree setup
  N_all <- tree_pruned$Nnode + length(tree_pruned$tip.label)
  t_parent <- torch_tensor(tree_pruned$edge[,1], dtype=torch_long(), device=device)
  t_child  <- torch_tensor(tree_pruned$edge[,2], dtype=torch_long(), device=device)
  tip_idx  <- torch_tensor(seq_len(n), dtype=torch_long(), device=device)
  node_ids <- torch_arange(1, N_all + 1, dtype=torch_long(), device=device)
  
  model <- TreeDAE_distill(p, N_all, hidden, K)$to(device=device)
  opt <- optim_adamw(model$parameters, lr=lr)
  
  best_val <- Inf; best_state <- NULL; patience <- 10
  
  for (epoch in 1:2000) {
    model$train(); opt$zero_grad()
    
    # Denoising corruption
    corrupt <- (torch_rand_like(t_X) < 0.3) & (t_mask == 1)
    X_in_t <- torch_where(corrupt, model$mask_token$expand_as(t_X), t_X)
    
    pred <- model(X_in_t, t_mask * (1 - corrupt$to(torch_float())), t_parent, t_child, tip_idx, node_ids)
    
    # Losses
    loss_rec <- nnf_mse_loss(pred[corrupt], t_X[corrupt])
    # Distillation: learn from BM teacher on missing cells
    loss_distill <- nnf_mse_loss(pred[t_mask == 0], t_bm[t_mask == 0])
    
    loss <- loss_rec + (lambda_distill * loss_distill)
    loss$backward(); nn_utils_clip_grad_norm_(model$parameters, 1.0); opt$step()
    
    if (epoch %% 100 == 0) {
      model$eval(); with_no_grad({
        # Inference with Refinement
        p0 <- model(t_X, t_mask, t_parent, t_child, tip_idx, node_ids)
        Xcurr <- t_X * t_mask + p0 * (1 - t_mask)
        for(r in 1:refine_steps) {
          p0 <- model(Xcurr, t_mask, t_parent, t_child, tip_idx, node_ids)
          Xcurr <- t_X * t_mask + p0 * (1 - t_mask)
        }
        val_rmse <- torch_sqrt(torch_mean((Xcurr[t_val_mask] - t_truth[t_val_mask])^2))$item()
      })
      if (val_rmse < best_val) {
        best_val <- val_rmse; best_state <- model$state_dict(); patience <- 10
      } else {
        patience <- patience - 1
      }
      if (patience <= 0) break
    }
  }
  model$load_state_dict(best_state)
  return(list(model=model, val_rmse=best_val))
}

# ----------------------------- execution & comparison -------------------------

cat("\n--- Running Distilled Tree-DAE ---\n")
results <- run_distilled_phylo_dae(lambda_distill=0.25, refine_steps=3)

# Final Inference
results$model$eval()
with_no_grad({
  N_all <- tree_pruned$Nnode + length(tree_pruned$tip.label)
  t_parent <- torch_tensor(tree_pruned$edge[,1], dtype=torch_long(), device=device)
  t_child  <- torch_tensor(tree_pruned$edge[,2], dtype=torch_long(), device=device)
  tip_idx  <- torch_tensor(seq_len(n), dtype=torch_long(), device=device)
  node_ids <- torch_arange(1, N_all + 1, dtype=torch_long(), device=device)
  
  X_pred_dae_tensor <- results$model(t_X, t_mask, t_parent, t_child, tip_idx, node_ids)
  X_pred_dae <- as.matrix(X_pred_dae_tensor$cpu())
})

rmse_test_dae <- rmse_vec(X_truth[test_idx], X_pred_dae[test_idx])

cat("\n=========================================\n")
cat(sprintf("FINAL TEST RMSE (BM):  %.4f\n", rmse_test_rphylo))
cat(sprintf("FINAL TEST RMSE (DAE): %.4f\n", rmse_test_dae))
cat("=========================================\n")

# Residual Diagnostic Plot
diff_mat <- X_pred_dae - X_pred_rphylo
df_res <- data.frame(resid = as.numeric(diff_mat[test_idx]))
ggplot(df_res, aes(x=resid)) + 
  geom_histogram(bins=50, fill="steelblue", alpha=0.7) +
  theme_minimal() + 
  labs(title="Residual Correction (DAE - BM)", subtitle="Ideally centered at 0 with small variance")