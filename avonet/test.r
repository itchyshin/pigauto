# ==============================================================================
# FAST Tree-DAE v2 + BM prior (AVONET benchmark)
# Key speedups:
#  - smaller model (hidden/K) by default
#  - RMSE computed in torch (no CPU matrix copy every eval)
#  - precompute node_ids (no arange in every forward)
#  - finite-grad check + nn_utils_clip_grad_norm_ (fast, no NA crash)
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
  m <- floor(missing_frac * length(all_idx))
  miss <- sample(all_idx, m)
  m_val <- floor(val_frac * length(miss))
  val_idx <- miss[seq_len(m_val)]
  test_idx <- miss[(m_val + 1):length(miss)]
  list(val_idx = val_idx, test_idx = test_idx)
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

# ----------------------------- device -----------------------------------------

device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")
# Optional: speed up CPU runs (tune to your machine)
# torch_set_num_threads(parallel::detectCores())

# ----------------------------- data -------------------------------------------

tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

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

t_X     <- torch_tensor(X_fill,  dtype=torch_float(), device=device)
t_mask  <- torch_tensor(mask_input, dtype=torch_float(), device=device)
t_truth <- torch_tensor(X_truth, dtype=torch_float(), device=device)

# boolean masks for fast RMSE-in-torch
val_mask_mat <- matrix(FALSE, nrow=n, ncol=p);  val_mask_mat[val_idx]  <- TRUE
test_mask_mat <- matrix(FALSE, nrow=n, ncol=p); test_mask_mat[test_idx] <- TRUE
t_val_mask  <- torch_tensor(val_mask_mat,  dtype=torch_bool(), device=device)
t_test_mask <- torch_tensor(test_mask_mat, dtype=torch_bool(), device=device)

# ----------------------------- baseline: Rphylopars ---------------------------

cat("\n--- Running Rphylopars (BM) ---\n")
df_rphylo <- data.frame(species = df_final$Species_Key, X_in)
res_rphylo <- phylopars(trait_data=df_rphylo, tree=tree_pruned, model="BM", pheno_error=FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]

rmse_val_rphylo  <- rmse_vec(X_truth[val_idx],  X_pred_rphylo[val_idx])
rmse_test_rphylo <- rmse_vec(X_truth[test_idx], X_pred_rphylo[test_idx])
cat(sprintf("Rphylopars RMSE | val: %.4f | test: %.4f\n", rmse_val_rphylo, rmse_test_rphylo))

# ----------------------------- tree edges -------------------------------------

N_all <- tree_pruned$Nnode + length(tree_pruned$tip.label)
stopifnot(N_all == max(tree_pruned$edge))

parent <- tree_pruned$edge[,1]
child  <- tree_pruned$edge[,2]
elen   <- tree_pruned$edge.length
if (is.null(elen)) stop("Tree has no branch lengths; BM penalty needs edge.length")

t_parent_all <- torch_tensor(parent, dtype=torch_long(), device=device)
t_child_all  <- torch_tensor(child,  dtype=torch_long(), device=device)
t_elen_all   <- torch_tensor(elen,   dtype=torch_float(), device=device)

tip_idx  <- torch_tensor(seq_len(n), dtype=torch_long(), device=device)
node_ids <- torch_arange(start=1, end=N_all + 1, dtype=torch_long(), device=device)  # precomputed

# ----------------------------- model ------------------------------------------

TreeDAE_fast <- nn_module(
  "TreeDAE_fast",
  initialize = function(n_traits, N_all, hidden=128, K=4, dropout=0.0) {
    self$n_traits <- n_traits
    self$N_all <- N_all
    self$K <- K
    
    self$in_fc <- nn_linear(n_traits + n_traits, hidden)
    self$in_ln <- nn_layer_norm(hidden)
    self$act   <- nn_gelu()
    self$drop  <- nn_dropout(dropout)
    
    self$node_emb <- nn_embedding(num_embeddings = N_all + 1, embedding_dim = hidden)
    self$node_scale <- nn_parameter(torch_tensor(0.05, dtype=torch_float()))
    
    self$msg_up   <- nn_linear(hidden, hidden)
    self$msg_down <- nn_linear(hidden, hidden)
    self$gate     <- nn_linear(hidden, hidden)
    self$lns      <- nn_module_list(lapply(seq_len(K), function(i) nn_layer_norm(hidden)))
    
    self$dec1 <- nn_linear(hidden, hidden)
    self$dec_ln1 <- nn_layer_norm(hidden)
    self$dec2 <- nn_linear(hidden, n_traits)
    
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits))
  },
  
  forward = function(x_tip, m_tip, parent_all, child_all, elen_all, tip_idx, node_ids) {
    
    # x_tip, m_tip are (n x p)
    z <- torch_cat(list(x_tip, m_tip), dim=2)
    h_tip <- self$in_fc(z)
    h_tip <- self$in_ln(h_tip)
    h_tip <- self$act(h_tip)
    h_tip <- self$drop(h_tip)
    
    # init all nodes from embeddings (learn internal nodes)
    h_all <- self$node_emb(node_ids) * self$node_scale
    
    # overwrite tips with data-driven embeddings (+ keep a bit of emb)
    h_all_tip <- h_all$index_select(1, tip_idx)
    h_all$index_copy_(dim=1, index=tip_idx, source=h_tip + h_all_tip)
    
    for (k in 1:self$K) {
      m_up   <- torch_zeros_like(h_all)
      m_down <- torch_zeros_like(h_all)
      
      src_p <- h_all$index_select(1, parent_all)
      m_down$index_add_(dim=1, index=child_all, source=self$msg_down(src_p))
      
      src_c <- h_all$index_select(1, child_all)
      m_up$index_add_(dim=1, index=parent_all, source=self$msg_up(src_c))
      
      m_total <- m_up + m_down
      gate <- torch_sigmoid(self$gate(h_all))
      
      h_all <- h_all + gate * m_total
      h_all <- self$lns[[k]](h_all)
      h_all <- self$act(h_all)
      h_all <- self$drop(h_all)
    }
    
    y_all <- self$dec1(h_all)
    y_all <- self$dec_ln1(y_all)
    y_all <- self$act(y_all)
    y_all <- self$dec2(y_all)
    
    y_tip <- y_all$index_select(1, tip_idx)
    list(y_tip=y_tip, y_all=y_all)
  }
)

bm_penalty_on_outputs <- function(y_all, parent_all, child_all, elen_all) {
  yp <- y_all$index_select(1, parent_all)
  yc <- y_all$index_select(1, child_all)
  elen_col <- elen_all$reshape(c(-1, 1))
  ((yc - yp)$pow(2) / elen_col)$mean()
}

rmse_torch <- function(pred, truth, mask_bool) {
  d <- (pred - truth)
  torch_sqrt(torch_mean(d[mask_bool]$pow(2)))
}

# ----------------------------- training --------------------------------------

run_tree_dae_fast <- function(
    hidden=128, K=4,
    lr=1e-3, weight_decay=1e-4,
    corrupt_p=0.85,
    lambda_bm=10,
    lambda_anchor=0.02,
    dropout=0.0,
    max_epochs=2500,          # << much smaller
    eval_every=200,           # << less frequent eval
    patience=10,
    clip_norm=1.0,
    seed=1,
    verbose=TRUE
) {
  set.seed(seed); torch_manual_seed(seed)
  
  model <- TreeDAE_fast(n_traits=p, N_all=N_all, hidden=hidden, K=K, dropout=dropout)
  model$to(device=device)
  opt <- optim_adamw(model$parameters, lr=lr, weight_decay=weight_decay)
  
  best_val <- Inf
  best_state <- NULL
  patience_left <- patience
  
  for (epoch in 1:max_epochs) {
    model$train()
    opt$zero_grad()
    
    # corruption mask (same shape as t_X)
    u <- torch_rand_like(t_X)
    corrupt <- (u < corrupt_p) & (t_mask == 1)
    if (as.numeric(corrupt$sum()$cpu()$item()) == 0) next
    corrupt_f <- corrupt$to(dtype=torch_float())
    
    X_in_t <- torch_where(corrupt, model$mask_token$to(device=device)$expand_as(t_X), t_X)
    mask_dyn <- t_mask * (1 - corrupt_f)
    
    out <- model(X_in_t, mask_dyn, t_parent_all, t_child_all, t_elen_all, tip_idx, node_ids)
    y_tip <- out$y_tip
    y_all <- out$y_all
    
    loss_rec <- nnf_mse_loss(y_tip[corrupt], t_X[corrupt])
    
    anchor_idx <- (mask_dyn == 1)   # observed AND not corrupted
    loss_anchor <- nnf_mse_loss(y_tip[anchor_idx], t_X[anchor_idx])
    
    loss_bm <- bm_penalty_on_outputs(y_all, t_parent_all, t_child_all, t_elen_all)
    
    loss <- loss_rec + lambda_anchor * loss_anchor + lambda_bm * loss_bm
    loss$backward()
    
    # Fast stability: if any grad non-finite, skip the step
    if (!all_grads_finite(model$parameters)) {
      opt$zero_grad()
      next
    }
    
    # Built-in clip is fast once grads are finite
    nn_utils_clip_grad_norm_(model$parameters, max_norm = clip_norm)
    opt$step()
    
    if (epoch %% eval_every == 0) {
      model$eval()
      with_no_grad({
        out0 <- model(t_X, t_mask, t_parent_all, t_child_all, t_elen_all, tip_idx, node_ids)
        pred0 <- out0$y_tip
        Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
        
        val_rmse_t <- rmse_torch(Xcurr, t_truth, t_val_mask)$cpu()$item()
        val_rmse <- as.numeric(val_rmse_t)
      })
      
      if (val_rmse + 1e-6 < best_val) {
        best_val <- val_rmse
        best_state <- model$state_dict()
        patience_left <- patience
      } else {
        patience_left <- patience_left - 1
      }
      
      if (verbose) {
        cat(sprintf("epoch %4d | rec %.4f | anc %.4f | bm %.4f | valRMSE %.4f | best %.4f | pat %d\n",
                    epoch, loss_rec$item(), loss_anchor$item(), loss_bm$item(),
                    val_rmse, best_val, patience_left))
      }
      
      if (patience_left <= 0) break
    }
  }
  
  if (!is.null(best_state)) model$load_state_dict(best_state)
  
  model$eval()
  with_no_grad({
    out0 <- model(t_X, t_mask, t_parent_all, t_child_all, t_elen_all, tip_idx, node_ids)
    pred0 <- out0$y_tip
    Xcurr <- t_X * t_mask + pred0 * (1 - t_mask)
    
    val_rmse <- as.numeric(rmse_torch(Xcurr, t_truth, t_val_mask)$cpu()$item())
    test_rmse <- as.numeric(rmse_torch(Xcurr, t_truth, t_test_mask)$cpu()$item())
    pred_mat <- as.matrix(Xcurr$cpu())
  })
  
  list(pred=pred_mat, val_rmse=val_rmse, test_rmse=test_rmse)
}

# ----------------------------- run (single fast config) -----------------------

cat("\n--- Running FAST Tree-DAE + BM prior ---\n")

# Start with this (fast). If you want more accuracy: increase hidden to 256 and K to 6.
best <- run_tree_dae_fast(
  hidden=128, K=4,
  lr=1e-3, corrupt_p=0.85,
  lambda_bm=10, lambda_anchor=0.02,
  max_epochs=2500, eval_every=200, patience=10,
  clip_norm=1.0, verbose=TRUE, seed=1
)

X_pred_tree_dae <- best$pred

cat("\n=========================================\n")
cat(" RESULTS (RMSE on Standardized Log Data)\n")
cat("=========================================\n")
cat(sprintf("Rphylopars (BM)        | val: %.4f | test: %.4f\n", rmse_val_rphylo, rmse_test_rphylo))
cat(sprintf("Tree-DAE fast (v2)     | val: %.4f | test: %.4f\n", best$val_rmse, best$test_rmse))
cat("=========================================\n")

# ----------------------------- visualization ---------------------------------

df_plot <- rbind(
  data.frame(Truth = X_truth[test_idx], Prediction = X_pred_rphylo[test_idx], Method = "Rphylopars (BM)"),
  data.frame(Truth = X_truth[test_idx], Prediction = X_pred_tree_dae[test_idx], Method = "Tree-DAE fast (v2)")
)

p_plot <- ggplot(df_plot, aes(x=Truth, y=Prediction)) +
  geom_point(alpha=0.25, size=1) +
  geom_abline(linetype="dashed") +
  facet_wrap(~Method, scales="free") +
  theme_minimal() +
  labs(
    title="Imputation on Held-out Cells (test split only)",
    subtitle=sprintf("Test RMSE: Tree-DAE=%.3f vs Rphylopars=%.3f", best$test_rmse, rmse_test_rphylo)
  )

print(p_plot)

##
# ==============================================================================
# Diagnostics: how close is Distilled Phylo-DAE to Rphylopars(BM)?
# + OU baseline
# + lambda_distill sweep
# ==============================================================================

cat("\n--- DIAGNOSTICS: DAE vs BM ---\n")

diff_mat <- X_pred_dae - X_pred_rphylo

summ_block <- function(name, idx=NULL) {
  if (is.null(idx)) {
    d <- as.numeric(diff_mat)
    t <- as.numeric(X_truth - X_pred_rphylo)   # BM error vs truth (all)
    td <- as.numeric(X_truth - X_pred_dae)     # DAE error vs truth (all)
  } else {
    d  <- as.numeric(diff_mat[idx])
    t  <- as.numeric(X_truth[idx] - X_pred_rphylo[idx])
    td <- as.numeric(X_truth[idx] - X_pred_dae[idx])
  }
  
  cat("\n", name, "\n", sep="")
  cat(sprintf("  |DAE - BM|: mean=%.6g  sd=%.6g  max=%.6g\n",
              mean(abs(d)), sd(d), max(abs(d))))
  cat(sprintf("  RMSE vs truth: BM=%.6g  DAE=%.6g  (DAE-BM)=%.6g\n",
              sqrt(mean(t^2)), sqrt(mean(td^2)), sqrt(mean((td - t)^2))))
  cat(sprintf("  Corr(pred_BM, pred_DAE)=%.6f\n",
              cor(as.numeric(if (is.null(idx)) X_pred_rphylo else X_pred_rphylo[idx]),
                  as.numeric(if (is.null(idx)) X_pred_dae   else X_pred_dae[idx]))))
  invisible(NULL)
}

summ_block("ALL CELLS")
summ_block("VAL CELLS",  val_idx)
summ_block("TEST CELLS", test_idx)

# How often are they *numerically identical* on held-out cells?
tol <- 1e-6
ident_val  <- mean(abs(diff_mat[val_idx])  < tol)
ident_test <- mean(abs(diff_mat[test_idx]) < tol)
cat(sprintf("\nFraction identical within tol=%.0e: val=%.3f  test=%.3f\n", tol, ident_val, ident_test))

# Quick residual plot (DAE - BM) on held-out cells
df_res <- rbind(
  data.frame(split="val",  resid = as.numeric(diff_mat[val_idx])),
  data.frame(split="test", resid = as.numeric(diff_mat[test_idx]))
)

p_res <- ggplot(df_res, aes(x=resid)) +
  geom_histogram(bins=60) +
  facet_wrap(~split, scales="free_y") +
  theme_minimal() +
  labs(title="Residual correction learned by DAE (DAE - BM)",
       subtitle="If this is tightly centered at 0, DAE is essentially emulating BM")
print(p_res)

# ==============================================================================
# Baseline check: does OU improve over BM in Rphylopars?
# ==============================================================================

cat("\n--- RPHYL0PARS: BM vs OU ---\n")

res_rphylo_OU <- phylopars(trait_data=df_rphylo, tree=tree_pruned, model="OU", pheno_error=FALSE)
X_pred_rphylo_OU <- res_rphylo_OU$anc_recon[df_rphylo$species, trait_cols]

rmse_val_OU  <- rmse_vec(X_truth[val_idx],  X_pred_rphylo_OU[val_idx])
rmse_test_OU <- rmse_vec(X_truth[test_idx], X_pred_rphylo_OU[test_idx])

cat(sprintf("BM | val %.4f | test %.4f\n", rmse_val_rphylo, rmse_test_rphylo))
cat(sprintf("OU | val %.4f | test %.4f\n", rmse_val_OU, rmse_test_OU))

# ==============================================================================
# Sweep lambda_distill (does relaxing distillation help or hurt?)
# - keep other settings fixed (adjust max_epochs down to keep it quick)
# ==============================================================================

cat("\n--- LAMBDA DISTILL SWEEP ---\n")

lambdas <- c(0.50, 0.20, 0.05, 0.00)

sweep <- lapply(lambdas, function(lam) {
  cat(sprintf("lambda_distill=%.2f ... ", lam))
  r <- run_distilled_phylo_dae(
    hidden=192, K=4, dropout=0.10,
    lr=2e-3, weight_decay=1e-5,
    corrupt_p=0.30,
    lambda_distill=lam,
    lambda_anchor=0.05,
    max_epochs=1200, eval_every=100, patience=8,
    refine_steps=6,
    seed=1, verbose=FALSE
  )
  cat(sprintf("val %.4f | test %.4f\n", r$val_rmse, r$test_rmse))
  data.frame(lambda_distill=lam, val_rmse=r$val_rmse, test_rmse=r$test_rmse)
})

sweep_df <- do.call(rbind, sweep)
print(sweep_df)

p_sweep <- ggplot(sweep_df, aes(x=lambda_distill, y=test_rmse)) +
  geom_point(size=2) + geom_line() +
  theme_minimal() +
  labs(title="Effect of distillation strength on test RMSE",
       subtitle="If RMSE worsens as lambda_distill decreases, BM already saturates the signal")
print(p_sweep)