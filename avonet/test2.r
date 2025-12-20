# ==============================================================================
# PHYLO-DAE CHALLENGER: RESIDUAL GCN + WEIGHTED LOSS
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(dplyr)
library(here)

# ── 1. Device & Data Setup ────────────────────────────────────────────────────
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

# Load Data
tree  <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])

common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

# Preprocessing
X_raw   <- as.matrix(df_final[, trait_cols])
X_log   <- log(X_raw)
X_truth <- scale(X_log) 

# ── 2. Tree Feature Engineering (With Caching) ───────────────────────────────
cache_path <- here("avonet", "phylo_cache.rds")

if (file.exists(cache_path)) {
  cat("Loading phylogenetic features from cache...\n")
  cache <- readRDS(cache_path)
  phylo_mat <- cache$phylo
  adj_mat   <- cache$adj
} else {
  cat("Computing tree features...\n")
  D <- cophenetic(tree_pruned)
  sigma <- median(D) * 0.35
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 1
  rs <- rowSums(A) + 1e-8
  Dinv <- diag(1 / sqrt(rs))
  adj_mat <- Dinv %*% A %*% Dinv
  
  A_spec <- A; diag(A_spec) <- 0
  L <- diag(rowSums(A_spec)) - A_spec
  eig <- eigen(L, symmetric = TRUE)
  k <- 16 
  phylo_mat <- eig$vectors[, (nrow(A)-k):(nrow(A)-1)]
  saveRDS(list(phylo = phylo_mat, adj = adj_mat), cache_path)
}

t_phylo <- torch_tensor(phylo_mat, dtype=torch_float(), device=device)
t_adj   <- torch_tensor(adj_mat, dtype=torch_float(), device=device)

# ── 3. Artificial Missingness (25%) ───────────────────────────────────────────
set.seed(555)
missing_frac <- 0.25
mask_true <- matrix(1, nrow=nrow(X_truth), ncol=ncol(X_truth))
missing_idx <- sample(length(X_truth), floor(missing_frac * length(X_truth)))
mask_true[missing_idx] <- 0
X_with_NA <- X_truth
X_with_NA[missing_idx] <- NA

# ── 4. Baseline: Rphylopars ───────────────────────────────────────────────────
cat("\n--- Running Rphylopars (ML Baseline) ---\n")
df_rphylo <- data.frame(species = df_final$Species_Key, X_with_NA)
res_rphylo <- phylopars(trait_data = df_rphylo, tree = tree_pruned, model = "BM")
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]
rmse_rphylo_total <- sqrt(mean((X_truth[missing_idx] - X_pred_rphylo[missing_idx])^2))
cat("Rphylopars Total RMSE:", round(rmse_rphylo_total, 4), "\n")

# ── 5. Enhanced Phylo-DAE Challenger (Bulletproof Version) ──────────────────
PhyloDAE_Challenger <- nn_module(
  "PhyloDAE_Challenger",
  initialize = function(n_traits, n_phylo, hidden_dim = 512) {
    
    # Force cast to integers for R and ensure dimensions are clean
    n_traits <- as.integer(n_traits)
    n_phylo  <- as.integer(n_phylo)
    hidden_dim <- as.integer(hidden_dim)
    
    # Input: Traits (4) + Mask (4) + Phylo (16)
    input_dim <- n_traits + n_traits + n_phylo
    
    self$fc1 <- nn_linear(input_dim, hidden_dim)
    self$bn1 <- nn_batch_norm1d(hidden_dim)
    
    # GCN Layers with Residuals
    self$gcn1 <- nn_linear(hidden_dim, hidden_dim)
    self$gcn2 <- nn_linear(hidden_dim, hidden_dim)
    
    self$fc2 <- nn_linear(hidden_dim, hidden_dim)
    self$out <- nn_linear(hidden_dim, n_traits)
    
    self$act  <- nn_gelu()
    self$drop <- nn_dropout(0.1)
    
    # Explicitly set dtype for the parameter to avoid the floating_point error
    self$mask_token <- nn_parameter(torch_randn(1, n_traits, dtype = torch_float()))
  },
  
  forward = function(x, m, p, adj) {
    # 1. Input Layer
    x_input <- x * m + (1 - m) * self$mask_token
    h <- torch_cat(list(x_input, m, p), dim = 2)
    
    h <- self$fc1(h) %>% self$bn1() %>% self$act() %>% self$drop()
    
    # 2. Residual Graph Convolutions
    # Force the neighbor message passing
    h_graph <- torch_matmul(adj, h)
    h <- h + self$act(self$gcn1(h_graph)) 
    
    h_graph2 <- torch_matmul(adj, h)
    h <- h + self$act(self$gcn2(h_graph2)) 
    
    # 3. Output
    h <- self$fc2(h) %>% self$act()
    self$out(h)
  }
)

# ── 6. Training with Explicit Type Casting ──────────────────────────────

# Explicitly clean the matrices
X_fill <- as.matrix(X_with_NA)
X_fill[is.na(X_fill)] <- 0

# Ensure device and dtype are strictly set
t_X            <- torch_tensor(X_fill, dtype = torch_float(), device = device)
t_mask         <- torch_tensor(as.matrix(mask_true), dtype = torch_float(), device = device)
t_truth_tensor <- torch_tensor(as.matrix(X_truth), dtype = torch_float(), device = device)
t_phylo        <- torch_tensor(as.matrix(phylo_mat), dtype = torch_float(), device = device)
t_adj          <- torch_tensor(as.matrix(adj_mat), dtype = torch_float(), device = device)

# Initialize the model using explicit integers
model <- PhyloDAE_Challenger(
  n_traits   = as.integer(4), 
  n_phylo    = as.integer(16), 
  hidden_dim = as.integer(512)
)$to(device = device)

optimizer <- optim_adamw(model$parameters, lr = 1e-3, weight_decay = 1e-3)
scheduler <- lr_one_cycle(optimizer, max_lr = 0.005, epochs = 10000, steps_per_epoch = 1)

# Trait weights: Mass (Col 1) is 3x importance
trait_weights <- torch_tensor(c(3.0, 1.0, 1.0, 1.0), device = device)$view(c(1, 4))

cat("\n--- Training Challenger Edition (Fixed Dtype) ---\n")

for (epoch in 1:10000) {
  model$train()
  optimizer$zero_grad()
  
  # Strategic Masking
  u <- torch_rand_like(t_X)
  probs <- torch_tensor(c(0.3, 0.1, 0.1, 0.1), device = device)
  
  # Ensure mask is a pure logical for indexing
  training_mask <- (u < probs) & (t_mask == 1)
  
  # Cast for the forward pass
  curr_input_mask <- t_mask * (1 - training_mask$to(torch_float()))
  
  preds <- model(t_X, curr_input_mask, t_phylo, t_adj)
  
  # Weighted Loss Logic
  diff_sq <- (preds - t_truth_tensor)^2
  weighted_diff <- diff_sq * trait_weights
  
  # Subsetting with a boolean mask
  loss <- weighted_diff[training_mask]$mean()
  
  loss$backward()
  optimizer$step()
  scheduler$step()
  
  if (epoch %% 1000 == 0) {
    cat(sprintf("Epoch %d | Loss: %.6f\n", epoch, as.numeric(loss$cpu())))
  }
}
# ── 7. Inference & Refinement ─────────────────────────────────────────────────
cat("\nRefining Predictions...\n")
model$eval()
with_no_grad({
  final_preds <- model(t_X, t_mask, t_phylo, t_adj)
  X_curr <- t_X * t_mask + final_preds * (1 - t_mask)
  for(i in 1:10) { # Increased refinement steps
    final_preds <- model(X_curr, t_mask, t_phylo, t_adj)
    X_curr <- t_X * t_mask + final_preds * (1 - t_mask)
  }
  X_pred_dae <- as.matrix(X_curr$cpu())
})

# ── 8. Detailed Comparison ────────────────────────────────────────────────────
# Calculate per-trait RMSE
rmse_results <- data.frame(Trait = trait_cols, Rphylopars = 0, PhyloDAE = 0)

for(i in 1:length(trait_cols)) {
  col_mask <- (row(X_truth) %in% 1:nrow(X_truth)) & (col(X_truth) == i)
  target_idx <- intersect(missing_idx, which(col_mask))
  
  rmse_results$Rphylopars[i] <- sqrt(mean((X_truth[target_idx] - X_pred_rphylo[target_idx])^2))
  rmse_results$PhyloDAE[i]   <- sqrt(mean((X_truth[target_idx] - X_pred_dae[target_idx])^2))
}

cat("\n=========================================\n")
print(rmse_results)
cat("=========================================\n")
cat(sprintf("Rphylopars Total RMSE: %.4f\n", rmse_rphylo_total))
cat(sprintf("Phylo-DAE Total RMSE : %.4f\n", sqrt(mean((X_truth[missing_idx] - X_pred_dae[missing_idx])^2))))
cat("=========================================\n")