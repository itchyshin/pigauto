# ==============================================================================
# PHYLO-DAE V3: THE CLADE-PRIORITY CHALLENGER (FULL SCRIPT)
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

# Match Tree and Data
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

# Log-transform and Standardize
X_raw   <- as.matrix(df_final[, trait_cols])
X_log   <- log(X_raw)
X_truth <- scale(X_log) 

# ── 2. Tree Feature Engineering (With Caching) ───────────────────────────────
cache_path <- here("avonet", "phylo_cache_v3.rds")

if (file.exists(cache_path)) {
  cat("Loading phylogenetic features from cache...\n")
  cache <- readRDS(cache_path)
  phylo_mat <- cache$phylo
  adj_mat   <- cache$adj
} else {
  cat("Computing tree features...\n")
  D <- cophenetic(tree_pruned)
  # Sharp Kernel: Focuses on local relatives (Clade Priority)
  sigma <- median(D) * 0.20 
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 1
  
  # Standard GCN Normalization
  rs <- rowSums(A) + 1e-8
  Dinv <- diag(1 / sqrt(rs))
  adj_mat <- Dinv %*% A %*% Dinv
  
  # Spectral Features (Higher resolution k=32)
  A_spec <- A; diag(A_spec) <- 0
  L <- diag(rowSums(A_spec)) - A_spec
  eig <- eigen(L, symmetric = TRUE)
  k <- 32 
  phylo_mat <- eig$vectors[, (nrow(A)-k):(nrow(A)-1)]
  saveRDS(list(phylo = phylo_mat, adj = adj_mat), cache_path)
}

t_phylo <- torch_tensor(as.matrix(phylo_mat), dtype=torch_float(), device=device)
t_adj   <- torch_tensor(as.matrix(adj_mat), dtype=torch_float(), device=device)

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

# ── 5. Clade-Priority Architecture ───────────────────────────────────────────
PhyloDAE_V3 <- nn_module(
  "PhyloDAE_V3",
  initialize = function(n_traits, n_phylo, hidden_dim = 512) {
    n_traits <- as.integer(n_traits)
    n_phylo  <- as.integer(n_phylo)
    hidden_dim <- as.integer(hidden_dim)
    
    # input: traits + mask + phylo
    input_dim <- as.integer(n_traits + n_traits + n_phylo)
    
    self$fc1 <- nn_linear(input_dim, hidden_dim)
    self$bn1 <- nn_batch_norm1d(hidden_dim)
    
    # GCN Gated Layers
    self$gcn1 <- nn_linear(hidden_dim, hidden_dim)
    self$gate <- nn_linear(hidden_dim, hidden_dim)
    
    self$out <- nn_linear(hidden_dim, n_traits)
    self$act <- nn_gelu()
    self$drop <- nn_dropout(0.05)
    
    # Explicitly float parameter to avoid the common dtype error
    self$mask_token <- nn_parameter(torch_randn(1, n_traits, dtype = torch_float()))
  },
  
  forward = function(x, m, p, adj) {
    # 1. Trait Injection
    x_input <- x * m + (1 - m) * self$mask_token
    h_init  <- torch_cat(list(x_input, m, p), dim = as.integer(2))
    
    # 2. Clade Aggregation (Early Message Passing)
    h_agg <- torch_matmul(adj, h_init)
    
    h <- self$fc1(h_agg) %>% self$bn1() %>% self$act()
    
    # 3. Residual Gating (Learns when to trust phylogeny vs. traits)
    h_res   <- h
    h_graph <- torch_matmul(adj, h)
    g       <- torch_sigmoid(self$gate(h))
    h       <- h_res + g * self$act(self$gcn1(h_graph))
    
    self$out(h)
  }
)

# ── 6. Training with Weighted Penalty ────────────────────────────────────────
X_fill <- as.matrix(X_with_NA); X_fill[is.na(X_fill)] <- 0
t_X            <- torch_tensor(X_fill, dtype = torch_float(), device = device)
t_mask         <- torch_tensor(as.matrix(mask_true), dtype = torch_float(), device = device)
t_truth_tensor <- torch_tensor(as.matrix(X_truth), dtype = torch_float(), device = device)

# Initialize with Dynamic Dimensions
model <- PhyloDAE_V3(
  n_traits   = ncol(X_truth), 
  n_phylo    = ncol(phylo_mat), 
  hidden_dim = 512
)$to(device = device)

optimizer <- optim_adamw(model$parameters, lr = 5e-4, weight_decay = 1e-2)
scheduler <- lr_one_cycle(optimizer, max_lr = 0.002, epochs = 10000, steps_per_epoch = 1)

# Heavily weight Mass (5x) to beat the linear baseline
trait_weights <- torch_tensor(c(5.0, 1.0, 1.0, 1.0), device = device)$view(as.integer(c(1, 4)))

cat("\n--- Training V3 (Clade Priority) ---\n")
for (epoch in 1:10000) {
  model$train()
  optimizer$zero_grad()
  
  u <- torch_rand_like(t_X)
  probs <- torch_tensor(c(0.3, 0.1, 0.1, 0.1), device = device)
  training_mask <- (u < probs) & (t_mask == 1)
  
  curr_input_mask <- t_mask * (1 - training_mask$to(torch_float()))
  preds <- model(t_X, curr_input_mask, t_phylo, t_adj)
  
  # Loss only on masked points with trait weights
  loss <- ((preds - t_truth_tensor)^2 * trait_weights)[training_mask]$mean()
  
  loss$backward()
  optimizer$step()
  scheduler$step()
  
  if (epoch %% 1000 == 0) cat(sprintf("Epoch %d | Loss: %.6f\n", epoch, as.numeric(loss$cpu())))
}

# ── 7. Inference & Final Refinement ───────────────────────────────────────────
cat("\nRefining Predictions...\n")
model$eval()
with_no_grad({
  final_preds <- model(t_X, t_mask, t_phylo, t_adj)
  X_pred_dae <- as.matrix(final_preds$cpu())
})

# ── 8. Results ────────────────────────────────────────────────────────────────
rmse_results <- data.frame(Trait = trait_cols, Rphylopars = 0, PhyloDAE = 0)
for(i in 1:length(trait_cols)) {
  col_mask   <- (col(X_truth) == i)
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