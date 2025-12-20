# ==============================================================================
# PHYLO-DAE V2: FULL BENCHMARK SCRIPT
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

# Load Data (Adjust paths as necessary)
tree  <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))

# Clean & Match
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)
trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])

common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

# Preprocessing
X_raw   <- as.matrix(df_final[, trait_cols])
X_log   <- log(X_raw)
X_truth <- scale(X_log) # Standardized for NN stability

# ── 2. Tree Feature Engineering (With Caching) ───────────────────────────────
cache_path <- here("avonet", "phylo_cache.rds")

if (file.exists(cache_path)) {
  cat("Loading phylogenetic features from cache...\n")
  cache <- readRDS(cache_path)
  phylo_mat <- cache$phylo
  adj_mat   <- cache$adj
} else {
  cat("Computing tree features (this may take a while)...\n")
  
  # A) Adjacency Matrix (Gaussian Kernel)
  D <- cophenetic(tree_pruned)
  sigma <- median(D) * 0.35
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 1
  rs <- rowSums(A) + 1e-8
  Dinv <- diag(1 / sqrt(rs))
  adj_mat <- Dinv %*% A %*% Dinv
  
  # B) Spectral Features (Laplacian Eigenmaps)
  A_spec <- A; diag(A_spec) <- 0
  L <- diag(rowSums(A_spec)) - A_spec
  eig <- eigen(L, symmetric = TRUE)
  k <- 16 
  phylo_mat <- eig$vectors[, (nrow(A)-k):(nrow(A)-1)]
  
  # Save for future tuning sessions
  saveRDS(list(phylo = phylo_mat, adj = adj_mat), cache_path)
  cat("Tree features cached to:", cache_path, "\n")
}

# Move to device
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
rmse_rphylo <- sqrt(mean((X_truth[missing_idx] - X_pred_rphylo[missing_idx])^2))
cat("Rphylopars RMSE:", round(rmse_rphylo, 4), "\n")

# ── 5. Phylo-DAE V2 Definition ────────────────────────────────────────────────
PhyloDAE_V2 <- nn_module(
  "PhyloDAE_V2",
  initialize = function(n_traits, n_phylo, hidden_dim = 256) {
    input_dim <- n_traits + n_traits + n_phylo 
    
    self$enc_in <- nn_linear(input_dim, hidden_dim)
    self$layer1 <- nn_linear(hidden_dim, hidden_dim)
    
    # Gated Message Passing (Captures evolutionary context)
    self$graph_weight <- nn_linear(hidden_dim, hidden_dim)
    self$graph_gate   <- nn_linear(hidden_dim, hidden_dim)
    
    self$dec_hidden <- nn_linear(hidden_dim * 2, hidden_dim)
    self$dec_out    <- nn_linear(hidden_dim, n_traits)
    
    self$act  <- nn_gelu()
    self$dropout <- nn_dropout(0.15)
    self$mask_token <- nn_parameter(torch_randn(1, n_traits))
  },
  
  forward = function(x, m, p, adj) {
    # Replace NAs with learnable mask token
    x_input <- x * m + (1 - m) * self$mask_token
    
    h_local <- self$enc_in(torch_cat(list(x_input, m, p), dim = 2)) %>% 
      self$act() %>% self$dropout()
    h_local <- self$layer1(h_local) %>% self$act()
    
    # Aggregate neighbors in phylogenetic space
    h_graph <- torch_matmul(adj, h_local)
    h_graph <- self$graph_weight(h_graph)
    
    # Gate decides reliance on phylogeny vs individual trait correlation
    gate <- torch_sigmoid(self$graph_gate(h_local))
    h_combined <- torch_cat(list(h_local, gate * h_graph), dim = 2)
    
    out <- self$dec_hidden(h_combined) %>% self$act() %>% self$dec_out()
    return(out)
  }
)

# ── 6. Training with OneCycleLR ───────────────────────────────────────────────
X_fill <- X_with_NA
X_fill[is.na(X_fill)] <- 0
t_X    <- torch_tensor(X_fill, dtype=torch_float(), device=device)
t_mask <- torch_tensor(mask_true, dtype=torch_float(), device=device)

model <- PhyloDAE_V2(n_traits=4, n_phylo=16, hidden_dim=512)$to(device=device)
optimizer <- optim_adamw(model$parameters, lr = 1e-3, weight_decay = 1e-4)

n_epochs <- 5000
scheduler <- lr_one_cycle(optimizer, max_lr = 0.005, epochs = n_epochs, steps_per_epoch = 1)

cat("\n--- Training Phylo-DAE V2 ---\n")
for (epoch in 1:n_epochs) {
  model$train()
  optimizer$zero_grad()
  
  # Dynamic Masking Strategy
  u <- torch_rand_like(t_X)
  training_mask <- (u < 0.2) & (t_mask == 1)
  current_input_mask <- t_mask * (1 - training_mask$to(dtype=torch_float()))
  
  preds <- model(t_X, current_input_mask, t_phylo, t_adj)
  loss  <- nnf_mse_loss(preds[training_mask], t_X[training_mask])
  
  loss$backward()
  optimizer$step()
  scheduler$step()
  
  if (epoch %% 500 == 0) cat(sprintf("Epoch %d | Loss: %.6f\n", epoch, loss$item()))
}

# ── 7. Inference & Refinement ─────────────────────────────────────────────────
cat("\nRefining Predictions...\n")
model$eval()
with_no_grad({
  final_preds <- model(t_X, t_mask, t_phylo, t_adj)
  X_refined <- t_X * t_mask + final_preds * (1 - t_mask)
  # Recurrent refinement loop
  for(i in 1:5) {
    final_preds <- model(X_refined, t_mask, t_phylo, t_adj)
    X_refined <- t_X * t_mask + final_preds * (1 - t_mask)
  }
  X_pred_dae <- as.matrix(X_refined$cpu())
})

# ── 8. Comparison ─────────────────────────────────────────────────────────────
rmse_dae <- sqrt(mean((X_truth[missing_idx] - X_pred_dae[missing_idx])^2))

cat("\n=========================================\n")
cat(sprintf("Rphylopars RMSE: %.4f\n", rmse_rphylo))
cat(sprintf("Phylo-DAE V2 RMSE: %.4f\n", rmse_dae))
cat("=========================================\n")

if (rmse_dae < rmse_rphylo) {
  cat(">> SUCCESS: Phylo-DAE outperformed Rphylopars! <<\n")
} else {
  cat(">> Rphylopars still leads. Increasing hidden_dim or epochs may help. <<\n")
}

