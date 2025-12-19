# ==============================================================================
# AVONET BENCHMARK: Rphylopars vs. Phylo-GNN-DAE (With Caching)
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(dplyr)
library(here)

# ── 1. Device & Data Loading ─────────────────────────────────────────────────
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

# Load Data
tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))

# Clean Names
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

# Select Traits (Multivariate)
trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])

# Match Tree
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

# Preprocessing (Log + Scale)
X_raw <- as.matrix(df_final[, trait_cols])
X_log <- log(X_raw)
X_mean <- colMeans(X_log)
X_sd   <- apply(X_log, 2, sd)
X_truth <- scale(X_log)

cat("Data Loaded: ", nrow(X_truth), "species x", ncol(X_truth), "traits.\n")

# ── 2. CACHING LOGIC: Tree Features ──────────────────────────────────────────
# We calculate these as R matrices first, save them, then convert to Tensors.

cache_file <- here("avonet", "phylo_cache.rds")

# Helper Functions (Modified to return R Matrices, not Tensors)
calc_adj_symnorm <- function(tree) {
  cat("  -> Calculating Cophenetic Distance (Slow)...\n")
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  cat("  -> Constructing Adjacency...\n")
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 1 # Self loops
  rs <- rowSums(A) + 1e-8
  Dinv_sqrt <- diag(1 / sqrt(rs))
  A_norm <- Dinv_sqrt %*% A %*% Dinv_sqrt
  return(as.matrix(A_norm))
}

calc_spectral_features <- function(tree, k = 8) {
  cat("  -> Calculating Laplacian Eigenmaps...\n")
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2)); diag(A) <- 0
  L <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric = TRUE)
  # Smallest non-zero eigenvectors
  coords <- eig$vectors[, (nrow(A)-k):(nrow(A)-1)] 
  return(as.matrix(coords))
}

if (file.exists(cache_file)) {
  cat("\n--- Loading Cached Phylogenetic Features ---\n")
  cached_data <- readRDS(cache_file)
  
  # Verify cache matches current tree size
  if (nrow(cached_data$phylo_emb) != nrow(df_final)) {
    stop("Cache file dimensions do not match current dataset! Delete 'phylo_cache.rds' and rerun.")
  }
  
  phylo_emb_mat <- cached_data$phylo_emb
  adj_mat       <- cached_data$adj
  cat("Loaded successfully from:", cache_file, "\n")
  
} else {
  cat("\n--- Computing Phylogenetic Features (First Run Only) ---\n")
  cat("This process involves large matrix decompositions and may take time...\n")
  
  phylo_emb_mat <- calc_spectral_features(tree_pruned, k = 8)
  adj_mat       <- calc_adj_symnorm(tree_pruned)
  
  # Save to disk
  saveRDS(list(phylo_emb = phylo_emb_mat, adj = adj_mat), cache_file)
  cat("Calculations complete. Saved to:", cache_file, "\n")
}

# Convert R Matrices to Torch Tensors
t_phylo <- torch_tensor(phylo_emb_mat, dtype = torch_float(), device = device)
t_adj   <- torch_tensor(adj_mat, dtype = torch_float(), device = device)


# ── 3. Artificial Missingness ────────────────────────────────────────────────
set.seed(123)
missing_frac <- 0.25
mask_true <- matrix(1, nrow=nrow(X_truth), ncol=ncol(X_truth))
missing_idx <- sample(length(X_truth), floor(missing_frac * length(X_truth)))
mask_true[missing_idx] <- 0
X_with_NA <- X_truth
X_with_NA[missing_idx] <- NA

# ── 4. Method 1: Rphylopars (Baseline) ───────────────────────────────────────
cat("\n--- Running Rphylopars ---\n")
df_rphylo <- data.frame(species = df_final$Species_Key, X_with_NA)
res_rphylo <- phylopars(trait_data = df_rphylo, tree = tree_pruned, model = "BM", pheno_error = FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]
rmse_rphylo <- sqrt(mean((X_truth[missing_idx] - X_pred_rphylo[missing_idx])^2))
cat("Rphylopars RMSE:", round(rmse_rphylo, 4), "\n")

# ── 5. Method 2: Phylo-DAE (GNN-DAE) ─────────────────────────────────────────
cat("\n--- Running Phylo-DAE (GNN-DAE) ---\n")

# Define Module
PhyloDAE <- nn_module(
  "PhyloDAE",
  initialize = function(n_traits, n_phylo, hidden_dim = 64) {
    input_dim <- n_traits + n_traits + n_phylo 
    self$enc1 <- nn_linear(input_dim, hidden_dim)
    self$enc2 <- nn_linear(hidden_dim, hidden_dim)
    self$msg   <- nn_linear(hidden_dim, hidden_dim)
    self$alpha <- nn_parameter(torch_tensor(0.1, device = device))
    self$dec1 <- nn_linear(hidden_dim, hidden_dim)
    self$dec2 <- nn_linear(hidden_dim, n_traits)
    self$act  <- nn_relu()
    self$drop <- nn_dropout(0.2)
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits, device = device))
  },
  forward = function(x, m, p, adj) {
    combined <- torch_cat(list(x, m, p), dim = 2)
    h <- self$enc1(combined); h <- self$act(h); h <- self$drop(h)
    h <- self$enc2(h);        h <- self$act(h)
    m_agg <- torch_matmul(adj, h)
    h <- h + self$alpha * self$msg(m_agg)
    h <- self$dec1(h); h <- self$act(h)
    self$dec2(h)
  }
)

# Prepare Tensors
X_fill <- X_with_NA
X_fill[is.na(X_fill)] <- 0
t_X          <- torch_tensor(X_fill, dtype = torch_float(), device = device)
t_mask_fixed <- torch_tensor(mask_true, dtype = torch_float(), device = device)

# Train Loop
model <- PhyloDAE(n_traits = ncol(X_truth), n_phylo = ncol(t_phylo), hidden_dim = 96)
model$to(device = device)
optimizer <- optim_adam(model$parameters, lr = 0.005)

cat("Training with Dynamic Corruption...\n")
n_epochs <- 1500
corruption_rate <- 0.3

for (i in 1:n_epochs) {
  model$train()
  optimizer$zero_grad()
  
  u <- torch_rand_like(t_X)
  is_corrupted <- (u < corruption_rate) & (t_mask_fixed == 1)
  
  if (as.numeric(is_corrupted$sum()$cpu()) > 0) {
    # Fix for R logical NOT: use (1 - float)
    is_corrupted_float <- is_corrupted$to(dtype = torch_float())
    mask_dynamic <- t_mask_fixed * (1 - is_corrupted_float)
    
    X_corrupted <- torch_where(is_corrupted, model$mask_token$expand_as(t_X), t_X)
    pred <- model(X_corrupted, mask_dynamic, t_phylo, t_adj)
    
    loss <- nnf_mse_loss(pred[is_corrupted], t_X[is_corrupted])
    loss$backward()
    optimizer$step()
    
    if (i %% 200 == 0) cat(sprintf("Epoch %d | Loss: %.4f\n", i, loss$item()))
  }
}

# Final Prediction
model$eval()
with_no_grad({
  final_pred <- model(t_X, t_mask_fixed, t_phylo, t_adj)
  X_pred_dae <- as.matrix(final_pred$cpu())
})

# ── 6. Comparison ────────────────────────────────────────────────────────────
rmse_dae <- sqrt(mean((X_truth[missing_idx] - X_pred_dae[missing_idx])^2))

cat("\n=========================================\n")
cat(" FINAL RESULTS (RMSE) \n")
cat("=========================================\n")
cat(sprintf("Rphylopars: %.4f\n", rmse_rphylo))
cat(sprintf("Phylo-DAE : %.4f\n", rmse_dae))
cat("=========================================\n")

# Visualization
df_plot <- data.frame(
  Truth = c(X_truth[missing_idx], X_truth[missing_idx]),
  Prediction = c(X_pred_rphylo[missing_idx], X_pred_dae[missing_idx]),
  Method = rep(c("Rphylopars", "Phylo-DAE"), each = length(missing_idx))
)

p <- ggplot(df_plot, aes(x = Truth, y = Prediction, color = Method)) +
  geom_point(alpha = 0.3) +
  geom_abline(color="black", linetype="dashed") +
  facet_wrap(~Method) +
  theme_minimal() +
  labs(title = "Imputation Performance: AVONET",
       subtitle = paste0("Phylo-DAE RMSE: ", round(rmse_dae, 3)))
print(p)