# ==============================================================================
# AVONET SOTA: Structure-Aware Graph Transformer (Phylo-Attention)
# ==============================================================================
# Key Innovation: Replaces fixed GCN with a Graph Transformer that uses
# the phylogenetic distance as an "Attention Bias".
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(dplyr)
library(here)

# ── 1. Setup ─────────────────────────────────────────────────────────────────
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

# Load Data
tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

# Select Continuous Traits
trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])

# Match Tree
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

# Preprocessing
X_raw <- as.matrix(df_final[, trait_cols])
X_log <- log(X_raw)
X_truth <- scale(X_log) # Standardized Log-Traits

cat("Data:", nrow(X_truth), "species. Traits:", colnames(X_truth), "\n")

# ── 2. Phylogenetic Attention Bias (The "Structure" in Transformer) ──────────
# Instead of a hard adjacency, we create a "Distance Bias" matrix.
# Closer species = Higher Bias.

cache_file <- here("avonet", "phylo_attn_cache.rds")

if (file.exists(cache_file)) {
  cat("Loading cached tree features...\n")
  cache <- readRDS(cache_file)
  if(nrow(cache$phylo) == nrow(df_final)) {
    t_phylo <- torch_tensor(cache$phylo, dtype=torch_float(), device=device)
    t_attn_bias <- torch_tensor(cache$attn_bias, dtype=torch_float(), device=device)
  } else { stop("Cache mismatch. Delete .rds and rerun.") }
} else {
  cat("Calculating Tree Attention Bias (may take a moment)...\n")
  
  # 1. Cophenetic Distance
  D <- cophenetic(tree_pruned)
  
  # 2. Attention Bias: Log-Space Gaussian Kernel
  # Transformer Attention is Softmax(Scores + Bias).
  # We want Bias ~ 0 for close relatives, -Inf for distant ones.
  # We use log( exp(-D^2 / sigma) ) = -D^2 / sigma
  sigma <- median(D) * 0.5
  Attn_Bias <- - (D^2) / (2 * sigma^2)
  
  # 3. Spectral Embeddings (Positional Encoding)
  # Still useful as "GPS" coordinates for the nodes
  A <- exp(Attn_Bias); diag(A) <- 0
  L <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric=TRUE)
  k <- 24 # More eigenvectors for transformer
  phylo_mat <- eig$vectors[, (nrow(A)-k):(nrow(A)-1)]
  
  saveRDS(list(phylo=phylo_mat, attn_bias=Attn_Bias), cache_file)
  
  t_phylo <- torch_tensor(phylo_mat, dtype=torch_float(), device=device)
  t_attn_bias <- torch_tensor(Attn_Bias, dtype=torch_float(), device=device)
}

# ── 3. Artificial Missingness ────────────────────────────────────────────────
set.seed(42)
missing_frac <- 0.25
mask_true <- matrix(1, nrow=nrow(X_truth), ncol=ncol(X_truth))
missing_idx <- sample(length(X_truth), floor(missing_frac * length(X_truth)))
mask_true[missing_idx] <- 0
X_with_NA <- X_truth
X_with_NA[missing_idx] <- NA

# ── 4. Baseline: Rphylopars ──────────────────────────────────────────────────
cat("\n--- Running Rphylopars (Baseline) ---\n")
df_rphylo <- data.frame(species = df_final$Species_Key, X_with_NA)
# BM is the theoretical limit if data is purely linear Brownian Motion
res_rphylo <- phylopars(trait_data = df_rphylo, tree = tree_pruned, model = "BM", pheno_error = FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]
rmse_rphylo <- sqrt(mean((X_truth[missing_idx] - X_pred_rphylo[missing_idx])^2))
cat("Rphylopars RMSE:", round(rmse_rphylo, 4), "\n")

# ── 5. The Model: Graph Transformer (Phylo-Attention) ────────────────────────
cat("\n--- Running Phylo-Graph-Transformer ---\n")

# Single Head Attention Layer with Structural Bias
GraphAttention <- nn_module(
  "GraphAttention",
  initialize = function(in_dim, out_dim) {
    self$query <- nn_linear(in_dim, out_dim)
    self$key   <- nn_linear(in_dim, out_dim)
    self$value <- nn_linear(in_dim, out_dim)
    self$scale <- 1 / sqrt(out_dim)
  },
  forward = function(h, attn_bias) {
    # h: [N, Dim]
    Q <- self$query(h)
    K <- self$key(h)
    V <- self$value(h)
    
    # Standard Self-Attention Score: Q @ K^T
    scores <- torch_matmul(Q, K$transpose(1, 2)) * self$scale
    
    # Inject Tree Structure: Add Bias
    # Nodes that are far apart in tree get a large negative number (effectively masked)
    # But the model can LEARN to overcome this if needed.
    scores <- scores + attn_bias
    
    attn_weights <- nnf_softmax(scores, dim = 2)
    
    out <- torch_matmul(attn_weights, V)
    return(out)
  }
)

PhyloTransformer <- nn_module(
  "PhyloTransformer",
  initialize = function(n_traits, n_phylo, hidden_dim = 128) {
    input_dim <- n_traits + n_traits + n_phylo
    
    # Encoder
    self$enc_in <- nn_linear(input_dim, hidden_dim)
    self$act    <- nn_gelu()
    self$drop   <- nn_dropout(0.1)
    
    # Transformer Blocks
    self$gat1 <- GraphAttention(hidden_dim, hidden_dim)
    self$ln1  <- nn_layer_norm(hidden_dim)
    self$ff1  <- nn_sequential(nn_linear(hidden_dim, hidden_dim*2), nn_gelu(), nn_linear(hidden_dim*2, hidden_dim))
    
    self$gat2 <- GraphAttention(hidden_dim, hidden_dim)
    self$ln2  <- nn_layer_norm(hidden_dim)
    self$ff2  <- nn_sequential(nn_linear(hidden_dim, hidden_dim*2), nn_gelu(), nn_linear(hidden_dim*2, hidden_dim))
    
    # Decoder
    self$head <- nn_linear(hidden_dim, n_traits)
    
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits, device=device))
  },
  
  forward = function(x, m, p, attn_bias) {
    combined <- torch_cat(list(x, m, p), dim=2)
    h <- self$enc_in(combined)
    h <- self$act(h); h <- self$drop(h)
    
    # Block 1 (Attention + Skip + Norm + FF)
    h_attn <- self$gat1(h, attn_bias)
    h <- self$ln1(h + h_attn) # Skip connection
    h_ff <- self$ff1(h)
    h <- self$ln1(h + h_ff)   # Skip connection
    
    # Block 2
    h_attn <- self$gat2(h, attn_bias)
    h <- self$ln2(h + h_attn)
    h_ff <- self$ff2(h)
    h <- self$ln2(h + h_ff)
    
    out <- self$head(h)
    return(out)
  }
)

# Setup
X_fill <- X_with_NA; X_fill[is.na(X_fill)] <- 0
t_X <- torch_tensor(X_fill, dtype=torch_float(), device=device)
t_mask <- torch_tensor(mask_true, dtype=torch_float(), device=device)

model <- PhyloTransformer(n_traits=4, n_phylo=24, hidden_dim=192)
model$to(device=device)
# AdamW + Weight Decay is crucial for Transformers
optimizer <- optim_adamw(model$parameters, lr=0.001, weight_decay=1e-4)

# Training with Self-Supervised Masking
cat("Training Graph Transformer...\n")
n_epochs <- 2000

for (i in 1:n_epochs) {
  model$train()
  optimizer$zero_grad()
  
  # Dynamic Corruption
  u <- torch_rand_like(t_X)
  is_corrupted <- (u < 0.30) & (t_mask == 1)
  
  if (as.numeric(is_corrupted$sum()$cpu()) > 0) {
    is_corrupted_f <- is_corrupted$to(dtype=torch_float())
    
    # Input: Corrupted replaced by Token
    X_in <- torch_where(is_corrupted, model$mask_token$expand_as(t_X), t_X)
    
    # Mask: 0 for Corrupted (telling model "Predict This")
    mask_dynamic <- t_mask * (1 - is_corrupted_f)
    
    pred <- model(X_in, mask_dynamic, t_phylo, t_attn_bias)
    
    loss <- nnf_mse_loss(pred[is_corrupted], t_X[is_corrupted])
    loss$backward()
    optimizer$step()
    
    if(i %% 200 == 0) cat(sprintf("Epoch %d | Loss: %.5f\n", i, loss$item()))
  }
}

# ── 6. Refinement Inference (The "Diffusion" Trick) ──────────────────────────
cat("\nRunning Refinement Inference...\n")
model$eval()
with_no_grad({
  # Initial Guess
  pred_0 <- model(t_X, t_mask, t_phylo, t_attn_bias)
  X_curr <- t_X * t_mask + pred_0 * (1 - t_mask)
  
  # Refinement Loop: Feed the model's own predictions back in
  for(k in 1:15) {
    # Notice we keep t_mask same (we still "don't know" the missing values)
    # But X_curr is now filled with our best guess.
    pred_k <- model(X_curr, t_mask, t_phylo, t_attn_bias)
    X_curr <- t_X * t_mask + pred_k * (1 - t_mask)
  }
  X_pred_dae <- as.matrix(X_curr$cpu())
})

# ── 7. Results ───────────────────────────────────────────────────────────────
rmse_dae <- sqrt(mean((X_truth[missing_idx] - X_pred_dae[missing_idx])^2))

cat("\n=========================================\n")
cat(" FINAL RESULTS (RMSE) \n")
cat("=========================================\n")
cat(sprintf("Rphylopars (Linear)     : %.4f\n", rmse_rphylo))
cat(sprintf("Phylo-Transformer (Deep): %.4f\n", rmse_dae))
cat("=========================================\n")

# Visualization
df_plot <- data.frame(
  Truth = c(X_truth[missing_idx], X_truth[missing_idx]),
  Prediction = c(X_pred_rphylo[missing_idx], X_pred_dae[missing_idx]),
  Method = rep(c("Rphylopars", "Phylo-Transformer"), each = length(missing_idx))
)

p <- ggplot(df_plot, aes(x = Truth, y = Prediction, color = Method)) +
  geom_point(alpha = 0.3, size=1) +
  geom_abline(color="black", linetype="dashed") +
  facet_wrap(~Method) +
  theme_minimal() +
  labs(
    title = "Imputation Challenge: AVONET",
    subtitle = paste0("Can Deep Learning beat Linear Evolution?\nTransformer RMSE: ", round(rmse_dae, 4))
  )
print(p)