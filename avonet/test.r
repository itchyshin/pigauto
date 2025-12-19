# ==============================================================================
# AVONET BENCHMARK: Pure Morphology (Can NN beat Linear Models?)
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(dplyr)
library(here)

# ── 1. Device & Data ──────────────────────────────────────────────────────────
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

# Load Data
tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))

# Clean & Match
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

# ── 2. Select ONLY Continuous Traits (No Ecology) ─────────────────────────────
# We use these to predict each other (e.g., Wing -> Mass)
trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

# Filter complete cases to establish Ground Truth
df_traits <- na.omit(avonet[, c("Species_Key", trait_cols)])

# Match Tree
common <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

cat("Dataset:", nrow(df_final), "species x", length(trait_cols), "traits\n")

# Preprocessing: Log-Transform (crucial for size data) & Scale
X_raw <- as.matrix(df_final[, trait_cols])
X_log <- log(X_raw)
X_mean <- colMeans(X_log)
X_sd   <- apply(X_log, 2, sd)
X_truth <- scale(X_log)

# ── 3. Tree Feature Engineering (Cached) ─────────────────────────────────────
# We use two views of the tree:
# 1. Spectral Embeddings (Global coordinates)
# 2. Adjacency Matrix (Local neighbor smoothing)

cache_file <- here("avonet", "phylo_cache_morpho.rds")

if (file.exists(cache_file)) {
  cat("Loading cached tree features...\n")
  cache <- readRDS(cache_file)
  # Simple check to ensure cache matches current data size
  if(nrow(cache$phylo) == nrow(df_final)) {
    t_phylo <- torch_tensor(cache$phylo, dtype=torch_float(), device=device)
    t_adj   <- torch_tensor(cache$adj, dtype=torch_float(), device=device) 
  } else {
    stop("Cache size mismatch! Delete .rds file and rerun.")
  }
} else {
  cat("Computing tree features (this takes time)... \n")
  
  # A) Adjacency (Gaussian Kernel on Cophenetic Distance)
  D <- cophenetic(tree_pruned)
  sigma <- median(D) * 0.35 # Sharper kernel for more local info
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 1 # Self-loops
  # Row-Normalize
  rs <- rowSums(A) + 1e-8
  Dinv <- diag(1 / sqrt(rs))
  adj_mat <- Dinv %*% A %*% Dinv
  
  # B) Spectral Features (Laplacian Eigenmaps)
  # We use the same A but zero-diagonal for Laplacian
  A_spec <- A; diag(A_spec) <- 0
  L <- diag(rowSums(A_spec)) - A_spec
  eig <- eigen(L, symmetric = TRUE)
  # Top K smallest non-zero eigenvectors
  k <- 16 
  phylo_mat <- eig$vectors[, (nrow(A)-k):(nrow(A)-1)]
  
  saveRDS(list(phylo=phylo_mat, adj=adj_mat), cache_file)
  
  t_phylo <- torch_tensor(phylo_mat, dtype=torch_float(), device=device)
  t_adj   <- torch_tensor(adj_mat, dtype=torch_float(), device=device)
}

# ── 4. Artificial Missingness ────────────────────────────────────────────────
set.seed(555)
missing_frac <- 0.25

# Mask: 0 = Missing, 1 = Observed
mask_true <- matrix(1, nrow=nrow(X_truth), ncol=ncol(X_truth))
missing_idx <- sample(length(X_truth), floor(missing_frac * length(X_truth)))
mask_true[missing_idx] <- 0
X_with_NA <- X_truth
X_with_NA[missing_idx] <- NA

# ── 5. Baseline: Rphylopars ──────────────────────────────────────────────────
cat("\n--- Running Rphylopars ---\n")
# "BM" is linear. "OU" is better but linear.
df_rphylo <- data.frame(species = df_final$Species_Key, X_with_NA)
start_t <- Sys.time()
res_rphylo <- phylopars(trait_data = df_rphylo, tree = tree_pruned, model = "BM", pheno_error = FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]
end_t <- Sys.time()

rmse_rphylo <- sqrt(mean((X_truth[missing_idx] - X_pred_rphylo[missing_idx])^2))
cat("Rphylopars RMSE:", round(rmse_rphylo, 4), "| Time:", round(difftime(end_t, start_t, units="secs"),2), "s\n")

# ── 6. Phylo-DAE (Optimized for Pure Morphology) ─────────────────────────────
cat("\n--- Running Phylo-DAE ---\n")

# Improved Architecture: Residual GNN
PhyloDAE <- nn_module(
  "PhyloDAE",
  initialize = function(n_traits, n_phylo, hidden_dim = 128) {
    
    # Input: Traits + Mask + Phylo
    input_dim <- n_traits + n_traits + n_phylo 
    
    self$enc1 <- nn_linear(input_dim, hidden_dim)
    self$bn1  <- nn_batch_norm1d(hidden_dim)
    self$enc2 <- nn_linear(hidden_dim, hidden_dim)
    
    # Graph Layer (with learnable gating)
    self$msg   <- nn_linear(hidden_dim, hidden_dim)
    self$gate  <- nn_linear(hidden_dim, hidden_dim) 
    
    self$dec1 <- nn_linear(hidden_dim, hidden_dim)
    self$dec2 <- nn_linear(hidden_dim, n_traits)
    
    self$act  <- nn_gelu() # GeLU is often better than ReLU for regression
    self$drop <- nn_dropout(0.1)
    
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits, device=device))
  },
  
  forward = function(x, m, p, adj) {
    # x: Data, m: Mask, p: Phylo
    combined <- torch_cat(list(x, m, p), dim = 2)
    
    # Encoder
    h <- self$enc1(combined)
    h <- self$bn1(h)
    h <- self$act(h)
    h <- self$drop(h)
    
    h_resid <- h # Save for skip connection
    h <- self$enc2(h); h <- self$act(h)
    
    # Graph Message Passing
    m_agg <- torch_matmul(adj, h) # Aggregate neighbors
    
    # Gated Residual Update: H = H + sigmoid(Gate(H)) * Message(Agg)
    gate_val <- torch_sigmoid(self$gate(h))
    h <- h_resid + gate_val * self$msg(m_agg)
    
    # Decoder
    h <- self$dec1(h); h <- self$act(h)
    out <- self$dec2(h)
    
    return(out)
  }
)

# Prepare Tensors
X_fill <- X_with_NA
X_fill[is.na(X_fill)] <- 0
t_X    <- torch_tensor(X_fill, dtype=torch_float(), device=device)
t_mask <- torch_tensor(mask_true, dtype=torch_float(), device=device)

model <- PhyloDAE(n_traits=4, n_phylo=16, hidden_dim=144) # Increased width
model$to(device=device)
optimizer <- optim_adamw(model$parameters, lr=0.002, weight_decay=1e-5)

# Training Loop
n_epochs <- 2500
cat("Training with Masked Autoencoding...\n")

for (i in 1:n_epochs) {
  model$train()
  optimizer$zero_grad()
  
  # Dynamic Corruption Strategy
  # We randomly hide 30% of the *observed* data
  u <- torch_rand_like(t_X)
  is_corrupted <- (u < 0.30) & (t_mask == 1)
  
  if (as.numeric(is_corrupted$sum()$cpu()) > 0) {
    is_corrupted_f <- is_corrupted$to(dtype=torch_float())
    
    # The network sees: 
    # 1. Corrupted Values -> Replaced by Mask Token
    X_in <- torch_where(is_corrupted, model$mask_token$expand_as(t_X), t_X)
    
    # 2. Mask -> Tells net "This is missing" (0)
    # Important: We tell the net the corrupted spots are missing (0)
    # so it learns to fill them.
    mask_dynamic <- t_mask * (1 - is_corrupted_f)
    
    pred <- model(X_in, mask_dynamic, t_phylo, t_adj)
    
    # Loss only on the corrupted spots (Self-Supervised)
    loss <- nnf_mse_loss(pred[is_corrupted], t_X[is_corrupted])
    
    loss$backward()
    optimizer$step()
    
    if (i %% 500 == 0) cat(sprintf("Epoch %d | Loss: %.5f\n", i, loss$item()))
  }
}

# ── 7. Inference with Iterative Refinement ───────────────────────────────────
cat("\nRefining Predictions...\n")
model$eval()
with_no_grad({
  # Pass 1: Initial Guess
  # Input: Known data + 0 for missing + Mask indicating which are missing
  pred_0 <- model(t_X, t_mask, t_phylo, t_adj)
  
  # Fill the gaps
  X_curr <- t_X * t_mask + pred_0 * (1 - t_mask)
  
  # Refinement Loop (10 steps)
  # We feed the *filled* matrix back in.
  # Crucially, we keep the mask same (0 for targets), so the net knows 
  # "I am still guessing these, but now I have context from my previous guess".
  for(k in 1:10) {
    pred_k <- model(X_curr, t_mask, t_phylo, t_adj)
    X_curr <- t_X * t_mask + pred_k * (1 - t_mask)
  }
  
  X_pred_dae <- as.matrix(X_curr$cpu())
})

# ── 8. Results ───────────────────────────────────────────────────────────────
rmse_dae <- sqrt(mean((X_truth[missing_idx] - X_pred_dae[missing_idx])^2))

cat("\n=========================================\n")
cat(" RESULTS (RMSE on Standardized Log Data)\n")
cat("=========================================\n")
cat(sprintf("Rphylopars: %.4f\n", rmse_rphylo))
cat(sprintf("Phylo-DAE : %.4f\n", rmse_dae))
cat("=========================================\n")
if (rmse_dae < rmse_rphylo) cat(">> Phylo-DAE WINS! <<\n") else cat(">> Rphylopars Wins. <<\n")

# Visualization
df_plot <- data.frame(
  Truth = c(X_truth[missing_idx], X_truth[missing_idx]),
  Prediction = c(X_pred_rphylo[missing_idx], X_pred_dae[missing_idx]),
  Method = rep(c("Rphylopars", "Phylo-DAE"), each = length(missing_idx))
)

p <- ggplot(df_plot, aes(x = Truth, y = Prediction, color = Method)) +
  geom_point(alpha = 0.3, size=1) +
  geom_abline(color="black", linetype="dashed") +
  facet_wrap(~Method) +
  theme_minimal() +
  labs(
    title = "Imputation: Neural Network vs Linear Model",
    subtitle = paste0("RMSE: DAE=", round(rmse_dae,3), " vs Rphylo=", round(rmse_rphylo,3))
  )
print(p)