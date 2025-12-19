# ==============================================================================
# AVONET CHALLENGE: Phylo-DAE (with Ecology) vs Rphylopars
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

tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

# ── 2. Data Engineering (The Secret Sauce) ───────────────────────────────────

# A) Continuous Traits (Targets)
trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

# B) Ecological Covariates (Predictors Rphylopars largely ignores)
# These explain "Why" the beak is shaped that way.
eco_cols <- c("Habitat", "Trophic.Niche", "Primary.Lifestyle")

# Filter complete cases for both Traits AND Ecology
df_main <- na.omit(avonet[, c("Species_Key", trait_cols, eco_cols)])

# Match Tree
common <- intersect(tree$tip.label, df_main$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_final <- df_main[match(tree_pruned$tip.label, df_main$Species_Key), ]

cat("Dataset Size:", nrow(df_final), "species.\n")

# C) Preprocessing
# 1. Traits: Log + Scale
X_raw <- as.matrix(df_final[, trait_cols])
X_log <- log(X_raw)
X_mean <- colMeans(X_log)
X_sd   <- apply(X_log, 2, sd)
X_truth <- scale(X_log)

# 2. Ecology: One-Hot Encoding
# Function to turn factors into binary matrix
make_one_hot <- function(df, cols) {
  out <- NULL
  for(col in cols) {
    # Convert to factor then model matrix
    f <- factor(df[[col]])
    mm <- model.matrix(~ f - 1)
    colnames(mm) <- paste0(col, "_", levels(f))
    if(is.null(out)) out <- mm else out <- cbind(out, mm)
  }
  return(out)
}

E_covariates <- make_one_hot(df_final, eco_cols)
# Standardize covariates slightly (helps convergence)
E_covariates <- scale(E_covariates, center=FALSE, scale=TRUE) 

cat("Covariates added:", ncol(E_covariates), "ecological features.\n")

# ── 3. Caching Phylogenetic Features ─────────────────────────────────────────
cache_file <- here("avonet", "phylo_cache.rds")

if (file.exists(cache_file) && nrow(readRDS(cache_file)$phylo) == nrow(df_final)) {
  cat("Loading cached tree features...\n")
  cache <- readRDS(cache_file)
  t_phylo <- torch_tensor(cache$phylo, dtype=torch_float(), device=device)
  t_adj   <- torch_tensor(cache$adj, dtype=torch_float(), device=device)
} else {
  cat("Calculating tree features (once)...\n")
  # 1. Adjacency
  D <- cophenetic(tree_pruned)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2))
  diag(A) <- 1
  rs <- rowSums(A) + 1e-8
  Dinv <- diag(1 / sqrt(rs))
  adj_mat <- Dinv %*% A %*% Dinv
  
  # 2. Spectral
  L <- diag(rowSums(A)) - A; diag(L) <- 0 # Laplacian approximation
  eig <- eigen(L, symmetric=TRUE)
  k <- 16 # Increased K to capture more detail
  phylo_mat <- eig$vectors[, (nrow(A)-k):(nrow(A)-1)]
  
  saveRDS(list(phylo=phylo_mat, adj=adj_mat), cache_file)
  t_phylo <- torch_tensor(phylo_mat, dtype=torch_float(), device=device)
  t_adj   <- torch_tensor(adj_mat, dtype=torch_float(), device=device)
}

# ── 4. Artificial Missingness ────────────────────────────────────────────────
set.seed(999) # New seed
missing_frac <- 0.25
mask_true <- matrix(1, nrow=nrow(X_truth), ncol=ncol(X_truth))
missing_idx <- sample(length(X_truth), floor(missing_frac * length(X_truth)))
mask_true[missing_idx] <- 0
X_with_NA <- X_truth
X_with_NA[missing_idx] <- NA

# ── 5. Baseline: Rphylopars ──────────────────────────────────────────────────
cat("\n--- Running Rphylopars ---\n")
# Note: We give Rphylopars ONLY traits + tree (Standard usage)
# Giving it 30+ categorical covariates makes it crash or run forever.
df_rphylo <- data.frame(species = df_final$Species_Key, X_with_NA)
res_rphylo <- phylopars(trait_data = df_rphylo, tree = tree_pruned, model = "BM", pheno_error = FALSE)
X_pred_rphylo <- res_rphylo$anc_recon[df_final$Species_Key, trait_cols]
rmse_rphylo <- sqrt(mean((X_truth[missing_idx] - X_pred_rphylo[missing_idx])^2))
cat("Rphylopars RMSE:", round(rmse_rphylo, 4), "\n")

# ── 6. The "Slayer" Model: Eco-Phylo-GNN ─────────────────────────────────────
cat("\n--- Running Eco-Phylo-GNN ---\n")

PhyloDAE <- nn_module(
  "PhyloDAE",
  initialize = function(n_traits, n_phylo, n_cov, hidden_dim = 128) {
    
    # Input: Traits + Mask + Phylo + Ecology
    input_dim <- n_traits + n_traits + n_phylo + n_cov
    
    self$enc1 <- nn_linear(input_dim, hidden_dim)
    self$bn1  <- nn_batch_norm1d(hidden_dim)
    self$enc2 <- nn_linear(hidden_dim, hidden_dim)
    
    # Graph Convolution 1
    self$msg1   <- nn_linear(hidden_dim, hidden_dim)
    self$alpha1 <- nn_parameter(torch_tensor(0.5, device=device)) # Stronger initial alpha
    
    # Graph Convolution 2 (Deeper propagation)
    self$msg2   <- nn_linear(hidden_dim, hidden_dim)
    self$alpha2 <- nn_parameter(torch_tensor(0.3, device=device))
    
    self$dec1 <- nn_linear(hidden_dim, hidden_dim)
    self$dec2 <- nn_linear(hidden_dim, n_traits)
    
    self$act  <- nn_relu()
    self$drop <- nn_dropout(0.15)
    
    self$mask_token <- nn_parameter(torch_zeros(1, n_traits, device=device))
  },
  
  forward = function(x, m, p, eco, adj) {
    # Combine all knowledge
    combined <- torch_cat(list(x, m, p, eco), dim = 2)
    
    # Encode
    h <- self$enc1(combined)
    h <- self$bn1(h)
    h <- self$act(h)
    h <- self$drop(h)
    
    h <- self$enc2(h); h <- self$act(h)
    
    # GNN Layer 1
    m_agg1 <- torch_matmul(adj, h)
    h <- h + self$alpha1 * self$msg1(m_agg1)
    
    # GNN Layer 2 (Skip connection included implicitly by +=)
    m_agg2 <- torch_matmul(adj, h)
    h <- h + self$alpha2 * self$msg2(m_agg2)
    
    # Decode
    h <- self$dec1(h); h <- self$act(h)
    out <- self$dec2(h)
    return(out)
  }
)

# Prepare Tensors
X_fill <- X_with_NA
X_fill[is.na(X_fill)] <- 0 # Initial guess is mean
t_X    <- torch_tensor(X_fill, dtype=torch_float(), device=device)
t_mask <- torch_tensor(mask_true, dtype=torch_float(), device=device)
t_eco  <- torch_tensor(E_covariates, dtype=torch_float(), device=device)

# Initialize Model
model <- PhyloDAE(n_traits=4, n_phylo=16, n_cov=ncol(E_covariates), hidden_dim=128)
model$to(device=device)
optimizer <- optim_adamw(model$parameters, lr=0.002, weight_decay=1e-4)

# Training with Harder Corruption
n_epochs <- 2000
cat("Training...\n")

for (i in 1:n_epochs) {
  model$train()
  optimizer$zero_grad()
  
  # Dynamic Corruption
  # We use a higher rate (40%) to force it to rely on Ecology/Phylo
  u <- torch_rand_like(t_X)
  is_corrupted <- (u < 0.40) & (t_mask == 1)
  
  if (as.numeric(is_corrupted$sum()$cpu()) > 0) {
    is_corrupted_f <- is_corrupted$to(dtype=torch_float())
    mask_dynamic <- t_mask * (1 - is_corrupted_f)
    
    # Replace corrupted with token
    X_in <- torch_where(is_corrupted, model$mask_token$expand_as(t_X), t_X)
    
    pred <- model(X_in, mask_dynamic, t_phylo, t_eco, t_adj)
    
    # MSE Loss only on corrupted data
    loss <- nnf_mse_loss(pred[is_corrupted], t_X[is_corrupted])
    
    loss$backward()
    optimizer$step()
    
    if(i %% 200 == 0) cat(sprintf("Epoch %d | Loss: %.4f\n", i, loss$item()))
  }
}

# ── 7. Inference with REFINE LOOP ────────────────────────────────────────────
cat("\nRunning Refinement Inference...\n")
model$eval()
with_no_grad({
  # Step 1: Initial Guess
  # We input 0s for missing, but mask=0 tells net "these are unknown"
  pred_0 <- model(t_X, t_mask, t_phylo, t_eco, t_adj)
  
  # Fill the missing spots with the first guess
  # We construct a "filled" matrix
  # obs_vals * mask + pred_vals * (1-mask)
  X_curr <- t_X * t_mask + pred_0 * (1 - t_mask)
  
  # Step 2: Refinement Loop (The "Diffusion-like" step)
  # Feed the filled matrix back in. The net now sees its own guesses.
  # We keep mask=0 so it knows which ones are still "targets" to improve.
  for(k in 1:10) {
    pred_k <- model(X_curr, t_mask, t_phylo, t_eco, t_adj)
    # Update only the missing slots
    X_curr <- t_X * t_mask + pred_k * (1 - t_mask)
  }
  
  X_pred_dae <- as.matrix(X_curr$cpu())
})

# ── 8. Final Scoreboard ──────────────────────────────────────────────────────
rmse_dae <- sqrt(mean((X_truth[missing_idx] - X_pred_dae[missing_idx])^2))

cat("\n=========================================\n")
cat(" FINAL RESULTS (RMSE) \n")
cat("=========================================\n")
cat(sprintf("Rphylopars: %.4f\n", rmse_rphylo))
cat(sprintf("Phylo-DAE : %.4f\n", rmse_dae))
cat("=========================================\n")

# Plot
df_plot <- data.frame(
  Truth = c(X_truth[missing_idx], X_truth[missing_idx]),
  Prediction = c(X_pred_rphylo[missing_idx], X_pred_dae[missing_idx]),
  Method = rep(c("Rphylopars", "Phylo-DAE"), each = length(missing_idx))
)

ggplot(df_plot, aes(x = Truth, y = Prediction, color = Method)) +
  geom_point(alpha = 0.3) +
  geom_abline(color="black", linetype="dashed") +
  facet_wrap(~Method) +
  theme_minimal() +
  labs(title = "AVONET Imputation Challenge",
       subtitle = paste0("Can we beat the baseline?\nDAE: ", round(rmse_dae, 4), " vs Rphylo: ", round(rmse_rphylo, 4)))