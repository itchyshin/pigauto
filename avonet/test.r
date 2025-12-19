# ==============================================================================
# AVONET & PHYLO-DAE: Real Data Benchmark
# ==============================================================================
# 1. Loads AVONET (Traits) and Hackett Tree (Phylogeny).
# 2. Selects continuous traits: Mass, Beak, Tarsus, Wing.
# 3. Artificially deletes data (25% missingness).
# 4. Compares imputation accuracy (RMSE) of Rphylopars vs Phylo-DAE.
# ==============================================================================

library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

# ── 1. Device Setup ──────────────────────────────────────────────────────────
device <- if (cuda_is_available()) {
  torch_device("cuda")
} else if (backends_mps_is_available()) {
  torch_device("mps") # Mac Metal
} else {
  torch_device("cpu")
}
cat("Using device:", as.character(device), "\n")

# ── 2. Data Loading & Preprocessing ──────────────────────────────────────────
cat("\n--- Loading Data ---\n")

# Load Tree
# Ensure the file is in an 'avonet' subfolder relative to your project root
tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))

# Load AVONET Data
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))

# Clean Species Names: AVONET uses "Genus species", Tree uses "Genus_species"
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

# Select Specific Continuous Traits
# Note: Using Log transformation as biological size data is usually log-normal
trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
cat("Selected Traits:", paste(trait_cols, collapse=", "), "\n")

# Filter Data: Select only relevant columns and Species Key
df_traits <- avonet[, c("Species_Key", trait_cols)]

# 1. Filter for complete cases FIRST (to have a Ground Truth for evaluation)
df_traits <- na.omit(df_traits)

# 2. Match Tree and Data
common_species <- intersect(tree$tip.label, df_traits$Species_Key)
tree_pruned <- keep.tip(tree, common_species)
df_final <- df_traits[match(tree_pruned$tip.label, df_traits$Species_Key), ]

cat("Species with complete data matching tree:", nrow(df_final), "\n")

# Prepare Ground Truth Matrix (Log Transformed & Scaled)
X_raw <- as.matrix(df_final[, trait_cols])
X_log <- log(X_raw) 
# Z-score scaling (crucial for Neural Networks)
X_mean <- colMeans(X_log)
X_sd   <- apply(X_log, 2, sd)
X_truth <- scale(X_log)

# ── 3. Introduce Artificial Missingness ──────────────────────────────────────
set.seed(123)
missing_frac <- 0.25 # 25% missing data
n_total <- length(X_truth)
n_missing <- floor(missing_frac * n_total)

# Create Mask (1 = Observed, 0 = Missing)
mask <- matrix(1, nrow=nrow(X_truth), ncol=ncol(X_truth))
missing_indices <- sample(n_total, n_missing)
mask[missing_indices] <- 0

# Create Input Data with NAs
X_with_NA <- X_truth
X_with_NA[missing_indices] <- NA

cat("Introduced", n_missing, "missing entries (", missing_frac*100, "%).\n")

# ── 4. Method 1: Rphylopars (Baseline) ───────────────────────────────────────
cat("\n--- Running Rphylopars (Baseline) ---\n")
start_time <- Sys.time()

# Prepare data frame for Rphylopars
# FIX: 'species' must be the FIRST column and must match tree tip labels
df_rphylo <- data.frame(species = df_final$Species_Key, X_with_NA)

# Run Imputation (Using Brownian Motion model for speed on large trees)
res_rphylo <- phylopars(trait_data = df_rphylo, tree = tree_pruned, model = "BM")

# Extract predicted values
# FIX: Match rows by species name to ensure alignment
trait_names <- colnames(X_with_NA)
X_pred_rphylo <- as.matrix(res_rphylo$anc_recon[df_final$Species_Key, trait_names])

end_time <- Sys.time()
cat("Rphylopars finished in:", round(difftime(end_time, start_time, units="secs"), 2), "s\n")

# ── 5. Method 2: Phylo-DAE (Neural Network) ──────────────────────────────────
cat("\n--- Running Phylo-DAE ---\n")

# A) Generate Spectral Features from Tree (The "Phylo" part)
get_spectral_features <- function(tree, k = 10) {
  # Fast calculation using sparse matrices for large trees
  n <- length(tree$tip.label)
  # Inverse phylogenetic covariance approximation (Graph Laplacian of the tree)
  # Using simple cophenetic distance for this demo (approximation)
  D <- cophenetic(tree)
  sigma <- median(D) * 0.5
  A <- exp(- (D^2) / (2 * sigma^2)); diag(A) <- 0
  L <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric = TRUE)
  # Return k smallest non-zero eigenvectors
  coords <- eig$vectors[, (n - k):(n - 1)] 
  return(as.matrix(coords))
}

cat("Generating spectral embeddings (this may take a moment for large trees)...\n")
# Using k=8 spectral dimensions
phylo_emb <- get_spectral_features(tree_pruned, k = 8)
phylo_emb <- scale(phylo_emb)

# B) Define the Denoising Autoencoder Module
PhyloDAE <- nn_module(
  "PhyloDAE",
  initialize = function(n_traits, n_phylo, n_hidden = 64) {
    # Encoder takes: [Traits (filled with 0), Mask (0/1), Phylo_Embeddings]
    input_dim <- n_traits + n_traits + n_phylo
    
    self$encoder <- nn_sequential(
      nn_linear(input_dim, n_hidden),
      nn_relu(),
      nn_dropout(0.2),
      nn_linear(n_hidden, n_hidden),
      nn_relu()
    )
    
    # Decoders for Mean and Log-Variance
    self$head_mu <- nn_linear(n_hidden, n_traits)
    self$head_var <- nn_linear(n_hidden, n_traits)
  },
  
  forward = function(x, m, p) {
    # Replace NAs in x with 0 for the network input
    x0 <- x$clone()
    x0[is.nan(x0)] <- 0 
    
    # Concatenate [Data, Mask, Phylo]
    combined <- torch_cat(list(x0, m, p), dim = 2)
    
    z <- self$encoder(combined)
    mu <- self$head_mu(z)
    logvar <- self$head_var(z)
    list(mu = mu, logvar = logvar)
  }
)

# C) Prepare Tensors
# Note: Input traits NA should be 0, but we also pass the mask so the net knows it's missing
X_in <- X_with_NA
X_in[is.na(X_in)] <- 0

t_traits <- torch_tensor(X_in, dtype = torch_float(), device = device)
t_mask   <- torch_tensor(mask, dtype = torch_float(), device = device)
t_phylo  <- torch_tensor(phylo_emb, dtype = torch_float(), device = device)
t_truth  <- torch_tensor(X_truth, dtype = torch_float(), device = device)

# D) Train Loop
model <- PhyloDAE(n_traits = ncol(X_truth), n_phylo = ncol(phylo_emb))
model$to(device = device)
optimizer <- optim_adam(model$parameters, lr = 0.005)

cat("Training Phylo-DAE...\n")
n_epochs <- 1000
loss_hist <- numeric(n_epochs)

for (i in 1:n_epochs) {
  optimizer$zero_grad()
  
  out <- model(t_traits, t_mask, t_phylo)
  
  # Heteroscedastic Loss (NLL)
  # Only calculate loss on OBSERVED data (mask == 1)
  # L = 0.5 * exp(-logvar) * (y - mu)^2 + 0.5 * logvar
  recon_loss <- 0.5 * torch_exp(-out$logvar) * (t_truth - out$mu)^2 + 0.5 * out$logvar
  loss <- (recon_loss * t_mask)$sum() / t_mask$sum()
  
  loss$backward()
  optimizer$step()
  
  loss_hist[i] <- as.numeric(loss$item())
  if (i %% 200 == 0) cat(sprintf("Epoch %d | Loss: %.4f\n", i, loss_hist[i]))
}

# E) Predict
model$eval()
with_no_grad({
  out_final <- model(t_traits, t_mask, t_phylo)
  X_pred_dae <- as.matrix(out_final$mu$cpu())
})

# ── 6. Evaluation & Comparison ───────────────────────────────────────────────

# Extract missing values for evaluation
truth_vals <- X_truth[missing_indices]
pred_rphylo_vals <- X_pred_rphylo[missing_indices]
pred_dae_vals    <- X_pred_dae[missing_indices]

# Calculate RMSE
rmse_rphylo <- sqrt(mean((truth_vals - pred_rphylo_vals)^2))
rmse_dae    <- sqrt(mean((truth_vals - pred_dae_vals)^2))

cat("\n=========================================\n")
cat(" RESULTS SUMMARY (RMSE on Missing Data) \n")
cat("=========================================\n")
cat("Dataset Size:   ", nrow(X_truth), "species\n")
cat("Traits:         ", colnames(X_truth), "\n")
cat("Missing %:      ", missing_frac * 100, "%\n\n")
cat(sprintf("Rphylopars RMSE: %.4f\n", rmse_rphylo))
cat(sprintf("Phylo-DAE  RMSE: %.4f\n", rmse_dae))
cat("=========================================\n")

# ── 7. Visualization ─────────────────────────────────────────────────────────

# Create plotting data frame
df_plot <- data.frame(
  Truth = c(truth_vals, truth_vals),
  Prediction = c(pred_rphylo_vals, pred_dae_vals),
  Method = rep(c("Rphylopars", "Phylo-DAE"), each = length(truth_vals))
)

p <- ggplot(df_plot, aes(x = Truth, y = Prediction, color = Method)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  facet_wrap(~Method) +
  theme_minimal() +
  scale_color_manual(values = c("Rphylopars" = "#E69F00", "Phylo-DAE" = "#56B4E9")) +
  labs(
    title = paste("Imputation Performance on AVONET"),
    subtitle = paste("RMSE: Rphylo =", round(rmse_rphylo,3), "| Phylo-DAE =", round(rmse_dae,3)),
    x = "True Value (Scaled Log-Transformed)",
    y = "Imputed Value"
  )

print(p)