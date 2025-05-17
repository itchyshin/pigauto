######
#####3

# ════════════════════════════════════════════════════════════════════════════
# 8)  Modified Implementation for Our Dataset --------------------------------
# ════════════════════════════════════════════════════════════════════════════

# Create checkpoint directory if needed
if (!dir.exists("checkpoints")) dir.create("checkpoints", recursive = TRUE)

# Run imputation with corrected parameters
system.time(
  res <- impute_phylo(
    trait_data   = as.data.frame(traits_df_miss),
    phylo_tree   = tree,
    env_data     = env_data,
    species_id   = tree$tip.label,
    latent_dim   = 128,        # Increased capacity for complex dataset
    hidden_depth = 2,          # Deeper network for better feature extraction
    epochs       = 4000,       # Increased training duration
    lr           = 0.01,       # Adjusted learning rate
    patience     = 200,        # More tolerant early stopping
    n_samples    = 100,        # MC samples for uncertainty estimation
    ckpt_path    = "checkpoints"
  )
)

# ════════════════════════════════════════════════════════════════════════════
# 9)  Evaluation and Visualization -------------------------------------------
# ════════════════════════════════════════════════════════════════════════════

# 1. Continuous Trait Evaluation
original_cts <- as.matrix(traits_df[, grep("cnt", names(traits_df))])
imputed_cts <- as.matrix(res$completed_data[, grep("cnt", names(res$completed_data))])
missing_mask <- is.na(traits_df_miss[, grep("cnt", names(traits_df_miss))])

cat("Continuous Traits RMSE:", 
    rmse(original_cts[missing_mask], imputed_cts[missing_mask]), "\n")

# 2. Categorical Trait Evaluation
evaluate_categorical <- function(original, imputed, missing_mask) {
  matches <- original[missing_mask] == imputed[missing_mask]
  sum(matches, na.rm = TRUE) / sum(missing_mask)
}

for (col in names(traits_df)[sapply(traits_df, is.factor)]) {
  acc <- evaluate_categorical(traits_df[[col]], 
                              res$completed_data[[col]], 
                              is.na(traits_df_miss[[col]]))
  cat(sprintf("%s accuracy: %.2f%%\n", col, acc * 100))
}

# 3. Training History Visualization
plot_training_history <- function(history) {
  ggplot(data.frame(epoch = seq_along(history), loss = history)) +
    geom_line(aes(x = epoch, y = loss), color = "steelblue") +
    labs(title = "Training Loss History", 
         x = "Epoch", y = "Loss") +
    theme_minimal()
}

print(plot_training_history(res$history))

# 4. Uncertainty Analysis
plot_uncertainty <- function(uncertainty) {
  uncertainty_long <- tidyr::gather(uncertainty, "trait", "uncertainty")
  ggplot(uncertainty_long, aes(x = uncertainty)) +
    geom_histogram(fill = "steelblue", bins = 30) +
    facet_wrap(~trait, scales = "free") +
    labs(title = "Imputation Uncertainty Distribution",
         x = "Uncertainty", y = "Count") +
    theme_minimal()
}

print(plot_uncertainty(res$uncertainty))

# ════════════════════════════════════════════════════════════════════════════
# 10) Final Output ----------------------------------------------------------
# ════════════════════════════════════════════════════════════════════════════

# Save results
write.csv(res$completed_data, "imputed_data.csv", row.names = FALSE)
saveRDS(res, "imputation_results.rds")

# Print final message
cat("\nImputation completed successfully!\n",
    "Results saved to:\n",
    "- imputed_data.csv\n", 
    "- imputation_results.rds\n",
    "Training stopped at epoch", res$stopped_epoch, "\n")


# grid search

library(purrr)

grid <- expand.grid(
  latent_dim   = c(16, 32, 64),
  hidden_mult  = c(1, 2),
  dropout_rate = c(0.0, 0.3),
  lr           = c(1e-3, 1e-2),
  weight_decay = c(0, 1e-4),
  epochs       = c(500, 1000),
  stringsAsFactors = FALSE
)

results <- pmap_dfr(grid, function(latent_dim, hidden_mult,
                                   dropout_rate, lr,
                                   weight_decay, epochs) {
  # retrain with these hyperparams:
  fit <- impute_phylo(
    trait_data            = as.data.frame(obs_train),
    phylo_tree            = tree,
    env_data              = as.data.frame(env_data),
    species_id            = tree$tip.label,
    latent_dim            = latent_dim,
    epochs                = epochs,
    n_samples             = 100,            # keep fixed for tuning
    dropout_rate          = dropout_rate,
    hidden_mult           = hidden_mult,
    lr                    = lr,
    weight_decay          = weight_decay,
    mask_obs_uncertainty  = TRUE,
    restore_observed      = TRUE,
    validation_mask       = val_mask       # see note below
  )
  
  # compute RMSE on the *validation* positions:
  pred  <- as.matrix(fit$completed_data)[val_mask]
  truth <- obs_traits[val_mask]
  rmse  <- rmse(truth, pred)
  
  tibble(
    latent_dim, hidden_mult, dropout_rate, lr,
    weight_decay, epochs, val_rmse = rmse
  )
})

best <- results[which.min(results$val_rmse), ]
print(best)

#####