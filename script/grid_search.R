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