# Internal torch helpers -- not exported

# Detect best available device (CUDA > MPS > CPU)
get_device <- function() {
  if (torch::cuda_is_available()) {
    torch::torch_device("cuda")
  } else if (torch::backends_mps_is_available()) {
    torch::torch_device("mps")
  } else {
    torch::torch_device("cpu")
  }
}

# RMSE computed entirely in torch (avoids per-epoch CPU copy)
rmse_torch <- function(pred, truth, mask_bool) {
  d <- pred - truth
  torch::torch_sqrt(torch::torch_mean(d[mask_bool]$pow(2)))
}

# Compute mixed-type loss on corrupted cells
# pred, truth: (n x p_latent) tensors
# corrupt_mask: (n x p_latent) bool tensor (TRUE = corrupted)
# trait_map: list of trait descriptors
compute_mixed_loss <- function(pred, truth, corrupt_mask, trait_map) {
  loss_total <- torch::torch_tensor(0.0, device = pred$device,
                                     requires_grad = TRUE)
  n_traits <- 0L

  for (tm in trait_map) {
    lc <- tm$latent_cols
    tp <- tm$type

    if (tp %in% c("continuous", "count", "ordinal")) {
      # MSE loss on corrupted cells for this trait's latent column
      m <- corrupt_mask[, lc]
      if (as.numeric(m$sum()$item()) == 0) next
      loss_t <- torch::nnf_mse_loss(pred[, lc][m], truth[, lc][m])
      loss_total <- loss_total + loss_t
      n_traits <- n_traits + 1L

    } else if (tp == "binary") {
      # BCE with logits on corrupted cells
      m <- corrupt_mask[, lc]
      if (as.numeric(m$sum()$item()) == 0) next
      loss_t <- torch::nnf_binary_cross_entropy_with_logits(
        pred[, lc][m], truth[, lc][m]
      )
      loss_total <- loss_total + loss_t
      n_traits <- n_traits + 1L

    } else if (tp == "categorical") {
      # Cross-entropy on corrupted species for this trait
      K <- tm$n_latent
      # Corruption mask for categorical: species is corrupted if ANY of K cols
      species_mask <- corrupt_mask[, lc[1]]
      n_corrupt <- as.numeric(species_mask$sum()$item())
      if (n_corrupt == 0) next

      # logits: (n_corrupt x K)
      logits_corrupt <- pred[species_mask, ][, lc]
      # truth one-hot: (n_corrupt x K) -> class indices
      truth_onehot <- truth[species_mask, ][, lc]
      targets <- truth_onehot$argmax(dim = 2L)

      loss_t <- torch::nnf_cross_entropy(logits_corrupt, targets)
      loss_total <- loss_total + loss_t
      n_traits <- n_traits + 1L
    }
  }

  if (n_traits > 0L) loss_total / n_traits else loss_total
}

# Compute composite validation metric (per-trait, averaged)
# Returns a scalar loss (lower = better)
composite_val_loss <- function(pred, truth, val_mask, trait_map) {
  losses <- numeric(0)

  for (tm in trait_map) {
    lc <- tm$latent_cols
    tp <- tm$type

    if (tp %in% c("continuous", "count", "ordinal")) {
      m <- val_mask[, lc]
      if (as.numeric(m$sum()$item()) == 0) next
      l <- as.numeric(torch::nnf_mse_loss(
        pred[, lc][m], truth[, lc][m]
      )$item())
      losses <- c(losses, l)

    } else if (tp == "binary") {
      m <- val_mask[, lc]
      if (as.numeric(m$sum()$item()) == 0) next
      l <- as.numeric(torch::nnf_binary_cross_entropy_with_logits(
        pred[, lc][m], truth[, lc][m]
      )$item())
      losses <- c(losses, l)

    } else if (tp == "categorical") {
      species_mask <- val_mask[, lc[1]]
      n_val <- as.numeric(species_mask$sum()$item())
      if (n_val == 0) next
      logits_val <- pred[species_mask, ][, lc]
      truth_oh   <- truth[species_mask, ][, lc]
      targets    <- truth_oh$argmax(dim = 2L)
      l <- as.numeric(torch::nnf_cross_entropy(logits_val, targets)$item())
      losses <- c(losses, l)
    }
  }

  if (length(losses) == 0) return(Inf)
  mean(losses)
}


# Cosine annealing LR with linear warmup.
# Returns the scheduled learning rate for the given epoch.
cosine_lr <- function(epoch, warmup, total, lr_max, lr_min = 1e-5) {
  if (warmup > 0 && epoch <= warmup) {
    return(lr_max * epoch / warmup)
  }
  if (total <= warmup) return(lr_max)
  progress <- (epoch - warmup) / (total - warmup)
  lr_min + 0.5 * (lr_max - lr_min) * (1 + cos(pi * progress))
}
