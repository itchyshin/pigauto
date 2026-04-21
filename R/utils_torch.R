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

# GPU memory snapshot helper for debugging large-n OOM. Returns a
# one-line summary like "23.4/42.9 GiB (alloc/reserved)" on CUDA,
# or the empty string otherwise (so it's safe to always call).
# Set PIGAUTO_DEBUG_GPU_MEM=1 to enable printing in fit_pigauto.
#
# Tries three strategies in order:
#   1. torch::cuda_memory_stats() with structured keys (standard
#      PyTorch shape: $allocated_bytes.all.current etc.)
#   2. torch::cuda_memory_stats() with flattened list fallback
#   3. nvidia-smi subprocess (works even if torch API doesn't
#      expose the bytes we need)
gpu_mem_str <- function() {
  if (!torch::cuda_is_available()) return("")

  # Strategy 1: standard PyTorch key names
  stats <- tryCatch(torch::cuda_memory_stats(), error = function(e) NULL)
  if (!is.null(stats)) {
    alloc    <- stats[["allocated_bytes.all.current"]]
    reserved <- stats[["reserved_bytes.all.current"]]
    if (!is.null(alloc) && !is.null(reserved) &&
        is.numeric(alloc) && is.numeric(reserved)) {
      return(sprintf("%.2f/%.2f GiB (alloc/reserved, via cuda_memory_stats)",
                      alloc / 1024^3, reserved / 1024^3))
    }
    # Strategy 2: pattern match any "current" bytes keys
    nm <- names(stats)
    a_key <- nm[grepl("allocated_bytes.*current", nm)][1]
    r_key <- nm[grepl("reserved_bytes.*current", nm)][1]
    if (!is.na(a_key) && !is.na(r_key)) {
      a_val <- tryCatch(stats[[a_key]], error = function(e) NULL)
      r_val <- tryCatch(stats[[r_key]], error = function(e) NULL)
      if (is.numeric(a_val) && is.numeric(r_val)) {
        return(sprintf("%.2f/%.2f GiB (alloc/reserved, matched keys: %s/%s)",
                        a_val / 1024^3, r_val / 1024^3, a_key, r_key))
      }
    }
  }

  # Strategy 3: nvidia-smi fallback
  smi <- tryCatch(
    system("nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader,nounits",
           intern = TRUE),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (!is.null(smi) && length(smi) >= 1L && nchar(smi[1])) {
    parts <- strsplit(smi[1], ",\\s*")[[1]]
    if (length(parts) >= 2L) {
      used  <- suppressWarnings(as.numeric(parts[1]))
      total <- suppressWarnings(as.numeric(parts[2]))
      if (is.finite(used) && is.finite(total)) {
        return(sprintf("%.2f/%.2f GiB (used/total, via nvidia-smi)",
                        used / 1024, total / 1024))
      }
    }
  }

  # All strategies failed -- return a marker that at least proves the
  # function got called (so we know instrumentation fires even if we
  # can't read memory numbers).
  "(memory query unavailable; CUDA alive)"
}

# Convenience logger: prints `[GPU @ tag] <mem>` when
# PIGAUTO_DEBUG_GPU_MEM=1 is set and CUDA is available. No-op
# otherwise.  Used by fit_pigauto() at phase transitions to
# localise where the ~40 GB predict-stage OOM originates.
gpu_mem_checkpoint <- function(tag) {
  if (!identical(Sys.getenv("PIGAUTO_DEBUG_GPU_MEM"), "1")) return(invisible())
  m <- gpu_mem_str()
  if (!nchar(m)) return(invisible())
  # cat to stdout so the SLURM log captures it unambiguously (message()
  # goes to stderr which some tee pipelines drop or reorder).
  cat(sprintf("[GPU @ %s] %s\n", tag, m))
  flush.console()
  invisible()
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

    if (tp %in% c("continuous", "count", "ordinal", "proportion")) {
      # MSE loss on corrupted cells for this trait's latent column
      m <- corrupt_mask[, lc]
      if (as.numeric(m$sum()$item()) == 0) next
      loss_t <- torch::nnf_mse_loss(pred[, lc][m], truth[, lc][m])
      loss_total <- loss_total + loss_t
      n_traits <- n_traits + 1L

    } else if (tp == "multi_proportion") {
      # MSE in CLR space, averaged across K latent columns.
      # Corruption treats the group as one trait: species is corrupted
      # if ANY component column is corrupted. The K CLR columns are
      # always corrupted together in preprocess (rows are complete or
      # not), so we just use lc[1] as the row mask.
      species_mask <- corrupt_mask[, lc[1]]
      n_corrupt <- as.numeric(species_mask$sum()$item())
      if (n_corrupt == 0) next
      pred_sub  <- pred[species_mask, ][, lc]
      truth_sub <- truth[species_mask, ][, lc]
      loss_t <- torch::nnf_mse_loss(pred_sub, truth_sub)
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

    } else if (tp == "zi_count") {
      # Dual loss: BCE on gate (col 1) + conditional MSE on magnitude (col 2)
      # Gate loss (BCE)
      m_gate <- corrupt_mask[, lc[1]]
      n_gate <- as.numeric(m_gate$sum()$item())
      loss_gate <- torch::torch_tensor(0.0, device = pred$device)
      if (n_gate > 0) {
        loss_gate <- torch::nnf_binary_cross_entropy_with_logits(
          pred[, lc[1]][m_gate], truth[, lc[1]][m_gate]
        )
      }
      # Magnitude loss (MSE, only for non-zero cells where truth is finite)
      m_mag <- corrupt_mask[, lc[2]]
      truth_mag <- truth[, lc[2]]
      # finite_mask: corrupted AND truth is finite (non-zero observations)
      finite_mask <- m_mag & torch::torch_isfinite(truth_mag)
      n_finite <- as.numeric(finite_mask$sum()$item())
      loss_mag <- torch::torch_tensor(0.0, device = pred$device)
      if (n_finite > 0) {
        loss_mag <- torch::nnf_mse_loss(
          pred[, lc[2]][finite_mask], truth_mag[finite_mask]
        )
      }
      # Average gate and magnitude losses (each counts as half a trait)
      if (n_gate > 0 || n_finite > 0) {
        loss_total <- loss_total + 0.5 * loss_gate + 0.5 * loss_mag
        n_traits <- n_traits + 1L
      }
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

    if (tp %in% c("continuous", "count", "ordinal", "proportion")) {
      m <- val_mask[, lc]
      if (as.numeric(m$sum()$item()) == 0) next
      l <- as.numeric(torch::nnf_mse_loss(
        pred[, lc][m], truth[, lc][m]
      )$item())
      losses <- c(losses, l)

    } else if (tp == "multi_proportion") {
      species_mask <- val_mask[, lc[1]]
      n_val <- as.numeric(species_mask$sum()$item())
      if (n_val == 0) next
      l <- as.numeric(torch::nnf_mse_loss(
        pred[species_mask, ][, lc], truth[species_mask, ][, lc]
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

    } else if (tp == "zi_count") {
      # Gate: BCE
      m_gate <- val_mask[, lc[1]]
      n_gate <- as.numeric(m_gate$sum()$item())
      l_gate <- 0
      if (n_gate > 0) {
        l_gate <- as.numeric(torch::nnf_binary_cross_entropy_with_logits(
          pred[, lc[1]][m_gate], truth[, lc[1]][m_gate]
        )$item())
      }
      # Magnitude: conditional MSE
      m_mag <- val_mask[, lc[2]]
      truth_mag <- truth[, lc[2]]
      finite_mask <- m_mag & torch::torch_isfinite(truth_mag)
      n_finite <- as.numeric(finite_mask$sum()$item())
      l_mag <- 0
      if (n_finite > 0) {
        l_mag <- as.numeric(torch::nnf_mse_loss(
          pred[, lc[2]][finite_mask], truth_mag[finite_mask]
        )$item())
      }
      if (n_gate > 0 || n_finite > 0) {
        losses <- c(losses, 0.5 * l_gate + 0.5 * l_mag)
      }
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
