# R/calibration_three_way.R
# Safety-floor (three-way mean-gate) calibration helpers.
# See specs/2026-04-23-safety-floor-mean-gate-design.md.

# Build the 3-simplex grid at resolution `step`.
#
# Returns an integer-divided grid of (r_bm, r_gnn, r_mean) points where
# each row sums to exactly 1 and each entry is a non-negative multiple
# of `step`. At step = 0.05 this is 231 rows; step = 0.25 is 15 rows.
#
# The (1, 0, 0), (0, 1, 0), and (0, 0, 1) corners are always included:
# (0, 0, 1) is the grand-mean baseline candidate that gives the safety
# invariant `calibrated_loss <= mean_loss` by construction.
#
# @param step numeric in (0, 1]. Default 0.05 (231 points).
# @return numeric matrix with columns c("r_bm", "r_gnn", "r_mean").
#
# @noRd
simplex_grid <- function(step = 0.05) {
  stopifnot(is.numeric(step), length(step) == 1L, step > 0, step <= 1)
  k   <- as.integer(round(1 / step))
  if (abs(k * step - 1) > 1e-8) {
    stop("step must evenly divide 1; got step = ", step)
  }
  n_rows <- (k + 1L) * (k + 2L) / 2L
  out <- matrix(NA_real_, nrow = n_rows, ncol = 3L)
  row <- 1L
  for (i in 0:k) {
    for (j in 0:(k - i)) {
      out[row, ] <- c(i, j, k - i - j) / k
      row <- row + 1L
    }
  }
  colnames(out) <- c("r_bm", "r_gnn", "r_mean")
  out
}

# Compute the safety-floor mean-baseline scalar for one latent column.
#
# Returns a single numeric on the same latent scale as the column's
# bm_val / gnn_val (z-score for continuous/count/ordinal/proportion,
# logit for binary/categorical, z-score of magnitude for zi_count).
#
# For continuous-family columns: plain training grand mean on the
# already-z-scored `x_col`. This is approximately 0 (not exactly 0
# because z-scoring uses ALL preprocessed observations including ones
# that later become held-out test / val).
#
# For binary / categorical sub-columns: class-1 frequency clipped to
# [0.01, 0.99], then `qlogis()` to move to logit scale.
#
# @param x_col       numeric, one latent column, length n_obs
# @param train_mask  logical, TRUE for training-observed rows
# @param trait_type  character, one of "continuous", "count", "ordinal",
#                    "proportion", "zi_mag", "binary", "categorical"
# @return scalar numeric
#
# @noRd
mean_baseline_scalar <- function(x_col, train_mask, trait_type) {
  stopifnot(length(x_col) == length(train_mask),
            is.logical(train_mask),
            length(trait_type) == 1L)
  keep <- train_mask & !is.na(x_col)
  if (sum(keep) == 0L) return(0)
  if (trait_type %in% c("continuous", "count", "ordinal",
                          "proportion", "zi_mag")) {
    mean(x_col[keep])
  } else if (trait_type %in% c("binary", "categorical")) {
    p <- mean(x_col[keep])
    p <- pmin(pmax(p, 0.01), 0.99)
    stats::qlogis(p)
  } else {
    stop("mean_baseline_scalar: unsupported trait_type = ", trait_type)
  }
}
