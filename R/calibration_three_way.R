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
