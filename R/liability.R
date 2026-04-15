#' Liability schema for a trait-map entry
#'
#' Returns the per-type liability contract used by the Level-C joint
#' baseline. Every trait, regardless of observed type, maps to one or
#' more continuous underlying liabilities that evolve jointly by BM on
#' the tree. Discrete types add a thresholding / argmax constraint on
#' top of the continuous liability.
#'
#' @param tm A single entry of `pigauto_data$trait_map`.
#' @return A list with components:
#'   \describe{
#'     \item{n_liability}{Integer, number of continuous liability
#'       dimensions this trait occupies in the joint vector.}
#'     \item{kind}{Character. `"continuous"` (observed = liability),
#'       `"threshold"` (binary), `"ordered_threshold"` (ordinal),
#'       `"argmax"` (categorical, K liabilities),
#'       `"sum_zero"` (multi_proportion, K CLR liabilities),
#'       `"mixed_2"` (zi_count: gate + magnitude).}
#'     \item{thresholds}{Numeric vector or `NULL`. Cut-points on the
#'       liability scale for threshold-style discrete observations.}
#'   }
#'
#' @keywords internal
#' @noRd
liability_info <- function(tm) {
  tp <- tm$type
  if (tp %in% c("continuous", "count", "proportion")) {
    list(n_liability = 1L, kind = "continuous", thresholds = NULL)
  } else if (tp == "binary") {
    list(n_liability = 1L, kind = "threshold", thresholds = 0)
  } else if (tp == "ordinal") {
    K <- length(tm$levels)
    list(n_liability = 1L, kind = "ordered_threshold",
         thresholds = seq_len(K - 1L) - 0.5)
  } else if (tp == "categorical") {
    list(n_liability = tm$n_latent, kind = "argmax", thresholds = NULL)
  } else if (tp == "zi_count") {
    list(n_liability = 2L, kind = "mixed_2", thresholds = list(gate = 0))
  } else if (tp == "multi_proportion") {
    list(n_liability = tm$n_latent, kind = "sum_zero", thresholds = NULL)
  } else {
    stop("No liability contract defined for type '", tp, "'", call. = FALSE)
  }
}

# Posterior mean and variance of a Gaussian liability truncated by a
# binary observation. Uses standard truncated-normal formulas.
# y = 1 means liability > 0; y = 0 means liability <= 0.
#
# @keywords internal
# @noRd
estep_liability_binary <- function(y, mu_prior, sd_prior) {
  # Standardise the truncation point
  alpha <- (0 - mu_prior) / sd_prior
  z     <- stats::dnorm(alpha) / stats::pnorm(alpha, lower.tail = (y == 0))
  sign  <- if (y == 1) 1 else -1
  # Truncated-normal mean (Johnson-Kotz-Balakrishnan):
  #   E[X | X > 0] = mu + sigma * phi(alpha) / (1 - Phi(alpha))
  #   E[X | X < 0] = mu - sigma * phi(alpha) / Phi(alpha)
  mean_post <- mu_prior + sign * sd_prior * z
  # Truncated-normal variance:
  #   Var = sigma^2 * (1 + alpha*z - z^2)   for y = 1 with alpha = (0-mu)/sigma
  # (alpha is negated for y = 0; sign above handles location).
  var_post <- sd_prior^2 * (1 + sign * alpha * z - z^2)
  var_post <- max(var_post, 0)  # numerical guard
  list(mean = mean_post, var = var_post)
}
