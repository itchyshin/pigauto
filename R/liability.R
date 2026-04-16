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

# Posterior mean and variance of a Gaussian liability constrained to
# an interval [t_{k-1}, t_k] by the observed ordinal class k.
# thresholds: length-(K-1) vector of ordered cut-points.
# k: observed class in 1..K.
#
# @keywords internal
# @noRd
estep_liability_ordinal <- function(k, thresholds, mu_prior, sd_prior) {
  K <- length(thresholds) + 1L
  lower <- if (k == 1L) -Inf else thresholds[k - 1L]
  upper <- if (k == K)   Inf else thresholds[k]
  # Standardised truncation bounds
  a <- (lower - mu_prior) / sd_prior
  b <- (upper - mu_prior) / sd_prior
  # Truncated-normal moments (Johnson-Kotz-Balakrishnan):
  #   Z = Phi(b) - Phi(a)
  #   E = mu + sigma * (phi(a) - phi(b)) / Z
  #   Var = sigma^2 * (1 + (a*phi(a) - b*phi(b))/Z - ((phi(a)-phi(b))/Z)^2)
  phi_a <- if (is.finite(a)) stats::dnorm(a) else 0
  phi_b <- if (is.finite(b)) stats::dnorm(b) else 0
  Z     <- stats::pnorm(b) - stats::pnorm(a)
  Z     <- max(Z, 1e-300)
  m     <- (phi_a - phi_b) / Z
  ap    <- if (is.finite(a)) a * phi_a else 0
  bp    <- if (is.finite(b)) b * phi_b else 0
  var_z <- 1 + (ap - bp) / Z - m^2
  list(mean = mu_prior + sd_prior * m,
       var  = max(sd_prior^2 * var_z, 0))
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

# Soft E-step for a binary liability given proportion p in [0, 1].
# Rao-Blackwell convex combination of truncated-Gaussian posteriors.
# At p in {0, 1} reduces to estep_liability_binary.
# At p = 0.5 with symmetric prior returns the prior mean.
#
# @param p numeric in [0, 1].
# @param mu_prior numeric scalar.
# @param sd_prior positive numeric scalar.
# @return list(mean, var).
# @keywords internal
# @noRd
estep_liability_binary_soft <- function(p, mu_prior, sd_prior) {
  if (isTRUE(p == 1)) return(estep_liability_binary(y = 1L, mu_prior = mu_prior, sd_prior = sd_prior))
  if (isTRUE(p == 0)) return(estep_liability_binary(y = 0L, mu_prior = mu_prior, sd_prior = sd_prior))
  r1 <- estep_liability_binary(y = 1L, mu_prior = mu_prior, sd_prior = sd_prior)
  r0 <- estep_liability_binary(y = 0L, mu_prior = mu_prior, sd_prior = sd_prior)
  m <- p * r1$mean + (1 - p) * r0$mean
  v <- p * r1$var + (1 - p) * r0$var + p * (r1$mean - m)^2 + (1 - p) * (r0$mean - m)^2
  list(mean = m, var = max(v, 0))
}

# Plug-in E-step approximation for categorical liabilities.
# Observed-class liability gets +1 SD above the others, then project to
# sum-zero (CLR-consistent). Exact Gibbs version lives in Phase 6 EM.
#
# @keywords internal
# @noRd
estep_liability_categorical <- function(k, mu_prior, sd_prior) {
  K <- length(mu_prior)
  stopifnot(K == length(sd_prior), k >= 1L, k <= K)
  m <- mu_prior
  m[k] <- m[k] + sd_prior[k]   # boost observed class
  m <- m - mean(m)             # sum-to-zero (CLR-consistent)
  list(mean = m, var = sd_prior^2)  # prior-scale variance; refined in Phase 6
}

# Dispatcher: given a trait_map entry, the observed value (in the latent
# representation produced by preprocess_traits), and the current prior on
# this trait's liability, return the E-step posterior mean and variance.
#
# For continuous-ish types (continuous/count/proportion/multi_proportion)
# the observed value IS the liability - posterior is a point mass.
# For discrete types (binary/ordinal/categorical) the liability is latent;
# posterior follows the appropriate truncated/argmax distribution.
#
# @keywords internal
# @noRd
estep_liability <- function(tm, observed, mu_prior, sd_prior) {
  if (is.na(observed[1])) {
    # Truly missing - return prior
    return(list(mean = mu_prior, var = sd_prior^2))
  }
  tp <- tm$type
  if (tp %in% c("continuous", "count", "proportion")) {
    # Observed value IS the liability (on the preprocessed latent scale).
    list(mean = as.numeric(observed), var = rep(0, length(observed)))
  } else if (tp == "binary") {
    estep_liability_binary(y = as.integer(observed),
                           mu_prior = mu_prior, sd_prior = sd_prior)
  } else if (tp == "ordinal") {
    info <- liability_info(tm)
    # `observed` is z-scored in X_scaled: preprocess_traits encoded it as
    # (integer_class - 1 - mean) / sd. Un-z-score before recovering the
    # class index. Clamp to 1..K to guard against rounding drift at the
    # extremes (e.g. an observed value many SDs from the mean).
    K <- length(tm$levels)
    k <- as.integer(round(observed * tm$sd + tm$mean)) + 1L
    k <- max(1L, min(K, k))
    estep_liability_ordinal(k = k,
                            thresholds = info$thresholds,
                            mu_prior = mu_prior, sd_prior = sd_prior)
  } else if (tp == "categorical") {
    # observed is one-hot in the latent matrix; extract the class index
    k <- which(observed == 1)[1]
    estep_liability_categorical(k = k, mu_prior = mu_prior,
                                 sd_prior = sd_prior)
  } else if (tp %in% c("multi_proportion", "zi_count")) {
    # Phase 5 handles these; for now return prior with zero variance
    # (matches continuous semantics) - upstream callers should skip
    # these types until Phase 5 lands.
    list(mean = as.numeric(observed), var = rep(0, length(observed)))
  } else {
    stop("estep_liability: unsupported type '", tp, "'", call. = FALSE)
  }
}
