# Internal: non-BM trait simulation for benchmarking
#
# Generates species trait data under models that deviate from Brownian Motion,
# providing scenarios where the GNN correction can add value over the BM
# baseline.
#
# Not exported -- used by script/benchmark_simulation.R only.

#' Simulate non-BM trait data for benchmarking
#'
#' Generates species trait data under models that deviate from Brownian Motion.
#' Supports OU (stabilising selection), regime shifts (clade-specific optima),
#' and non-linear inter-trait correlations.
#'
#' @param tree Object of class \code{"phylo"}.
#' @param n_traits Integer. Number of traits to simulate (default 4).
#' @param scenario Character. One of \code{"OU"}, \code{"regime_shift"},
#'   \code{"nonlinear"}.
#' @param alpha Numeric. OU pull strength (only for \code{"OU"}).
#' @param theta Numeric. OU optimum (only for \code{"OU"}).
#' @param sigma Numeric. BM diffusion rate.
#' @param shift_magnitude Numeric. Regime shift size in SD units
#'   (only for \code{"regime_shift"}).
#' @param seed Integer seed or \code{NULL}.
#' @return A data.frame with species as rownames.
#' @importFrom ape rtree rTraitCont
#' @export
simulate_non_bm <- function(tree, n_traits = 4,
                            scenario = c("OU", "regime_shift", "nonlinear"),
                            alpha = 2.0, theta = 0, sigma = 1.0,
                            shift_magnitude = 2.0, seed = NULL) {

  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)

  n <- length(tree$tip.label)
  sp <- tree$tip.label

  traits <- switch(scenario,
    OU             = simulate_ou(tree, n_traits, alpha, theta, sigma),
    regime_shift   = simulate_regime_shift(tree, n_traits, sigma, shift_magnitude),
    nonlinear      = simulate_nonlinear(tree, n_traits, sigma)
  )

  rownames(traits) <- sp
  traits
}


# ---- OU: Ornstein-Uhlenbeck ------------------------------------------------
#
# Traits evolve under stabilising selection toward an optimum (theta).
# Higher alpha = stronger pull, more deviation from BM expectations.
# The BM baseline assumes no pull, so it over-imputes variation.

simulate_ou <- function(tree, n_traits, alpha, theta, sigma) {
  df <- data.frame(row.names = tree$tip.label)
  for (j in seq_len(n_traits)) {
    # ape::rTraitCont supports OU natively
    trait <- ape::rTraitCont(
      tree,
      model     = "OU",
      sigma     = sigma,
      alpha     = alpha,
      theta     = theta + (j - 1) * 0.5,  # slight offset per trait
      root.value = theta + (j - 1) * 0.5
    )
    df[[paste0("trait", j)]] <- as.numeric(trait[tree$tip.label])
  }
  df
}


# ---- Regime shift: clade-specific optima ------------------------------------
#
# Simulate under BM, then add a constant shift to each of the two major
# clades descending from the root.  This produces bimodal trait distributions
# that BM (which assumes a single Gaussian process) cannot capture.

simulate_regime_shift <- function(tree, n_traits, sigma, shift_magnitude) {
  n <- length(tree$tip.label)

  # Identify the two clades from root children
  root_node <- ape::Ntip(tree) + 1L
  root_children <- tree$edge[tree$edge[, 1] == root_node, 2]

  # Get tip indices for each clade
  clade_tips <- lapply(root_children, function(node) {
    if (node <= ape::Ntip(tree)) {
      # It's a tip
      return(node)
    }
    # Get all descendants
    desc <- ape::extract.clade(tree, node)
    match(desc$tip.label, tree$tip.label)
  })

  df <- data.frame(row.names = tree$tip.label)
  for (j in seq_len(n_traits)) {
    # Start with BM
    trait <- ape::rTraitCont(tree, model = "BM", sigma = sigma,
                             root.value = 0)
    trait <- as.numeric(trait[tree$tip.label])

    # Add regime shift: clade 1 gets +shift, clade 2 gets -shift
    shift <- shift_magnitude * (1 + 0.2 * (j - 1))  # vary slightly per trait
    trait[clade_tips[[1]]] <- trait[clade_tips[[1]]] + shift
    if (length(clade_tips) >= 2L) {
      trait[clade_tips[[2]]] <- trait[clade_tips[[2]]] - shift
    }
    df[[paste0("trait", j)]] <- trait
  }
  df
}


# ---- Non-linear: quadratic and exponential correlations --------------------
#
# Generate latent BM traits, then create response variables via non-linear
# transforms.  BM assumes linear covariance structure, so the baseline
# cannot capture these relationships.

simulate_nonlinear <- function(tree, n_traits, sigma) {
  n <- length(tree$tip.label)
  n_latent <- max(2L, n_traits - 1L)  # at least 2 latent BM traits

  # Generate latent BM traits
  latent <- matrix(NA_real_, n, n_latent)
  for (j in seq_len(n_latent)) {
    trait <- ape::rTraitCont(tree, model = "BM", sigma = sigma,
                             root.value = 0)
    latent[, j] <- as.numeric(trait[tree$tip.label])
  }

  df <- data.frame(row.names = tree$tip.label)

  # First trait: quadratic function of latent[,1]
  df$trait1 <- latent[, 1]^2 + stats::rnorm(n, 0, 0.3)

  # Second trait: exponential of latent[,2]
  df$trait2 <- exp(latent[, 2] / 2) + stats::rnorm(n, 0, 0.3)

  # Third trait: interaction of latent[,1] and latent[,2]
  if (n_traits >= 3L) {
    df$trait3 <- latent[, 1] * latent[, 2] + stats::rnorm(n, 0, 0.3)
  }

  # Remaining traits: raw BM (for comparison)
  if (n_traits >= 4L) {
    for (j in 4:n_traits) {
      idx <- min(j - 1L, n_latent)
      df[[paste0("trait", j)]] <- latent[, idx]
    }
  }

  df
}


# ---- Binary traits: latent BM → threshold → factor --------------------------
#
# Generates binary traits with controllable phylogenetic signal.
# A latent BM process provides phylogenetic structure; independent noise
# dilutes it according to `signal`.  The latent value is thresholded at a
# quantile (default 0.5 = balanced classes).
#
# @param tree   phylo object.
# @param n_traits Integer.
# @param signal  Numeric in (0, 1].  1 = pure phylogenetic; 0.2 = mostly noise.
# @param threshold_quantile Numeric in (0, 1).  0.5 = balanced; 0.9 = rare "yes".
# @param seed   Integer or NULL.
# @return data.frame with factor columns (levels "no", "yes").
# @noRd
simulate_binary_traits <- function(tree, n_traits = 4L, signal = 0.6,
                                   threshold_quantile = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- tree$tip.label
  n  <- length(sp)
  df <- data.frame(row.names = sp)

  for (j in seq_len(n_traits)) {
    latent <- ape::rTraitCont(tree, model = "BM", sigma = 1.0,
                              root.value = (j - 1) * 0.3)
    latent <- as.numeric(latent[sp])
    # Add noise to control phylogenetic signal
    noise_sd <- sqrt((1 - signal) / max(signal, 1e-6)) * stats::sd(latent)
    noisy <- latent + stats::rnorm(n, 0, noise_sd)
    # Threshold at the specified quantile
    thresh <- stats::quantile(noisy, probs = threshold_quantile)
    vals <- ifelse(noisy > thresh, "yes", "no")
    df[[paste0("bin", j)]] <- factor(vals, levels = c("no", "yes"))
  }
  df
}


# ---- Ordinal traits: latent BM → noise → bin into levels --------------------
#
# @param tree   phylo object.
# @param n_traits Integer.
# @param n_levels Integer.  Number of ordered categories (e.g. 3, 5, 7, 10).
# @param signal  Numeric in (0, 1].
# @param seed   Integer or NULL.
# @return data.frame with ordered factor columns.
# @noRd
simulate_ordinal_traits <- function(tree, n_traits = 3L, n_levels = 5L,
                                    signal = 0.6, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- tree$tip.label
  n  <- length(sp)
  df <- data.frame(row.names = sp)
  levs <- paste0("L", seq_len(n_levels))

  for (j in seq_len(n_traits)) {
    latent <- ape::rTraitCont(tree, model = "BM", sigma = 1.0,
                              root.value = (j - 1) * 0.2)
    latent <- as.numeric(latent[sp])
    noise_sd <- sqrt((1 - signal) / max(signal, 1e-6)) * stats::sd(latent)
    noisy <- latent + stats::rnorm(n, 0, noise_sd)
    # Bin into n_levels equally-spaced quantiles
    breaks <- stats::quantile(noisy,
                              probs = seq(0, 1, length.out = n_levels + 1L))
    breaks[1] <- breaks[1] - 1  # ensure min is captured
    breaks[n_levels + 1L] <- breaks[n_levels + 1L] + 1
    bin_idx <- as.integer(cut(noisy, breaks = breaks, include.lowest = TRUE))
    bin_idx <- pmin(pmax(bin_idx, 1L), n_levels)
    df[[paste0("ord", j)]] <- ordered(levs[bin_idx], levels = levs)
  }
  df
}


# ---- Count traits: latent BM → exp → Poisson / NegBin -----------------------
#
# @param tree   phylo object.
# @param n_traits Integer.
# @param mean_count Numeric.  Target median count (~).
# @param overdispersion Numeric or NULL.  NULL = Poisson; numeric = NegBin size.
# @param seed   Integer or NULL.
# @return data.frame with integer columns.
# @noRd
simulate_count_traits <- function(tree, n_traits = 3L, mean_count = 20,
                                  overdispersion = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- tree$tip.label
  n  <- length(sp)
  df <- data.frame(row.names = sp)

  # Calibrate BM sigma and root value so median lambda ~ mean_count
  root_val <- log(mean_count)
  bm_sigma <- 0.5  # moderate phylogenetic variation on log scale

  for (j in seq_len(n_traits)) {
    latent <- ape::rTraitCont(tree, model = "BM", sigma = bm_sigma,
                              root.value = root_val + (j - 1) * 0.2)
    lambda <- exp(as.numeric(latent[sp]))
    lambda <- pmax(lambda, 1e-4)  # floor to avoid zero rates

    if (is.null(overdispersion)) {
      counts <- stats::rpois(n, lambda)
    } else {
      counts <- stats::rnbinom(n, mu = lambda, size = overdispersion)
    }
    df[[paste0("count", j)]] <- as.integer(counts)
  }
  df
}


# ---- Categorical traits: K latent BMs → argmax → factor ---------------------
#
# @param tree   phylo object.
# @param n_traits Integer.
# @param n_levels Integer.  Number of unordered categories K.
# @param signal  Numeric in (0, 1].
# @param seed   Integer or NULL.
# @return data.frame with factor columns.
# @noRd
simulate_categorical_traits <- function(tree, n_traits = 2L, n_levels = 3L,
                                        signal = 0.6, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- tree$tip.label
  n  <- length(sp)
  df <- data.frame(row.names = sp)
  levs <- LETTERS[seq_len(n_levels)]

  for (j in seq_len(n_traits)) {
    # Generate K latent BM processes
    latent_mat <- matrix(NA_real_, n, n_levels)
    for (k in seq_len(n_levels)) {
      lat <- ape::rTraitCont(tree, model = "BM", sigma = 1.0,
                             root.value = 0)
      latent_mat[, k] <- as.numeric(lat[sp])
    }
    # Add noise to control phylogenetic signal
    noise_sd <- sqrt((1 - signal) / max(signal, 1e-6))
    noisy <- latent_mat + matrix(stats::rnorm(n * n_levels, 0, noise_sd),
                                 n, n_levels)
    # Argmax across K → category
    cat_idx <- apply(noisy, 1, which.max)
    df[[paste0("cat", j)]] <- factor(levs[cat_idx], levels = levs)
  }
  df
}


# ---- Proportion traits: latent BM → logistic → (0,1) ------------------------
#
# @param tree   phylo object.
# @param n_traits Integer.
# @param signal  Numeric in (0, 1].
# @param boundary_frac Numeric.  Fraction of values pushed toward 0 or 1.
# @param seed   Integer or NULL.
# @return data.frame with numeric columns in (0, 1).
# @noRd
simulate_proportion_traits <- function(tree, n_traits = 4L, signal = 0.6,
                                       boundary_frac = 0.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- tree$tip.label
  n  <- length(sp)
  df <- data.frame(row.names = sp)

  for (j in seq_len(n_traits)) {
    latent <- ape::rTraitCont(tree, model = "BM", sigma = 1.0,
                              root.value = (j - 1) * 0.3)
    latent <- as.numeric(latent[sp])
    noise_sd <- sqrt((1 - signal) / max(signal, 1e-6)) * stats::sd(latent)
    noisy <- latent + stats::rnorm(n, 0, noise_sd)
    # Logistic transform to (0, 1)
    props <- stats::plogis(noisy)

    # Optionally push some values toward boundaries
    if (boundary_frac > 0) {
      n_boundary <- max(1L, round(n * boundary_frac))
      boundary_idx <- sample.int(n, n_boundary)
      # Push half toward 0, half toward 1
      half <- ceiling(n_boundary / 2)
      props[boundary_idx[seq_len(half)]] <-
        stats::runif(half, 0.001, 0.05)
      if (n_boundary > half) {
        props[boundary_idx[(half + 1):n_boundary]] <-
          stats::runif(n_boundary - half, 0.95, 0.999)
      }
    }
    # Clamp to strict (0, 1) for logit safety
    props <- pmin(pmax(props, 0.001), 0.999)
    df[[paste0("prop", j)]] <- props
  }
  df
}


# ---- Zero-inflated count traits: BM + zero gate → ZI Poisson/NegBin ---------
#
# @param tree   phylo object.
# @param n_traits Integer.
# @param zero_frac Numeric in (0, 1).  Target fraction of zeros.
# @param mean_nz Numeric.  Target mean of non-zero counts.
# @param overdispersion Numeric or NULL.  NULL = Poisson; numeric = NegBin size.
# @param seed   Integer or NULL.
# @return data.frame with integer columns (many zeros).
# @noRd
simulate_zi_count_traits <- function(tree, n_traits = 3L, zero_frac = 0.5,
                                     mean_nz = 20, overdispersion = NULL,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- tree$tip.label
  n  <- length(sp)
  df <- data.frame(row.names = sp)

  for (j in seq_len(n_traits)) {
    # Zero-gate: latent BM → standardise → threshold at qnorm(zero_frac)
    gate_latent <- ape::rTraitCont(tree, model = "BM", sigma = 1.0,
                                   root.value = 0)
    gate_latent <- as.numeric(gate_latent[sp])
    # Standardise to N(0,1) so qnorm threshold controls the zero fraction
    gate_z <- (gate_latent - mean(gate_latent)) / max(stats::sd(gate_latent), 1e-6)
    is_zero <- gate_z < stats::qnorm(zero_frac)

    # Count magnitude: independent latent BM → exp → Poisson/NegBin
    count_latent <- ape::rTraitCont(tree, model = "BM", sigma = 0.5,
                                    root.value = log(mean_nz) + (j - 1) * 0.1)
    lambda <- exp(as.numeric(count_latent[sp]))
    lambda <- pmax(lambda, 1e-4)

    if (is.null(overdispersion)) {
      counts <- stats::rpois(n, lambda)
    } else {
      counts <- stats::rnbinom(n, mu = lambda, size = overdispersion)
    }
    counts <- pmax(counts, 1L)  # ensure non-zero counts are >= 1
    counts[is_zero] <- 0L

    df[[paste0("zi", j)]] <- as.integer(counts)
  }
  df
}


# ---- Multi-proportion traits: correlated BM on CLR space → softmax ---------
#
# Simulates compositional (K-category) data. Evolution is modelled as
# correlated Brownian motion on K-1 dimensional CLR space (the K CLR
# coordinates are constrained to sum to zero, so we simulate K independent
# BMs and project onto the sum-zero hyperplane, which is equivalent).
#
# @param tree   phylo object.
# @param K      Number of categories (integer >= 2).
# @param signal Numeric in (0, 1): proportion of variance from phylogeny.
# @param sigma  BM rate on CLR scale.
# @param seed   Integer or NULL.
# @return data.frame with K numeric columns (names cat1, ..., catK),
#         rows summing to 1. Rownames = tree$tip.label.
# @noRd
simulate_multi_proportion_traits <- function(tree, K = 5L, signal = 0.6,
                                             sigma = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- tree$tip.label
  n  <- length(sp)

  # Simulate K independent BMs
  clr_mat <- matrix(0, n, K)
  for (k in seq_len(K)) {
    tr <- ape::rTraitCont(tree, model = "BM", sigma = sigma,
                          root.value = 0)
    clr_mat[, k] <- as.numeric(tr[sp])
  }

  # Add Gaussian noise to control phylogenetic signal
  if (signal < 1) {
    noise_sd <- sqrt((1 - signal) / max(signal, 1e-6)) * stats::sd(clr_mat)
    clr_mat <- clr_mat + matrix(stats::rnorm(n * K, 0, noise_sd), n, K)
  }

  # Project to sum-zero (CLR constraint)
  clr_mat <- clr_mat - rowMeans(clr_mat)

  # Softmax to simplex
  ex <- exp(clr_mat - apply(clr_mat, 1, max))
  prop <- ex / rowSums(ex)

  df <- as.data.frame(prop)
  names(df) <- paste0("cat", seq_len(K))
  rownames(df) <- sp
  df
}
