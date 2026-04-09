# Internal: non-BM trait simulation for benchmarking
#
# Generates species trait data under models that deviate from Brownian Motion,
# providing scenarios where the GNN's residual correction can add value over
# the BM baseline.
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
