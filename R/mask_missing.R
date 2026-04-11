#' Create an observed/missing mask matrix
#'
#' Returns a logical matrix of the same dimensions as \code{X}, with
#' \code{TRUE} where values are observed (not \code{NA}).
#'
#' @param X numeric matrix (species x traits or species x latent columns).
#' @return Logical matrix, \code{TRUE} = observed.
#' @examples
#' X <- matrix(c(1, NA, 3, 4), nrow = 2)
#' mask_missing(X)
#' @export
mask_missing <- function(X) {
  !is.na(X)
}

#' Split cells into train/val/test for imputation evaluation
#'
#' Randomly designates a fraction of cells as "missing" and splits them into
#' validation and test sets.  When a \code{trait_map} is supplied, masking
#' operates at the **original trait level** -- all latent columns belonging
#' to one trait are held out together (important for categorical traits).
#'
#' The returned index vectors use linear (column-major) indexing.  Both
#' original-trait-space and latent-space indices are returned when a
#' \code{trait_map} is present.
#'
#' @param X numeric matrix (species x latent columns from
#'   \code{\link{preprocess_traits}}).  Used only for dimensions.
#' @param missing_frac numeric. Fraction of all (species, trait) cells to
#'   designate as missing (default \code{0.25}).
#' @param val_frac numeric. Fraction of missing cells for validation
#'   (default \code{0.25}); the rest become the test set.
#' @param seed integer. Random seed for reproducibility (default \code{555}).
#' @param trait_map list of trait descriptors (from \code{pigauto_data}).
#'   If \code{NULL}, masking is applied per latent column (v0.1 behaviour).
#' @param mechanism character. Missingness mechanism:
#'   \code{"MCAR"} (default, uniform random),
#'   \code{"MAR_trait"} (trait-dependent),
#'   \code{"MAR_phylo"} (clade-structured), or
#'   \code{"MNAR"} (value-dependent).
#' @param mechanism_args named list of mechanism-specific parameters:
#'   \describe{
#'     \item{For \code{"MAR_trait"}:}{
#'       \code{driver_col} (integer, column index in \code{X} that drives
#'         missingness; default 1),
#'       \code{beta} (numeric, severity; default 2.0).
#'     }
#'     \item{For \code{"MAR_phylo"}:}{
#'       \code{n_clades} (integer, number of high-missingness clades;
#'         default 2),
#'       \code{p_clade} (numeric, within-clade missingness probability;
#'         default 0.7),
#'       \code{p_base} (numeric, background missingness probability;
#'         default 0.1).
#'     }
#'     \item{For \code{"MNAR"}:}{
#'       \code{beta} (numeric, severity; default 2.0).
#'     }
#'   }
#' @param tree object of class \code{"phylo"}.  Required for
#'   \code{mechanism = "MAR_phylo"}, ignored otherwise.
#' @return A list with:
#'   \describe{
#'     \item{val_idx}{Integer vector of linear indices (latent space).}
#'     \item{test_idx}{Integer vector of linear indices (latent space).}
#'     \item{val_idx_trait}{Integer vector in original-trait space (if
#'       \code{trait_map} supplied).}
#'     \item{test_idx_trait}{Integer vector in original-trait space (if
#'       \code{trait_map} supplied).}
#'     \item{n}{Number of species (rows).}
#'     \item{p}{Number of latent columns.}
#'     \item{n_traits}{Number of original traits.}
#'     \item{mask}{Logical matrix (n x p_latent). \code{TRUE} = observed.}
#'     \item{mechanism}{Character string of the mechanism used.}
#'   }
#' @examples
#' X <- matrix(rnorm(100), nrow = 20)
#' splits <- make_missing_splits(X, missing_frac = 0.25, seed = 1)
#' length(splits$val_idx)
#'
#' # MAR: missingness depends on another trait
#' splits_mar <- make_missing_splits(X, mechanism = "MAR_trait",
#'   mechanism_args = list(driver_col = 1, beta = 2))
#' @export
make_missing_splits <- function(X, missing_frac = 0.25, val_frac = 0.25,
                                seed = 555, trait_map = NULL,
                                mechanism = c("MCAR", "MAR_trait",
                                              "MAR_phylo", "MNAR"),
                                mechanism_args = list(),
                                tree = NULL) {
  mechanism <- match.arg(mechanism)
  if (!is.matrix(X)) stop("'X' must be a matrix.")
  if (missing_frac <= 0 || missing_frac >= 1)
    stop("'missing_frac' must be in (0, 1).")
  if (val_frac <= 0 || val_frac >= 1)
    stop("'val_frac' must be in (0, 1).")

  n <- nrow(X)
  p <- ncol(X)

  if (mechanism == "MAR_phylo" && is.null(tree)) {
    stop("'tree' is required for mechanism = 'MAR_phylo'.")
  }

  if (is.null(trait_map)) {
    # ---- v0.1 behaviour: mask per latent column ----------------------------
    # Only sample from cells that are actually observed (non-NA).
    observed_idx <- which(!is.na(X))
    if (length(observed_idx) == 0L) {
      stop("No observed cells in X; cannot create splits.")
    }
    set.seed(seed)
    n_miss <- floor(missing_frac * length(observed_idx))

    # Generate per-cell weights (uniform for MCAR, non-uniform for MAR/MNAR)
    probs <- generate_missing_probs(X, n, p, mechanism, mechanism_args,
                                    tree, trait_map = NULL)
    prob_vec <- if (is.null(probs)) NULL else probs[observed_idx]

    miss     <- sample(observed_idx, n_miss, prob = prob_vec)
    n_val    <- floor(val_frac * length(miss))
    val_idx  <- miss[seq_len(n_val)]
    test_idx <- miss[(n_val + 1L):length(miss)]

    mask <- matrix(TRUE, nrow = n, ncol = p)
    mask[!is.na(X)] <- TRUE  # observed cells stay TRUE
    mask[c(val_idx, test_idx)] <- FALSE

    return(list(
      val_idx   = val_idx,
      test_idx  = test_idx,
      n         = n,
      p         = p,
      n_traits  = p,
      mask      = mask,
      mechanism = mechanism
    ))
  }

  # ---- Trait-level masking (mixed types) -----------------------------------
  # Build a (n_species x n_traits) observation mask: a cell is observed
  # when ALL its latent columns are non-NA in X.
  n_traits <- length(trait_map)
  obs_trait <- matrix(TRUE, nrow = n, ncol = n_traits)
  for (j in seq_len(n_traits)) {
    tm <- trait_map[[j]]
    for (lc in tm$latent_cols) {
      obs_trait[, j] <- obs_trait[, j] & !is.na(X[, lc])
    }
  }
  observed_trait_idx <- which(obs_trait)
  if (length(observed_trait_idx) == 0L) {
    stop("No observed (species, trait) cells; cannot create splits.")
  }

  set.seed(seed)
  n_miss <- floor(missing_frac * length(observed_trait_idx))

  # Generate per-cell weights in (n x n_traits) space
  probs <- generate_missing_probs(X, n, p, mechanism, mechanism_args,
                                  tree, trait_map)
  prob_vec <- if (is.null(probs)) NULL else probs[observed_trait_idx]

  miss     <- sample(observed_trait_idx, n_miss, prob = prob_vec)
  n_val    <- floor(val_frac * length(miss))
  val_trait   <- miss[seq_len(n_val)]
  test_trait  <- miss[(n_val + 1L):length(miss)]

  # Expand original-trait indices to latent-column indices
  val_latent  <- expand_trait_idx_to_latent(val_trait, n, trait_map)
  test_latent <- expand_trait_idx_to_latent(test_trait, n, trait_map)

  mask <- matrix(TRUE, nrow = n, ncol = p)
  mask[c(val_latent, test_latent)] <- FALSE

  list(
    val_idx        = val_latent,
    test_idx       = test_latent,
    val_idx_trait  = val_trait,
    test_idx_trait = test_trait,
    n              = n,
    p              = p,
    n_traits       = n_traits,
    mask           = mask,
    mechanism      = mechanism
  )
}


#' Expand original-trait-space linear indices to latent-space indices
#'
#' @param idx integer vector of linear indices in (n x n_traits) space.
#' @param n integer, number of species.
#' @param trait_map list of trait descriptors.
#' @return integer vector of linear indices in (n x p_latent) space.
#' @keywords internal
expand_trait_idx_to_latent <- function(idx, n, trait_map) {
  # Convert linear index to (row, trait_number)
  row_i   <- ((idx - 1L) %% n) + 1L
  trait_j <- ((idx - 1L) %/% n) + 1L

  latent_idx <- integer(0)
  for (i in seq_along(idx)) {
    tm <- trait_map[[trait_j[i]]]
    for (lc in tm$latent_cols) {
      latent_idx <- c(latent_idx, row_i[i] + (lc - 1L) * n)
    }
  }
  latent_idx
}


#' Generate per-cell missingness probability weights
#'
#' Returns a numeric vector (or matrix) of sampling weights that
#' \code{make_missing_splits} uses in place of uniform probabilities.
#' For \code{"MCAR"}, returns \code{NULL} (uniform sampling).
#'
#' @param X numeric matrix (n x p_latent).
#' @param n integer, number of rows (species).
#' @param p integer, number of latent columns.
#' @param mechanism character, one of \code{"MCAR"}, \code{"MAR_trait"},
#'   \code{"MAR_phylo"}, \code{"MNAR"}.
#' @param args named list of mechanism-specific parameters.
#' @param tree phylo object (needed for \code{"MAR_phylo"}).
#' @param trait_map list of trait descriptors, or \code{NULL}.
#' @return A numeric vector of per-cell weights (same length as the
#'   appropriate index space), or \code{NULL} for MCAR.
#' @keywords internal
#' @noRd
generate_missing_probs <- function(X, n, p, mechanism, args, tree,
                                   trait_map) {
  if (mechanism == "MCAR") return(NULL)

  # Determine the working space: (n x n_traits) if trait_map, else (n x p)
  if (!is.null(trait_map)) {
    n_cols <- length(trait_map)
  } else {
    n_cols <- p
  }

  if (mechanism == "MAR_trait") {
    # ---- Trait-dependent MAR ------------------------------------------------
    # P(cell missing) increases with the value of a driver column.
    # The driver is specified as a column index in X (latent space).
    driver_col <- args$driver_col %||% 1L
    beta       <- args$beta %||% 2.0

    driver <- X[, driver_col]
    # Z-score the driver (use observed values only)
    obs <- !is.na(driver)
    if (sum(obs) < 2L) {
      # Not enough observed values in driver — fall back to MCAR
      return(NULL)
    }
    z <- (driver - mean(driver[obs], na.rm = TRUE)) /
         stats::sd(driver[obs], na.rm = TRUE)
    z[is.na(z)] <- 0  # missing driver values get neutral weight

    # Logistic transform: higher driver values → higher miss probability
    row_probs <- stats::plogis(beta * z)

    # Build (n x n_cols) probability matrix — same row weight for all traits
    prob_mat <- matrix(row_probs, nrow = n, ncol = n_cols)

  } else if (mechanism == "MAR_phylo") {
    # ---- Phylogenetic MAR ---------------------------------------------------
    # Pick random internal nodes; their descendant tips get elevated
    # missingness (p_clade), everyone else gets p_base.
    if (is.null(tree)) stop("'tree' required for MAR_phylo.")
    n_clades <- args$n_clades %||% 2L
    p_clade  <- args$p_clade %||% 0.7
    p_base   <- args$p_base %||% 0.1

    tips     <- tree$tip.label
    n_tips   <- length(tips)
    # Internal nodes in ape are numbered (n_tips + 1):(n_tips + tree$Nnode)
    internal <- (n_tips + 1L):(n_tips + tree$Nnode)

    # Pick n_clades internal nodes (not the root, which would select all)
    root_node <- n_tips + 1L
    candidates <- setdiff(internal, root_node)
    if (length(candidates) < n_clades) {
      n_clades <- length(candidates)
    }
    selected <- sample(candidates, n_clades)

    # Find descendant tips for each selected node
    clade_tips <- character(0)
    for (nd in selected) {
      desc <- tips_from_node(tree, nd)
      clade_tips <- union(clade_tips, desc)
    }

    # Row-level probabilities
    row_probs <- rep(p_base, n)
    # Match tip labels to rows of X (rownames of X if available)
    rn <- rownames(X)
    if (!is.null(rn)) {
      clade_rows <- which(rn %in% clade_tips)
    } else {
      # Assume rows are in tree$tip.label order
      clade_rows <- which(tips %in% clade_tips)
    }
    row_probs[clade_rows] <- p_clade

    prob_mat <- matrix(row_probs, nrow = n, ncol = n_cols)

  } else if (mechanism == "MNAR") {
    # ---- Missing Not At Random ----------------------------------------------
    # P(cell missing) increases with the absolute value of the cell itself.
    # Extreme values (both tails) are more likely to be missing.
    beta <- args$beta %||% 2.0

    if (!is.null(trait_map)) {
      # Work in (n x n_traits) space: use first latent col of each trait
      prob_mat <- matrix(0, nrow = n, ncol = n_cols)
      for (j in seq_len(n_cols)) {
        tm  <- trait_map[[j]]
        col <- tm$latent_cols[1]
        v   <- X[, col]
        obs <- !is.na(v)
        if (sum(obs) < 2L) {
          prob_mat[, j] <- 1  # uniform fallback
          next
        }
        z <- (v - mean(v[obs], na.rm = TRUE)) /
             stats::sd(v[obs], na.rm = TRUE)
        z[is.na(z)] <- 0
        prob_mat[, j] <- stats::plogis(beta * abs(z))
      }
    } else {
      # Per latent column
      prob_mat <- matrix(0, nrow = n, ncol = n_cols)
      for (j in seq_len(n_cols)) {
        v   <- X[, j]
        obs <- !is.na(v)
        if (sum(obs) < 2L) {
          prob_mat[, j] <- 1
          next
        }
        z <- (v - mean(v[obs], na.rm = TRUE)) /
             stats::sd(v[obs], na.rm = TRUE)
        z[is.na(z)] <- 0
        prob_mat[, j] <- stats::plogis(beta * abs(z))
      }
    }
  } else {
    stop("Unknown mechanism: ", mechanism)
  }

  # Apply a floor to prevent degenerate zero-weight cells
  prob_mat <- pmax(prob_mat, 0.01)

  # Return as a vector (column-major, matching linear indexing)
  as.vector(prob_mat)
}


#' Get descendant tip labels from an internal node
#'
#' @param tree phylo object.
#' @param node integer, internal node number.
#' @return character vector of tip labels.
#' @keywords internal
#' @noRd
tips_from_node <- function(tree, node) {
  n_tips <- length(tree$tip.label)
  # BFS/DFS to find all descendant tips
  stack <- node
  tips  <- integer(0)
  while (length(stack) > 0L) {
    current <- stack[1L]
    stack   <- stack[-1L]
    children <- tree$edge[tree$edge[, 1] == current, 2]
    for (ch in children) {
      if (ch <= n_tips) {
        tips <- c(tips, ch)
      } else {
        stack <- c(stack, ch)
      }
    }
  }
  tree$tip.label[tips]
}
