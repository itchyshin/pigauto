#' Fit a BACE (Bayesian Augmentation using Chained Equations) baseline
#'
#' Uses \pkg{BACE} to provide a phylogenetically-informed Bayesian baseline
#' for all trait types.  Returns imputed means and between-imputation SEs
#' in latent scale, matching the interface of \code{\link{fit_baseline}}.
#'
#' @details
#' BACE runs chained MCMCglmm imputation: each trait is modelled as a
#' response with all others as predictors, cycling through multiple MCMC
#' runs.  This provides a fully phylogenetic baseline for binary,
#' categorical, and count traits --- not just continuous ones.
#'
#' The returned \code{mu} matrix is the mean across imputation runs;
#' \code{se} is the between-imputation SD (capturing imputation
#' uncertainty).  Both are in latent scale (same as
#' \code{pigauto_data$X_scaled}), so they can be passed directly to
#' \code{\link{fit_pigauto}} as the \code{baseline} argument.
#'
#' @section Which BACE datasets \code{mu} / \code{se} are built from:
#' Two paths are available and they differ in what the returned
#' \code{se} means.
#'
#' \describe{
#'   \item{\code{final_imp = FALSE} (default)}{\code{mu} / \code{se} are
#'     computed from the \code{runs} datasets of the \code{bace_imp()}
#'     chain.  Those are successive sweeps of one chained-equations
#'     chain, so they are autocorrelated by construction and still carry
#'     convergence transient.  The resulting \code{se} is a cheap
#'     dispersion summary, not a calibrated imputation standard error,
#'     and it has no coverage guarantee.}
#'   \item{\code{final_imp = TRUE}}{\code{BACE::bace_final_imp()} is run
#'     on the fitted chain object and \code{mu} / \code{se} are computed
#'     from its \code{n_final} final datasets.  Each of those runs starts
#'     independently from the converged chain, so they are proper
#'     multiple-imputation draws and \code{se} is a between-imputation SD
#'     in the Rubin (1987) sense.  Expect this \code{se} to be
#'     \emph{larger} than the default path's.  Cost is \code{n_final}
#'     extra MCMC fits per trait on top of the chain.}
#' }
#'
#' The default is \code{FALSE} so that existing calls are unchanged.
#' Note that this is a choice about which BACE datasets pigauto
#' summarises; it is not a correction to the default path's arithmetic.
#'
#' \code{BACE::bace_final_imp()} is only available in recent BACE
#' versions.  If it is missing, \code{final_imp = TRUE} errors with an
#' upgrade hint rather than silently falling back.
#'
#' @param data object of class \code{"pigauto_data"}.
#' @param tree object of class \code{"phylo"}.
#' @param splits list (output of \code{\link{make_missing_splits}}) or
#'   \code{NULL}.
#' @param runs integer. Number of BACE chained imputation iterations
#'   (default \code{5L}).
#' @param nitt integer. MCMC iterations per model (default \code{4000L}).
#' @param burnin integer. Burn-in iterations (default \code{1000L}).
#' @param thin integer. Thinning rate (default \code{10L}).
#' @param final_imp logical. If \code{TRUE}, append
#'   \code{BACE::bace_final_imp()} after the chain and build \code{mu} /
#'   \code{se} from its final datasets instead of the chain datasets.
#'   Default \code{FALSE} (unchanged behaviour).
#' @param n_final integer. Number of final imputation draws when
#'   \code{final_imp = TRUE} (default \code{15L}, matching the BACE
#'   simulation study).  Ignored when \code{final_imp = FALSE}.  Small
#'   values undercover; BACE's own default is 50.
#' @param verbose logical.
#' @return A list with:
#'   \describe{
#'     \item{mu}{Numeric matrix (n x p_latent), baseline means.}
#'     \item{se}{Numeric matrix (n x p_latent), between-imputation SEs.}
#'   }
#' @examples
#' \dontrun{
#' bl_bace <- fit_baseline_bace(pd, tree, splits = spl)
#' fit <- fit_pigauto(pd, tree, splits = spl, baseline = bl_bace)
#' }
#' @export
fit_baseline_bace <- function(data, tree, splits = NULL,
                              runs = 5L, nitt = 4000L,
                              burnin = 1000L, thin = 10L,
                              final_imp = FALSE, n_final = 15L,
                              verbose = TRUE) {

  if (!requireNamespace("BACE", quietly = TRUE)) {
    stop("BACE package required but not installed.\n",
         "  Install from local source: devtools::install('path/to/BACE')")
  }
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object.")
  }

  if (!is.logical(final_imp) || length(final_imp) != 1L || is.na(final_imp)) {
    stop("'final_imp' must be a single TRUE or FALSE.")
  }
  if (final_imp) {
    n_final <- as.integer(n_final)
    if (length(n_final) != 1L || is.na(n_final) || n_final < 2L) {
      stop("'n_final' must be a single integer >= 2 when final_imp = TRUE.")
    }
    if (!"bace_final_imp" %in% getNamespaceExports("BACE")) {
      stop("final_imp = TRUE requires a BACE version exporting ",
           "bace_final_imp().\n",
           "  Update BACE, or call with final_imp = FALSE.")
    }
  }

  trait_map   <- data$trait_map
  trait_names <- data$trait_names
  sp_names    <- data$species_names

  if (is.null(trait_map)) {
    stop("trait_map required. Use preprocess_traits() to produce pigauto_data.")
  }

  # ---- Reconstruct original-scale data and mask val/test cells ---------------
  df <- data$X_original
  n  <- nrow(df)

  if (!is.null(splits)) {
    if (!is.null(splits$val_idx_trait) && !is.null(splits$test_idx_trait)) {
      idx_trait <- c(splits$val_idx_trait, splits$test_idx_trait)
    } else {
      # Convert latent indices to original-trait indices
      idx_trait <- latent_idx_to_trait_idx(
        c(splits$val_idx, splits$test_idx), n, trait_map
      )
    }
    row_i <- ((idx_trait - 1L) %% n) + 1L
    col_j <- ceiling(idx_trait / n)
    for (k in seq_along(idx_trait)) {
      if (col_j[k] >= 1L && col_j[k] <= ncol(df)) {
        df[row_i[k], col_j[k]] <- NA
      }
    }
  }

  # ---- Add species column and build formula ----------------------------------
  df$species <- sp_names

  formula_str <- paste(
    trait_names[1], "~",
    paste(trait_names[-1], collapse = " + ")
  )

  # ---- Run BACE imputation ---------------------------------------------------
  fit <- BACE::bace_imp(
    fixformula     = formula_str,
    ran_phylo_form = "~ 1 | species",
    phylo          = tree,
    data           = df,
    runs           = as.integer(runs),
    nitt           = as.integer(nitt),
    burnin         = as.integer(burnin),
    thin           = as.integer(thin),
    verbose        = verbose
  )

  # ---- Choose which imputed datasets mu / se summarise -----------------------
  if (final_imp) {
    # Proper-MI path: n_final independent draws from the converged chain.
    final_fit <- BACE::bace_final_imp(
      bace_object    = fit,
      fixformula     = formula_str,
      ran_phylo_form = "~ 1 | species",
      phylo          = tree,
      nitt           = as.integer(nitt),
      burnin         = as.integer(burnin),
      thin           = as.integer(thin),
      n_final        = n_final,
      verbose        = verbose
    )
    imp_datasets <- final_fit$all_datasets
    if (!length(imp_datasets)) {
      stop("BACE::bace_final_imp() returned no datasets.")
    }
  } else {
    # Default path: chain datasets (skip Initial_Data).
    n_iter      <- length(fit$data) - 1L
    imp_datasets <- fit$data[seq(2, n_iter + 1)]
  }

  # ---- Re-encode each imputed dataset to latent scale ------------------------
  latent_mats <- lapply(imp_datasets, function(imp_df) {
    imp_df <- as.data.frame(imp_df)

    # Remove non-trait columns
    keep <- intersect(trait_names, names(imp_df))
    imp_df <- imp_df[, keep, drop = FALSE]

    # Coerce column types to match trait_map expectations
    for (tm in trait_map) {
      nm <- tm$name
      if (!nm %in% names(imp_df)) next

      if (tm$type == "binary") {
        if (!is.factor(imp_df[[nm]])) {
          imp_df[[nm]] <- factor(
            ifelse(as.numeric(imp_df[[nm]]) >= 0.5,
                   tm$levels[2], tm$levels[1]),
            levels = tm$levels
          )
        } else {
          imp_df[[nm]] <- factor(imp_df[[nm]], levels = tm$levels)
        }

      } else if (tm$type == "categorical") {
        if (!is.factor(imp_df[[nm]])) {
          imp_df[[nm]] <- factor(as.character(imp_df[[nm]]),
                                 levels = tm$levels)
        } else {
          imp_df[[nm]] <- factor(imp_df[[nm]], levels = tm$levels)
        }

      } else if (tm$type == "ordinal") {
        if (!is.ordered(imp_df[[nm]])) {
          if (is.numeric(imp_df[[nm]])) {
            # BACE ordinal: integer 1..K
            int_val <- pmin(pmax(round(imp_df[[nm]]), 1L),
                            length(tm$levels))
            imp_df[[nm]] <- ordered(tm$levels[int_val],
                                    levels = tm$levels)
          } else {
            imp_df[[nm]] <- ordered(imp_df[[nm]], levels = tm$levels)
          }
        }

      } else if (tm$type == "count") {
        imp_df[[nm]] <- as.integer(pmax(round(as.numeric(imp_df[[nm]])), 0L))

      } else if (tm$type == "continuous") {
        imp_df[[nm]] <- as.numeric(imp_df[[nm]])
      }
    }

    rownames(imp_df) <- sp_names
    encode_to_latent(imp_df, trait_map)$X
  })

  # ---- Compute mu (mean) and se (between-imputation SD) ----------------------
  p_latent <- ncol(latent_mats[[1]])
  mu <- Reduce("+", latent_mats) / length(latent_mats)

  se <- matrix(0, nrow = n, ncol = p_latent)
  if (length(latent_mats) > 1L) {
    for (j in seq_len(p_latent)) {
      vals <- sapply(latent_mats, function(m) m[, j])
      if (is.matrix(vals)) {
        se[, j] <- apply(vals, 1, stats::sd, na.rm = TRUE)
      }
    }
  }

  rownames(mu) <- sp_names
  colnames(mu) <- colnames(latent_mats[[1]])
  rownames(se) <- sp_names
  colnames(se) <- colnames(latent_mats[[1]])

  # Replace any NaN from encoding with 0
 mu[!is.finite(mu)] <- 0
  se[!is.finite(se)] <- 0

  list(mu = mu, se = se)
}


# ---- Internal: convert latent indices to original-trait indices ----------------

latent_idx_to_trait_idx <- function(latent_idx, n, trait_map) {
  row_i <- ((latent_idx - 1L) %% n) + 1L
  col_l <- ceiling(latent_idx / n)

  # Map latent column -> original trait number
  p_latent <- sum(vapply(trait_map, "[[", integer(1), "n_latent"))
  col_to_trait <- integer(p_latent)
  for (i in seq_along(trait_map)) {
    for (lc in trait_map[[i]]$latent_cols) {
      col_to_trait[lc] <- i
    }
  }

  trait_j <- col_to_trait[col_l]

  # Deduplicate (categorical traits have K latent cols but 1 original col)
  pairs <- unique(data.frame(row = row_i, trait = trait_j))
  pairs$row + (pairs$trait - 1L) * n
}
