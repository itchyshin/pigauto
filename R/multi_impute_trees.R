#' Propagate phylogenetic uncertainty through multiple imputation
#'
#' Run pigauto's imputation pipeline on each of `T` posterior phylogenies,
#' generating `m_per_tree` stochastic completions per tree for a total of
#' `T * m_per_tree` completed datasets.
#'
#' The key insight (Nakagawa & de Villemereuil 2019) is that phylogenetic
#' uncertainty can be folded into the multiple-imputation framework by
#' treating each (tree, imputation) pair as one completed dataset and
#' applying Rubin's rules across all of them. The resulting pooled standard
#' errors propagate *both* imputation uncertainty (missing trait values) and
#' phylogenetic uncertainty (unknown tree topology and branch lengths).
#'
#' @param traits data.frame. Same format as [multi_impute()] and [impute()].
#' @param trees list of `phylo` objects (class `multiPhylo` or plain list).
#'   Each tree must contain the species in `traits` as tips. Posterior
#'   samples from BirdTree.org (Jetz et al. 2012) are ideal; the bundled
#'   [trees300] dataset provides 50 posterior trees for [avonet300].
#' @param m_per_tree integer. Number of MC-dropout imputations per tree
#'   (default `5`). Total datasets = `length(trees) * m_per_tree`.
#' @param species_col character or `NULL`. See [impute()].
#' @param log_transform logical. Auto-log positive continuous columns
#'   (default `TRUE`).
#' @param missing_frac numeric. Fraction held out for validation/test
#'   during training (default `0.25`).
#' @param covariates data.frame or matrix of environmental covariates
#'   (fully observed, numeric). Passed through to [impute()].
#'   Default `NULL` (no covariates).
#' @param epochs integer. Maximum GNN training epochs per tree (default
#'   `2000`).
#' @param verbose logical. Print progress (default `TRUE`).
#' @param seed integer. Base random seed; each tree uses `seed + t - 1`
#'   so results are reproducible (default `1`).
#' @param ... additional arguments forwarded to [fit_pigauto()] via
#'   [impute()].
#'
#' @return An object of class `"pigauto_mi_trees"`, inheriting from
#'   `"pigauto_mi"`, with components:
#'   \describe{
#'     \item{`datasets`}{List of `T * m_per_tree` completed data.frames.
#'       Observed cells are preserved; missing cells are filled with
#'       imputation draws. Compatible with [with_imputations()].}
#'     \item{`m`}{Total number of datasets (`T * m_per_tree`).}
#'     \item{`n_trees`}{Number of posterior trees used.}
#'     \item{`m_per_tree`}{Imputations per tree.}
#'     \item{`tree_index`}{Integer vector of length `m`; element `i` gives
#'       the tree index (1..T) for dataset `i`.}
#'     \item{`pooled_point`}{Single data.frame averaging across all
#'       `T * m_per_tree` datasets. For reporting, not inference.}
#'     \item{`se`}{Matrix of per-cell pooled SEs (NA if not available).}
#'     \item{`imputed_mask`}{Logical matrix; `TRUE` where a cell was
#'       originally missing.}
#'     \item{`fits`}{List of `T` \code{pigauto_fit} objects (one per tree),
#'       for diagnostics.}
#'     \item{`trees`}{The input posterior trees.}
#'     \item{`species_col`}{Passed-through species column name.}
#'   }
#'
#' @details
#' For each tree the function runs the full pigauto pipeline
#' (preprocess -> baseline -> GNN -> predict). Topologies and branch
#' lengths vary across trees, so the phylogenetic graph, baseline
#' covariance, and GNN weights differ for each tree. Each tree produces
#' `m_per_tree` stochastic completions via MC dropout.
#'
#' Downstream usage is identical to [multi_impute()]: pass the result
#' to [with_imputations()] to fit a model on each dataset, then to
#' [pool_mi()] for Rubin's-rules pooling. The pooled standard errors
#' will be wider than those from a single tree because they incorporate
#' the extra between-tree variance.
#'
#' **Variance decomposition.** The between-imputation variance from
#' Rubin's rules has two sources: (1) within-tree sampling variance
#' (MC-dropout noise), and (2) between-tree variance (phylogenetic
#' uncertainty). The fraction of missing information (FMI) reported by
#' [pool_mi()] reflects both. To decompose them, compare FMI from
#' `multi_impute()` (single tree) with FMI from `multi_impute_trees()`.
#'
#' **Computation time.** Each tree requires a full retrain (~3 s for
#' 300 species). 50 trees x 5 imputations takes ~2.5 minutes on a
#' modern laptop. Increase `m_per_tree` if the within-tree FMI is large
#' (many missing cells); increase the number of trees if the between-tree
#' FMI is large (deep phylogenetic uncertainty).
#'
#' @references
#' Nakagawa S, de Villemereuil P (2019). "A general method for
#' simultaneously accounting for phylogenetic and species sampling
#' uncertainty via Rubin's rules in comparative analysis."
#' \emph{Systematic Biology} 68(4): 632-641.
#'
#' Jetz W, Thomas GH, Joy JB, Hartmann K, Mooers AO (2012). "The
#' global diversity of birds in space and time." \emph{Nature}
#' 491(7424): 444-448.
#'
#' @seealso [multi_impute()] for single-tree MI, [with_imputations()],
#'   [pool_mi()], [trees300]
#'
#' @examples
#' \dontrun{
#' library(pigauto)
#' data(avonet300, trees300)
#' df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
#'
#' # 50 trees x 5 imputations = 250 completed datasets
#' mi_trees <- multi_impute_trees(df, trees300, m_per_tree = 5)
#' print(mi_trees)
#'
#' # Downstream: phylogenetic GLS, pooled with Rubin's rules
#' fits <- with_imputations(mi_trees, function(d) {
#'   d$species <- rownames(d)
#'   nlme::gls(
#'     log(Mass) ~ log(Wing.Length),
#'     correlation = ape::corBrownian(
#'       phy = trees300[[1]], form = ~species),
#'     data = d, method = "ML"
#'   )
#' })
#' pool_mi(fits)
#' }
#'
#' @export
multi_impute_trees <- function(traits, trees, m_per_tree = 5L,
                               species_col = NULL,
                               log_transform = TRUE,
                               missing_frac = 0.25,
                               covariates = NULL,
                               epochs = 2000L, verbose = TRUE,
                               seed = 1L, ...) {

  # ---- Validate inputs -------------------------------------------------------
  m_per_tree <- as.integer(m_per_tree)
  if (!is.finite(m_per_tree) || m_per_tree < 1L) {
    stop("`m_per_tree` must be a positive integer. Got ", m_per_tree, ".",
         call. = FALSE)
  }

  if (!is.list(trees)) {
    stop("`trees` must be a list of phylo objects (multiPhylo or plain list).",
         call. = FALSE)
  }
  T_trees <- length(trees)
  if (T_trees < 2L) {
    stop("`trees` must contain at least 2 phylogenies for tree uncertainty ",
         "propagation. Got ", T_trees, ". For single-tree MI, use ",
         "`multi_impute()` instead.", call. = FALSE)
  }

  # Verify all elements are phylo objects
  is_phylo <- vapply(trees, inherits, logical(1), "phylo")
  if (!all(is_phylo)) {
    bad <- which(!is_phylo)
    stop("All elements of `trees` must be 'phylo' objects. Non-phylo at ",
         "indices: ", paste(head(bad, 5), collapse = ", "),
         if (length(bad) > 5) ", ...", call. = FALSE)
  }

  M_total <- T_trees * m_per_tree

  if (verbose) {
    cat(sprintf("multi_impute_trees: %d trees x %d imputations = %d datasets\n",
                T_trees, m_per_tree, M_total))
  }

  # ---- Run impute() on each tree ---------------------------------------------
  all_datasets  <- vector("list", M_total)
  tree_index    <- integer(M_total)
  all_fits      <- vector("list", T_trees)
  imputed_mask  <- NULL
  se_sum        <- NULL
  pooled_sum    <- NULL

  for (t in seq_len(T_trees)) {
    t_seed <- as.integer(seed + t - 1L)

    if (verbose) {
      cat(sprintf("  Tree %d/%d (seed=%d)...", t, T_trees, t_seed))
      t_start <- proc.time()
    }

    # Run full pipeline
    res <- impute(
      traits        = traits,
      tree          = trees[[t]],
      species_col   = species_col,
      log_transform = log_transform,
      missing_frac  = missing_frac,
      n_imputations = m_per_tree,
      covariates    = covariates,
      epochs        = as.integer(epochs),
      verbose       = FALSE,
      seed          = t_seed,
      ...
    )

    pred <- res$prediction

    # Verify we got the right number of imputation datasets
    if (is.null(pred$imputed_datasets) || length(pred$imputed_datasets) != m_per_tree) {
      stop("Tree ", t, ": predict() did not return ", m_per_tree,
           " imputed datasets. This is an internal error.", call. = FALSE)
    }

    # Build completed data.frames
    for (k in seq_len(m_per_tree)) {
      idx <- (t - 1L) * m_per_tree + k
      completed_info <- build_completed(traits, pred$imputed_datasets[[k]],
                                        species_col)
      all_datasets[[idx]] <- completed_info$completed
      tree_index[idx]     <- t
    }

    # Store the fit for diagnostics
    all_fits[[t]] <- res$fit

    # Accumulate imputed_mask (should be same across all trees)
    if (is.null(imputed_mask)) {
      imputed_mask <- res$imputed_mask
    }

    # Accumulate for pooled point estimate and SE
    trait_cols <- setdiff(names(traits), species_col)
    if (is.null(pooled_sum)) {
      pooled_sum <- res$completed
      # Initialise numeric columns to 0 for summing
      for (nm in trait_cols) {
        if (is.numeric(pooled_sum[[nm]])) {
          pooled_sum[[nm]] <- res$completed[[nm]]
        }
      }
    } else {
      for (nm in trait_cols) {
        if (is.numeric(res$completed[[nm]])) {
          pooled_sum[[nm]] <- pooled_sum[[nm]] + res$completed[[nm]]
        }
      }
    }

    if (verbose) {
      elapsed <- (proc.time() - t_start)[3]
      cat(sprintf(" %.1fs\n", elapsed))
    }
  }

  # ---- Pooled point estimate (average across trees) --------------------------
  pooled_point <- pooled_sum
  for (nm in trait_cols) {
    if (is.numeric(pooled_point[[nm]])) {
      pooled_point[[nm]] <- pooled_point[[nm]] / T_trees
    }
  }
  # For non-numeric columns, use the result from the first tree (majority vote
  # would be more principled, but for the point estimate this is secondary)

  # ---- Assemble output -------------------------------------------------------
  structure(
    list(
      datasets     = all_datasets,
      m            = M_total,
      n_trees      = T_trees,
      m_per_tree   = m_per_tree,
      tree_index   = tree_index,
      pooled_point = pooled_point,
      se           = NULL,   # per-cell SE not meaningful across trees
      imputed_mask = imputed_mask,
      fits         = all_fits,
      trees        = trees,
      species_col  = species_col
    ),
    class = c("pigauto_mi_trees", "pigauto_mi", "list")
  )
}


#' @export
print.pigauto_mi_trees <- function(x, ...) {
  T_trees    <- x$n_trees
  m_per_tree <- x$m_per_tree
  M_total    <- x$m

  # Count species/traits from the first dataset
  d1 <- x$datasets[[1]]
  trait_cols <- setdiff(names(d1), x$species_col)
  n_sp <- nrow(d1)

  total_cells <- length(x$imputed_mask)
  n_imp_cells <- sum(x$imputed_mask)
  pct <- if (total_cells > 0) 100 * n_imp_cells / total_cells else 0

  cat("pigauto multiple imputation with tree uncertainty\n")
  cat(sprintf("  Trees     : %d posterior phylogenies\n", T_trees))
  cat(sprintf("  Per tree  : %d imputations (MC dropout)\n", m_per_tree))
  cat(sprintf("  Total     : %d completed datasets\n", M_total))
  cat(sprintf("  Species   : %d\n", n_sp))
  cat(sprintf("  Traits    : %d -- %s\n", length(trait_cols),
              paste(trait_cols, collapse = ", ")))
  cat(sprintf("  Cells     : %d imputed / %d total (%.1f%%)\n",
              n_imp_cells, total_cells, pct))

  cat("\n  Access imputation draws:  mi$datasets[[i]]\n")
  cat("  Tree index for draw i:    mi$tree_index[i]\n")
  cat("  Fit downstream models:    with_imputations(mi, fit_fun)\n")
  cat("  Pool with Rubin's rules:  pool_mi(fits)\n")
  invisible(x)
}
