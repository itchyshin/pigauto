# Pick the reference tree used to train the GNN when `share_gnn = TRUE`.
#
# - `reference_tree` non-NULL: return it unchanged.
# - `reference_tree` NULL + phangorn installed: compute MCC via
#   `phangorn::maxCladeCred(trees)`.
# - phangorn missing: warn and fall back to `trees[[1]]`.
resolve_reference_tree <- function(trees, reference_tree = NULL) {
  if (!is.null(reference_tree)) {
    if (!inherits(reference_tree, "phylo")) {
      stop("`reference_tree` must be a 'phylo' object. Got ",
           class(reference_tree)[1], ".", call. = FALSE)
    }
    return(reference_tree)
  }
  if (requireNamespace("phangorn", quietly = TRUE)) {
    tr_multi <- trees
    if (!inherits(tr_multi, "multiPhylo")) class(tr_multi) <- "multiPhylo"
    return(phangorn::maxCladeCred(tr_multi))
  }
  warning("phangorn required for MCC tree selection; falling back to ",
          "trees[[1]]. Install with install.packages('phangorn') for ",
          "MCC-based selection.", call. = FALSE)
  trees[[1]]
}

#' Tree-aware multiple imputation (step 1 of 2)
#'
#' Run pigauto's full imputation pipeline on each of `T` posterior
#' phylogenies, generating `m_per_tree` stochastic completions per tree
#' for a total of `T * m_per_tree` completed datasets. Each completed
#' dataset is conditional on a specific posterior tree
#' (recorded in `mi$tree_index`).
#'
#' **This is step 1 of the two-step workflow for propagating tree
#' uncertainty.** Step 2 — varying the *downstream analysis tree* in
#' lockstep with the imputation tree — is the user's responsibility
#' and is described in the Examples section and in
#' `vignette("tree-uncertainty")`.
#'
#' `multi_impute_trees()` handles the imputation half (step 1) cleanly:
#' every completed dataset carries a different tree's signal so that
#' between-tree variation propagates into the pooled standard errors.
#' Step 2 is where Nakagawa & de Villemereuil (2019) enters: for each
#' completed dataset, fit the downstream model (e.g. `nlme::gls()` with
#' a `corBrownian` on `trees[[ mi$tree_index[i] ]]`), then pool the
#' T × M fits with [pool_mi()].
#'
#' @section When to use this:
#'
#' pigauto provides two multiple-imputation functions. Pick based on how
#' many trees you have:
#'
#' * **One tree** (single published phylogeny, single time-calibrated tree):
#'   use [multi_impute()]. The `m` MC-dropout imputations capture model
#'   uncertainty.
#' * **Multiple posterior trees** (BirdTree samples, BEAST posterior, etc.):
#'   use [multi_impute_trees()]. Between-tree variation is added to the
#'   pooled SEs via Rubin's rules (Nakagawa & de Villemereuil 2019).
#'
#' The two functions share the same downstream API — both return objects
#' compatible with [with_imputations()] and [pool_mi()].
#'
#' @importFrom utils head
#'
#' @param traits data.frame. Same format as [multi_impute()] and [impute()].
#' @param trees list of `phylo` objects (class `multiPhylo` or plain list).
#'   Each tree must contain the species in `traits` as tips. Posterior
#'   samples from BirdTree.org (Jetz et al. 2012) are ideal; the bundled
#'   [trees300] dataset provides 50 posterior trees for [avonet300].
#' @param m_per_tree integer. Number of MC-dropout imputations per tree
#'   (default `1`). Total datasets = `length(trees) * m_per_tree`. The
#'   canonical Nakagawa & de Villemereuil (2019) workflow uses T = 50
#'   posterior trees × m_per_tree = 1 = M = 50 total datasets.
#' @param species_col character or `NULL`. See [impute()].
#' @param trait_types named character vector overriding auto-detected
#'   trait types. Required for `"proportion"` and `"zi_count"`. See
#'   [impute()] / [preprocess_traits()]. Default `NULL` (auto-detect).
#' @param multi_proportion_groups named list declaring compositional
#'   trait groups (rows summing to 1), forwarded to [impute()] /
#'   [preprocess_traits()]. Default `NULL`.
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
#' @param share_gnn logical. If `TRUE` (default), fit the GNN once on a
#'   reference tree and reuse it across all posterior trees, recomputing
#'   only the BM baseline per tree. Gives a ~10-15x speedup at n=10k.
#'   See the "Share-GNN" section below for tree-uncertainty propagation
#'   details. Set `FALSE` to fit from scratch on every tree (the pre-v0.9.1
#'   behaviour) when you need exact tree-by-tree model independence.
#' @param reference_tree optional `phylo` used as the training tree when
#'   `share_gnn = TRUE`. Default `NULL` selects the maximum-clade-credibility
#'   tree via `phangorn::maxCladeCred(trees)`. If `phangorn` is not
#'   installed, falls back to `trees[[1]]` with a warning.
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
#'     \item{`share_gnn`}{Logical; `TRUE` if the shared-GNN path was used.}
#'     \item{`fit`}{Single \code{pigauto_fit} trained on the reference
#'       tree when `share_gnn = TRUE`; `NULL` otherwise.}
#'     \item{`fits`}{List of `T` \code{pigauto_fit} objects (one per tree)
#'       when `share_gnn = FALSE`; `NULL` when `share_gnn = TRUE`.}
#'     \item{`reference_tree`}{The reference `phylo` used for GNN training
#'       when `share_gnn = TRUE`; `NULL` otherwise.}
#'     \item{`trees`}{The input posterior trees.}
#'     \item{`species_col`}{Passed-through species column name.}
#'   }
#'
#' @section Share-GNN (tree-sharing) mode:
#'
#' Under `share_gnn = TRUE` the GNN weights and spectral features are
#' trained once on the reference tree (MCC by default). For each
#' posterior tree the BM / joint-MVN baseline is recomputed, and the
#' prediction is the blend `(1 - r_cal) * baseline_t + r_cal * gnn_shared`.
#' Because `r_cal` is calibrated once on held-out data at the reference
#' tree and applied uniformly, the tree-uncertainty contribution is:
#'
#' * Fully preserved when the gate is closed (r_cal near 0): the GNN
#'   contributes nothing, and the baseline varies per tree.
#' * Partially preserved when the gate is open: the baseline portion
#'   still varies, but the GNN portion is a tree-invariant constant —
#'   this slightly under-estimates tree variance in the GNN channel.
#' * Lost in the GNN channel when the gate is fully open (rare on real
#'   data; the baseline channel still carries tree variation).
#'
#' On every real dataset benchmarked in the v0.9.0 campaign the gate
#' closed partially or fully, so `share_gnn = TRUE` is cheap AND honest.
#' Set `share_gnn = FALSE` if you need exact per-tree model independence.
#'
#' @details
#' For each tree the function runs the full pigauto pipeline
#' (preprocess -> baseline -> GNN -> predict) when `share_gnn = FALSE`.
#' With the default `share_gnn = TRUE`, the GNN is trained once and only
#' the baseline is recomputed per tree. Topologies and branch lengths vary
#' across trees, so the phylogenetic baseline covariance differs for each
#' tree.
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
#' uncertainty at the imputation step). The fraction of missing
#' information (FMI) reported by [pool_mi()] reflects both. To decompose
#' them, compare FMI from [multi_impute()] (single tree) with FMI from
#' `multi_impute_trees()`.
#'
#' **Computation time.** With `share_gnn = TRUE` (default): one GNN fit
#' + T cheap baseline passes. Rough budget on a modern CPU laptop:
#'
#' \tabular{rrrr}{
#'   Species n \tab 1 fit \tab T = 50 share_gnn=TRUE \tab T = 50 share_gnn=FALSE \cr
#'   300 \tab ~30-60 s \tab ~3-5 min \tab 25-50 min \cr
#'   5,000 \tab ~5-10 min \tab ~10-20 min \tab 4-8 hr \cr
#'   10,000 \tab ~20-40 min \tab ~30-60 min \tab 17-33 hr
#' }
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
#' @section Safety floor + share_gnn (v0.9.1.9002+):
#'   When \code{share_gnn = TRUE} with \code{safety_floor = TRUE}, the
#'   grand-mean baseline \code{mean_baseline_per_col} and the three
#'   calibrated weights (\code{r_cal_bm}, \code{r_cal_gnn},
#'   \code{r_cal_mean}) are computed ONCE on the reference tree and
#'   reused across all posterior trees.  They are properties of the
#'   observed training traits, not of the tree topology.  Each posterior
#'   tree only recomputes the BM baseline; the GNN delta and the three
#'   weights stay fixed.  This preserves the Nakagawa & de Villemereuil
#'   (2019) tree-uncertainty integration story without re-calibrating
#'   the safety floor per tree, and keeps the shared-GNN speedup intact.
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
#' # ---- Step 1: tree-aware imputation (canonical N&dV 2019 workflow) --
#' # 50 trees x 1 imputation = 50 completed datasets (fast with share_gnn=TRUE)
#' mi <- multi_impute_trees(df, trees300, m_per_tree = 1L)
#' print(mi)
#'
#' # ---- Step 2: tree-aware analysis (Nakagawa & de Villemereuil 2019)
#' # For each completed dataset, fit the downstream model using the SAME
#' # tree that produced that dataset. `mi$tree_index[i]` gives the tree
#' # index (1..T) for dataset `i`.
#' fits <- Map(
#'   function(dat, t_idx) {
#'     dat$species <- rownames(dat)
#'     nlme::gls(
#'       log(Mass) ~ log(Wing.Length),
#'       correlation = ape::corBrownian(
#'         phy = trees300[[t_idx]], form = ~species),
#'       data = dat, method = "ML"
#'     )
#'   },
#'   mi$datasets,
#'   mi$tree_index
#' )
#'
#' # Rubin's rules: pooled SEs include both trait-imputation and
#' # phylogenetic-tree uncertainty.
#' pool_mi(fits)
#' }
#'
#' @export
multi_impute_trees <- function(traits, trees, m_per_tree = 1L,
                               species_col = NULL,
                               trait_types = NULL,
                               multi_proportion_groups = NULL,
                               log_transform = TRUE,
                               missing_frac = 0.25,
                               covariates = NULL,
                               epochs = 2000L, verbose = TRUE,
                               seed = 1L,
                               share_gnn = TRUE,
                               reference_tree = NULL,
                               ...) {

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

  # Warn if total datasets M = T * m_per_tree is small (Rubin's rules wobble)
  if (T_trees * m_per_tree < 10L) {
    warning("Low total imputations (T * m_per_tree = ",
            T_trees * m_per_tree,
            "). Consider m_per_tree = 5L for stable Rubin's rules pooling ",
            "when T < 20.", call. = FALSE)
  }

  M_total <- T_trees * m_per_tree

  if (verbose) {
    cat(sprintf("multi_impute_trees: %d trees x %d imputations = %d datasets\n",
                T_trees, m_per_tree, M_total))
  }

  # ---- Dispatch on share_gnn ------------------------------------------------
  if (isTRUE(share_gnn)) {
    out <- run_shared_gnn(
      traits = traits, trees = trees, m_per_tree = m_per_tree,
      species_col = species_col, trait_types = trait_types,
      multi_proportion_groups = multi_proportion_groups,
      log_transform = log_transform, missing_frac = missing_frac,
      covariates = covariates, epochs = as.integer(epochs),
      verbose = verbose, seed = seed,
      reference_tree = resolve_reference_tree(trees, reference_tree),
      ...
    )
  } else {
    out <- run_per_tree(
      traits = traits, trees = trees, m_per_tree = m_per_tree,
      species_col = species_col, trait_types = trait_types,
      multi_proportion_groups = multi_proportion_groups,
      log_transform = log_transform, missing_frac = missing_frac,
      covariates = covariates, epochs = as.integer(epochs),
      verbose = verbose, seed = seed, ...
    )
  }

  out$share_gnn      <- isTRUE(share_gnn)
  out$reference_tree <- if (isTRUE(share_gnn)) resolve_reference_tree(trees, reference_tree) else NULL
  out$trees          <- trees
  out$species_col    <- species_col
  class(out) <- c("pigauto_mi_trees", "pigauto_mi")
  out
}

# Internal: per-tree loop (pre-v0.9.1 behaviour, opt-in via share_gnn=FALSE)
run_per_tree <- function(traits, trees, m_per_tree,
                         species_col, trait_types, multi_proportion_groups,
                         log_transform, missing_frac, covariates,
                         epochs, verbose, seed, ...) {
  T_trees      <- length(trees)
  M_total      <- T_trees * m_per_tree
  all_datasets <- vector("list", M_total)
  tree_index   <- integer(M_total)
  all_fits     <- vector("list", T_trees)
  imputed_mask <- NULL
  pooled_sum   <- NULL
  trait_cols   <- setdiff(names(traits), species_col)

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
      trait_types   = trait_types,
      multi_proportion_groups = multi_proportion_groups,
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

    # Build completed data.frames.
    # When n_imputations=1 predict() returns pred$imputed directly (not a list).
    # When n_imputations>1 it returns pred$imputed_datasets (a list of m).
    imp_list <- if (!is.null(pred$imputed_datasets)) {
      pred$imputed_datasets
    } else {
      list(pred$imputed)
    }
    if (length(imp_list) != m_per_tree) {
      stop("Tree ", t, ": predict() did not return ", m_per_tree,
           " imputed datasets. This is an internal error.", call. = FALSE)
    }

    # `res$data$input_row_order` re-aligns internal-order imputations back
    # to the user's input row order (multi-obs reordering fix, 2026-04-26).
    # Each tree gets its own pigauto_data object, so we read this fresh
    # per tree.
    input_row_order <- res$data$input_row_order
    for (k in seq_len(m_per_tree)) {
      idx <- (t - 1L) * m_per_tree + k
      completed_info <- build_completed(traits, imp_list[[k]], species_col,
                                          input_row_order = input_row_order)
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
    if (is.null(pooled_sum)) {
      pooled_sum <- res$completed
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

  # Pooled point estimate (average across trees)
  pooled_point <- pooled_sum
  for (nm in trait_cols) {
    if (is.numeric(pooled_point[[nm]])) {
      pooled_point[[nm]] <- pooled_point[[nm]] / T_trees
    }
  }

  list(
    datasets     = all_datasets,
    m            = M_total,
    n_trees      = T_trees,
    m_per_tree   = m_per_tree,
    tree_index   = tree_index,
    pooled_point = pooled_point,
    se           = NULL,
    imputed_mask = imputed_mask,
    fits         = all_fits,
    fit          = NULL
  )
}

# Internal: shared-GNN path — fit once on reference_tree, recompute only
# the BM baseline per posterior tree, reuse the trained GNN.
run_shared_gnn <- function(traits, trees, m_per_tree,
                           species_col, trait_types, multi_proportion_groups,
                           log_transform, missing_frac, covariates,
                           epochs, verbose, seed, reference_tree, ...) {
  T_trees <- length(trees)
  M_total <- T_trees * m_per_tree
  if (verbose) {
    cat(sprintf("multi_impute_trees (share_gnn=TRUE): %d trees x %d imputations = %d datasets\n",
                T_trees, m_per_tree, M_total))
    cat("  Training GNN once on reference tree...\n")
  }

  # Fit on the reference tree — this is the expensive step, done ONCE.
  res_ref <- impute(
    traits = traits, tree = reference_tree,
    species_col = species_col, trait_types = trait_types,
    multi_proportion_groups = multi_proportion_groups,
    log_transform = log_transform, missing_frac = missing_frac,
    n_imputations = m_per_tree, covariates = covariates,
    epochs = as.integer(epochs), verbose = FALSE, seed = as.integer(seed), ...
  )
  fit_ref    <- res_ref$fit
  data_ref   <- res_ref$data
  splits_ref <- res_ref$splits
  graph_ref  <- fit_ref$graph   # graph is stored on the fit, not on pigauto_result

  all_datasets <- vector("list", M_total)
  tree_index   <- integer(M_total)
  imputed_mask <- NULL
  se_sum       <- NULL
  pooled_sum   <- NULL
  idx <- 0L

  for (t in seq_len(T_trees)) {
    if (verbose) cat(sprintf("  Tree %d/%d: baseline only...", t, T_trees))
    t0 <- proc.time()

    baseline_t <- fit_baseline(data_ref, trees[[t]], splits = splits_ref,
                               graph = graph_ref)
    pred_t <- stats::predict(fit_ref, return_se = TRUE,
                              n_imputations = m_per_tree,
                              baseline_override = baseline_t)

    # Build completed data.frames (reuse build_completed).  In the share-GNN
    # path, data_ref / fit_ref are tied to the FIRST tree's preprocess pass;
    # since input row order is independent of the tree (it depends only on
    # the input data.frame's species column), we can reuse data_ref's
    # input_row_order across all trees.  (Multi-obs reordering fix, 2026-04-26.)
    input_row_order <- data_ref$input_row_order
    for (k in seq_len(m_per_tree)) {
      idx <- idx + 1L
      one <- if (!is.null(pred_t$imputed_datasets)) pred_t$imputed_datasets[[k]]
             else pred_t$imputed
      info <- build_completed(traits, one, species_col,
                                input_row_order = input_row_order)
      all_datasets[[idx]] <- info$completed
      tree_index[idx]     <- t
      if (is.null(imputed_mask)) imputed_mask <- info$imputed_mask
    }

    # Running sums for pooled point / se
    if (!is.null(pred_t$imputed)) {
      m_imp <- as.matrix(pred_t$imputed[, vapply(pred_t$imputed, is.numeric, logical(1)), drop = FALSE])
      pooled_sum <- if (is.null(pooled_sum)) m_imp else pooled_sum + m_imp
    }
    if (!is.null(pred_t$se)) {
      m_se <- as.matrix(pred_t$se)
      se_sum <- if (is.null(se_sum)) m_se else se_sum + m_se
    }

    if (verbose) {
      elapsed <- (proc.time() - t0)[["elapsed"]]
      cat(sprintf(" done (%.1fs)\n", elapsed))
    }
  }

  pooled_point <- if (!is.null(pooled_sum)) as.data.frame(pooled_sum / T_trees) else NULL
  pooled_se    <- if (!is.null(se_sum))    se_sum / T_trees                   else NULL

  list(
    datasets     = all_datasets,
    m            = M_total,
    n_trees      = T_trees,
    m_per_tree   = m_per_tree,
    tree_index   = tree_index,
    pooled_point = pooled_point,
    se           = pooled_se,
    imputed_mask = imputed_mask,
    fit          = fit_ref,
    fits         = NULL
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
