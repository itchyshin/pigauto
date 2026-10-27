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

#' Posterior-tree prediction sensitivity
#'
#' Run pigauto's full imputation pipeline on each of `T` posterior
#' phylogenies, generating `m_per_tree` datasets per tree, either as point
#' predictions or stochastic completions. This returns `T * m_per_tree`
#' datasets. Each completed
#' dataset is conditional on a specific posterior tree
#' (recorded in `mi$tree_index`).
#'
#' **This is an experimental prediction-sensitivity check.** Use it to see
#' whether point predictions change across the trees you supply. It does not
#' provide inferential uncertainty for downstream model estimates. See
#' `vignette("tree-uncertainty")` for a step-by-step example.
#'
#' Each dataset reflects the tree assigned to it. Comparing datasets across
#' trees shows prediction sensitivity to the supplied tree sample. When
#' multiple draws are requested per tree, they can also vary within a tree,
#' depending on the selected draw method. These results are diagnostic, not
#' analysis-ready imputations. The downstream
#' fitting and pooling helpers reject them, so do not use them for inference.
#' This workflow does not combine tree-to-tree variation with imputation
#' variation into an inferential estimate.
#'
#' @section When to use this:
#'
#' Choose the prediction-diagnostic function based on how many trees you have:
#'
#' * **One tree** (single published phylogeny, single time-calibrated tree):
#'   use [multi_impute()] and select a draw method explicitly when needed.
#' * **Multiple posterior trees** (BirdTree samples, BEAST posterior, etc.):
#'   use [multi_impute_trees()] only to examine prediction sensitivity. Its
#'   results cannot be pooled for downstream inference.
#'
#' @importFrom utils head
#'
#' @param traits data.frame. Same format as [multi_impute()] and [impute()].
#' @param trees list of `phylo` objects (class `multiPhylo` or plain list).
#'   Each tree must contain the species in `traits` as tips. Posterior
#'   samples from BirdTree.org (Jetz et al. 2012) are ideal; the bundled
#'   [trees300] dataset provides 50 posterior trees for [avonet300].
#' @param m_per_tree integer. Number of completed datasets returned per
#'   tree (default `1`). Total datasets = `length(trees) * m_per_tree`.
#'   With `draws_method = "mc_dropout"` and the default of one, the
#'   prediction is deterministic conditional on each fitted model; this
#'   setting does not estimate within-tree variation.
#' @param species_col character or `NULL`. See [impute()].
#' @param trait_types named character vector overriding auto-detected
#'   trait types. Required for `"proportion"` and `"zi_count"`. See
#'   [impute()] / [preprocess_traits()]. Default `NULL` (auto-detect).
#' @param multi_proportion_groups named list declaring compositional
#'   trait groups (rows summing to 1), forwarded to [impute()] /
#'   [preprocess_traits()]. Default `NULL`.
#' @param log_transform logical. Log-transform every continuous trait whose
#'   observed values are all positive
#'   (default `TRUE`), then back-transform the imputations. Set it to
#'   `FALSE` when traits are already on a log scale (for example log body
#'   mass): otherwise they are logged a second time and imputed on a
#'   log(log) scale, which can bias analyses done on the log scale. With a
#'   mix of raw and already-logged traits, log the raw ones yourself and set
#'   `FALSE`. In a real-data check (PanTHERIA, 4,027 mammal species, four
#'   already-logged traits, posterior draws), the default shifted
#'   Rubin-pooled regression slopes by up to 3.9 standard errors of the
#'   complete-data slope, in a direction that depended on which values were
#'   missing; with `FALSE` the largest shift was 1.8.
#' @param missing_frac numeric. Fraction held out for validation/test
#'   during training (default `0.25`).
#' @param covariates data.frame or matrix of environmental covariates
#'   (fully observed, numeric). Passed through to [impute()].
#'   Default `NULL` (no covariates).
#' @param epochs integer. Maximum GNN training epochs per tree (default
#'   `2000`).
#' @param verbose logical. Print progress (default `TRUE`).
#' @param seed optional integer. Base random seed. How it is used depends on
#'   `share_gnn` and `draws_method`; see the "Reproducibility" section.
#'   The default `NULL` uses the current RNG stream.
#' @param share_gnn logical. With `gnn = TRUE`, `TRUE` (the default) trains
#'   one GNN on a reference tree and reuses it while recomputing the
#'   baseline for each posterior tree. Set `FALSE` to refit the GNN on each
#'   tree. With the default `gnn = FALSE`, no GNN is trained in either mode.
#'   See the "Share-GNN" section below.
#' @param reference_tree optional `phylo` used for the reference fit when
#'   `share_gnn = TRUE`. Default `NULL` selects the maximum-clade-credibility
#'   tree via `phangorn::maxCladeCred(trees)`. If `phangorn` is not
#'   installed, falls back to `trees[[1]]` with a warning.
#' @param draws_method character. How missing cells are completed. One of
#'   `"mc_dropout"` (default) or `"conformal"`. The methods differ in
#'   whether they add stochastic variation; see the section below. A single
#'   `mc_dropout` completion is a deterministic prediction, while
#'   `"conformal"` runs its sampling step even when `m_per_tree = 1`.
#' @param gnn logical. Passed through to [impute()] / [fit_pigauto()] for
#'   every tree. When `TRUE`, the usual GNN correction is trained
#'   (once, on the reference tree, when `share_gnn = TRUE`). When `FALSE`,
#'   no GNN is used anywhere in the pipeline: baseline-only fits
#'   throughout (the default since 0.11.0.9001; see [impute()]). With
#'   `share_gnn = TRUE`, each tree's baseline is fit without held-out
#'   splits (`splits = NULL`), matching the production
#'   `baseline_full` semantics [fit_pigauto()] uses.
#' @param ... additional arguments forwarded to [fit_pigauto()] via
#'   [impute()].
#'
#' @return An object of class `"pigauto_mi_trees"`, inheriting from
#'   `"pigauto_mi"`, with components:
#'   \describe{
#'     \item{`datasets`}{List of `T * m_per_tree` data.frames. Observed
#'       cells are preserved. Missing cells are filled with tree-conditional
#'       predictions or draws where available. A cell marked for imputation
#'       can remain `NA` if no prediction is available. These datasets are
#'       for sensitivity diagnostics.}
#'     \item{`m`}{Total number of datasets (`T * m_per_tree`).}
#'     \item{`n_trees`}{Number of posterior trees used.}
#'     \item{`m_per_tree`}{Number of prediction or completion datasets per tree.}
#'     \item{`tree_index`}{Integer vector of length `m`; element `i` gives
#'       the tree index (1..T) for dataset `i`.}
#'     \item{`pooled_point`}{Single data.frame averaging numeric trait columns
#'       across all `T * m_per_tree` datasets; an `NA` in a numeric trait cell
#'       propagates to its average. Factor and ordered-factor columns are
#'       taken from the first dataset. Compare the
#'       `datasets` directly to describe nonnumeric prediction sensitivity.
#'       For reporting, not inference.}
#'     \item{`se`}{On the shared-GNN path, the arithmetic average of the
#'       per-tree prediction uncertainty summaries (see [predict.pigauto_fit()]).
#'       Descriptive reporting only, not a pooled inferential SE. `NULL`
#'       when unavailable, including the separate-fit path.}
#'     \item{`imputed_mask`}{Logical matrix aligned to the input rows and
#'       trait columns; `TRUE` marks an originally missing target in a row
#'       included in the fit, not necessarily a successful fill. A marked
#'       cell can remain `NA` if no prediction was available. When
#'       `species_col` is supplied, its column is excluded.}
#'     \item{`share_gnn`}{Logical; `TRUE` if the shared-GNN path was used.}
#'     \item{`draws_method`}{Character; `"mc_dropout"` or `"conformal"`,
#'       echoing the argument this call used.}
#'     \item{`mi_workflow`}{`"pigauto_tree_sensitivity_diagnostic"`,
#'       recording that these are prediction-sensitivity completions and
#'       cannot be passed to [with_imputations()] or [pool_mi()].}
#'     \item{`fit`}{Single \code{pigauto_fit} fitted on the reference
#'       tree when `share_gnn = TRUE`; `NULL` otherwise.}
#'     \item{`fits`}{List of `T` \code{pigauto_fit} objects (one per tree)
#'       when `share_gnn = FALSE`; `NULL` when `share_gnn = TRUE`.}
#'     \item{`reference_tree`}{The reference `phylo` used for the shared fit
#'       when `share_gnn = TRUE`; `NULL` otherwise.}
#'     \item{`trees`}{The input posterior trees.}
#'     \item{`species_col`}{Passed-through species column name.}
#'     \item{`lambda_per_trait_by_tree`}{Returned only when
#'       `share_gnn = TRUE`. A list indexed in the same order as `trees`;
#'       element `i` is the named `lambda_per_trait` vector returned by
#'       [fit_baseline()] for `trees[[i]]`. Values can be fixed rather than
#'       estimated. Categorical one-vs-rest lambda estimates are not included
#'       in this vector, so a categorical column's default value of 1 is not
#'       an estimated class-specific lambda.}
#'   }
#'
#' @section Reproducibility:
#'
#' Pigauto routes the base seed separately to tree fits and sampling steps,
#' so the offset depends on the fitting mode. Here, \eqn{t} is a tree's
#' position in the supplied list. With `share_gnn = FALSE`, the separate fit
#' for tree \eqn{t} receives
#' `seed + t - 1`. With `share_gnn = TRUE`, the reference-tree fit receives
#' `seed`, and each tree's baseline refit receives `seed + t`. Stochastic
#' predictions reuse the stored seed from the reference fit; deterministic
#' single-pass predictions do not set a new seed. If
#' `draws_method = "conformal"`, the tree-specific completion sampling
#' step receives `seed + t - 1`; each completion within a tree also receives
#' its own seed offset. With `seed = NULL`, these steps use the current random-number stream.
#'
#' @section Share-GNN (tree-sharing) mode:
#'
#' When `gnn = TRUE, share_gnn = TRUE`, neural weights, the reference
#' graph, spectral features, and calibrated blend weights are reused from
#' the reference fit (MCC by default). The baseline is recomputed for each
#' posterior tree and passed into prediction again. Because baseline values
#' also enter the neural inputs, shared weights do not imply a constant neural
#' prediction across trees. No direction or magnitude of tree-variance bias
#' has been established for this approximation.
#'
#' Set `share_gnn = FALSE` with `gnn = TRUE` to fit a separate GNN for each
#' tree. With the default `gnn = FALSE`, there is no GNN channel.
#'
#' @details
#' For each tree the function runs preprocessing, baseline fitting and
#' prediction. A GNN is also fitted only when `gnn = TRUE`; with
#' `share_gnn = TRUE`, that GNN is trained once on a reference tree.
#' Topologies and branch lengths vary across trees, so the phylogenetic
#' baseline covariance differs for each
#' tree.
#'
#' The returned datasets may be compared descriptively to assess sensitivity
#' of point imputations to the tree sample. No calibrated downstream standard
#' error or Rubin-pooling claim is made for this path.
#'
#' **Computation time.** The following historical planning estimates apply
#' to `gnn = TRUE`, with one shared GNN fit plus T baseline passes when
#' `share_gnn = TRUE`. They are not timings for the current GNN-off default;
#' estimate the intended tree and sample setting before a campaign.
#'
#' \ifelse{html}{\out{
#' <table>
#' <thead><tr><th scope="col">Species n</th><th scope="col">1 fit</th><th scope="col">T = 50 share_gnn=TRUE</th><th scope="col">T = 50 share_gnn=FALSE</th></tr></thead>
#' <tbody>
#' <tr><td>300</td><td>~30-60 s</td><td>~3-5 min</td><td>25-50 min</td></tr>
#' <tr><td>5,000</td><td>~5-10 min</td><td>~10-20 min</td><td>4-8 hr</td></tr>
#' <tr><td>10,000</td><td>~20-40 min</td><td>~30-60 min</td><td>17-33 hr</td></tr>
#' </tbody>
#' </table>
#' }}{
#' \tabular{rrrr}{
#'   Species n \tab 1 fit \tab T = 50 share_gnn=TRUE \tab T = 50 share_gnn=FALSE \cr
#'   300 \tab ~30-60 s \tab ~3-5 min \tab 25-50 min \cr
#'   5,000 \tab ~5-10 min \tab ~10-20 min \tab 4-8 hr \cr
#'   10,000 \tab ~20-40 min \tab ~30-60 min \tab 17-33 hr
#' }
#' }
#'
#' @section How the completion methods vary (and their limitations):
#' With `draws_method = "mc_dropout"` (the default), one completion per
#' tree is a deterministic prediction conditional on that tree's fitted
#' model. When `m_per_tree > 1`, pigauto adds independent Gaussian
#' perturbations to latent cells, scaled by their baseline standard
#' deviations. If `gnn = TRUE`, dropout can add further variation. These
#' perturbations are not joint Brownian-motion posterior draws. Variation
#' can be zero when a source has no uncertainty or its contribution is
#' weighted out; decoding can also map different latent values to the same
#' category.
#'
#' With `draws_method = "conformal"`, missing-value sampling occurs even
#' when `m_per_tree = 1`. For continuous, count, ordinal and proportion
#' traits, it uses conformal widths where available and reported uncertainty.
#' When no finite conformal score is available for these traits, the route uses the reported
#' `pred$se` as the fallback scale input. For log-transformed
#' continuous traits, log1p-transformed counts, and logit-transformed
#' proportions, `pred$se` is reported on the original scale, but the draw is
#' made on the transformed scale. Pigauto does not fully convert `pred$se`
#' between these scales, so diagnostic draw spread may be too wide or too
#' narrow. Fitted point predictions are unchanged. These tree-based conformal
#' completions use global trait scores; they do not carry
#' per-cell near/far Mondrian scores into the draws, even when prediction
#' intervals use Mondrian calibration. Binary
#' and categorical traits are sampled from predicted probabilities;
#' zero-inflated counts combine a Bernoulli gate with a magnitude draw;
#' compositions use latent standard deviations. A fixed seed makes the
#' draws reproducible. Neither method validates downstream inference.
#'
#' When `share_gnn = TRUE`, conformal scores and calibrated blend weights
#' come from the reference fit and are reused across trees. With
#' `share_gnn = FALSE`, each tree uses scores and weights from its own fit.
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
#' @section Safety floor with a shared fit:
#'   When \code{share_gnn = TRUE} with \code{safety_floor = TRUE}, the
#'   grand-mean baseline \code{mean_baseline_per_col} and the three
#'   calibrated weights (\code{r_cal_bm}, \code{r_cal_gnn},
#'   \code{r_cal_mean}) are computed ONCE on the reference tree and
#'   reused across all posterior trees.  They are properties of the
#'   reference fit and are not recalibrated for each topology. Each posterior
#'   tree recomputes the baseline. The three weights stay fixed, while the
#'   GNN delta can change with its per-tree baseline inputs. This preserves
#'   prediction sensitivity to the
#'   per-tree baseline without re-calibrating the safety floor; it does not
#'   establish downstream inferential validity.
#'
#' @seealso [multi_impute()] for single-tree prediction diagnostics or,
#'   with `draws_method = "posterior"`, single-tree inference on continuous
#'   traits; [multi_impute_analysis()] for analysis-aware MI; and [trees300]
#'
#' @examples
#' \donttest{
#' data(avonet300, tree300)
#' tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
#' trees <- structure(list(tree, tree), class = "multiPhylo")
#' df <- avonet300[match(tree$tip.label, avonet300$Species_Key),
#'                 c("Mass", "Wing.Length"), drop = FALSE]
#' rownames(df) <- tree$tip.label
#' df$Mass[seq_len(3L)] <- NA_real_
#' mi <- multi_impute_trees(df, trees, m_per_tree = 1L, epochs = 5L,
#'                          verbose = FALSE)
#' print(mi)
#' mass_by_tree <- vapply(mi$datasets, function(dat) dat$Mass,
#'                        numeric(nrow(df)))
#' apply(mass_by_tree, 1L, stats::sd)
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
                               seed = NULL,
                               share_gnn = TRUE,
                               reference_tree = NULL,
                               draws_method = c("mc_dropout", "conformal"),
                               gnn = FALSE,
                               ...) {

  draws_method <- match.arg(draws_method)

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

  # Warn when a tree-sensitivity summary is based on very few datasets.
  if (T_trees * m_per_tree < 10L) {
    warning("Few tree-sensitivity datasets (T * m_per_tree = ",
            T_trees * m_per_tree,
            "). Descriptive sensitivity summaries may be unstable.",
            call. = FALSE)
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
      draws_method = draws_method,
      gnn = gnn,
      ...
    )
  } else {
    out <- run_per_tree(
      traits = traits, trees = trees, m_per_tree = m_per_tree,
      species_col = species_col, trait_types = trait_types,
      multi_proportion_groups = multi_proportion_groups,
      log_transform = log_transform, missing_frac = missing_frac,
      covariates = covariates, epochs = as.integer(epochs),
      verbose = verbose, seed = seed, draws_method = draws_method,
      gnn = gnn, ...
    )
  }

  out$share_gnn      <- isTRUE(share_gnn)
  out$reference_tree <- if (isTRUE(share_gnn)) resolve_reference_tree(trees, reference_tree) else NULL
  out$trees          <- trees
  out$species_col    <- species_col
  out$mi_workflow    <- "pigauto_tree_sensitivity_diagnostic"
  class(out) <- c("pigauto_tree_sensitivity_diagnostic", "pigauto_mi_trees",
                  "pigauto_mi")
  out
}

# Internal: per-tree loop (pre-v0.9.1 behaviour, opt-in via share_gnn=FALSE)
run_per_tree <- function(traits, trees, m_per_tree,
                         species_col, trait_types, multi_proportion_groups,
                         log_transform, missing_frac, covariates,
                         epochs, verbose, seed, draws_method = "mc_dropout",
                         gnn = FALSE, ...) {
  T_trees      <- length(trees)
  M_total      <- T_trees * m_per_tree
  all_datasets <- vector("list", M_total)
  tree_index   <- integer(M_total)
  all_fits     <- vector("list", T_trees)
  imputed_mask <- NULL
  pooled_sum   <- NULL
  pooled_n     <- 0L
  trait_cols   <- setdiff(names(traits), species_col)
  use_conformal <- identical(draws_method, "conformal")

  for (t in seq_len(T_trees)) {
    t_seed <- if (is.null(seed)) NULL else as.integer(seed + t - 1L)

    if (verbose) {
      cat(sprintf("  Tree %d/%d%s...", t, T_trees,
                  if (is.null(t_seed)) "" else paste0(" (seed=", t_seed, ")")))
      t_start <- proc.time()
    }

    # Run full pipeline. draws_method = "conformal" fits a single
    # deterministic pass (n_imputations = 1) and draws the m_per_tree
    # completions from the calibrated conformal scores via .conformal_draws()
    # -- the same sampling helper multi_impute() uses. draws_method =
    # "mc_dropout" (default) asks predict() for m_per_tree passes. A single
    # pass is a deterministic point prediction conditional on the fitted
    # state; multiple passes add baseline perturbations, with GNN dropout
    # contributing further variation when the GNN is enabled.
    res <- impute(
      traits        = traits,
      tree          = trees[[t]],
      species_col   = species_col,
      trait_types   = trait_types,
      multi_proportion_groups = multi_proportion_groups,
      log_transform = log_transform,
      missing_frac  = missing_frac,
      n_imputations = if (use_conformal) 1L else m_per_tree,
      covariates    = covariates,
      epochs        = as.integer(epochs),
      verbose       = FALSE,
      seed          = t_seed,
      gnn           = gnn,
      ...
    )

    pred <- res$prediction

    # `res$data$input_row_order` re-aligns internal-order imputations back
    # to the user's input row order (multi-obs reordering fix, 2026-04-26).
    # Each tree gets its own pigauto_data object, so we read this fresh
    # per tree.
    input_row_order <- res$data$input_row_order

    if (use_conformal) {
      completed_list <- .conformal_draws(traits, pred, res$imputed_mask,
                                         res$fit$trait_map, m_per_tree,
                                         species_col = species_col,
                                         seed = t_seed,
                                         input_row_order = input_row_order)
    } else {
      # When n_imputations=1 predict() returns pred$imputed directly (not a
      # list). When n_imputations>1 it returns pred$imputed_datasets (a list
      # of m).
      imp_list <- if (!is.null(pred$imputed_datasets)) {
        pred$imputed_datasets
      } else {
        list(pred$imputed)
      }
      if (length(imp_list) != m_per_tree) {
        stop("Tree ", t, ": predict() did not return ", m_per_tree,
             " imputed datasets. This is an internal error.", call. = FALSE)
      }
      completed_list <- lapply(imp_list, function(imp_df) {
        build_completed(traits, imp_df, species_col,
                        input_row_order = input_row_order)$completed
      })
    }

    for (k in seq_len(m_per_tree)) {
      idx <- (t - 1L) * m_per_tree + k
      completed_i <- completed_list[[k]]
      all_datasets[[idx]] <- completed_i
      tree_index[idx]     <- t

      if (is.null(imputed_mask)) {
        imputed_mask <- res$imputed_mask
      }

      if (is.null(pooled_sum)) {
        pooled_sum <- completed_i
        for (nm in trait_cols) {
          if (is.numeric(completed_i[[nm]])) {
            pooled_sum[[nm]] <- as.numeric(completed_i[[nm]])
          }
        }
      } else {
        for (nm in trait_cols) {
          if (is.numeric(completed_i[[nm]])) {
            pooled_sum[[nm]] <- pooled_sum[[nm]] + as.numeric(completed_i[[nm]])
          }
        }
      }
      pooled_n <- pooled_n + 1L
    }

    # Store the fit for diagnostics
    all_fits[[t]] <- res$fit

    if (verbose) {
      elapsed <- (proc.time() - t_start)[3]
      cat(sprintf(" %.1fs\n", elapsed))
    }
  }

  pooled_point <- pooled_sum
  if (!is.null(pooled_point) && pooled_n > 0L) {
    for (nm in trait_cols) {
      if (is.numeric(pooled_point[[nm]])) {
        pooled_point[[nm]] <- pooled_point[[nm]] / pooled_n
      }
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
    fit          = NULL,
    draws_method = draws_method
  )
}

# Internal: shared-GNN path — fit once on reference_tree, recompute only
# the BM baseline per posterior tree, reuse the trained GNN.
run_shared_gnn <- function(traits, trees, m_per_tree,
                           species_col, trait_types, multi_proportion_groups,
                           log_transform, missing_frac, covariates,
                           epochs, verbose, seed, reference_tree,
                           draws_method = "mc_dropout", gnn = FALSE, ...) {
  T_trees <- length(trees)
  M_total <- T_trees * m_per_tree
  dots <- list(...)
  use_conformal <- identical(draws_method, "conformal")
  if (verbose) {
    cat(sprintf("multi_impute_trees (share_gnn=TRUE): %d trees x %d imputations = %d datasets\n",
                T_trees, m_per_tree, M_total))
    if (isTRUE(gnn)) {
      cat("  Training GNN once on reference tree...\n")
    } else {
      cat("  GNN disabled; fitting the reference-tree baseline only...\n")
    }
  }

  # Fit on the reference tree — this is the expensive step, done ONCE.
  res_ref <- impute(
    traits = traits, tree = reference_tree,
    species_col = species_col, trait_types = trait_types,
    multi_proportion_groups = multi_proportion_groups,
    log_transform = log_transform, missing_frac = missing_frac,
    n_imputations = m_per_tree, covariates = covariates,
    epochs = as.integer(epochs), verbose = FALSE,
    seed = if (is.null(seed)) NULL else as.integer(seed),
    gnn = gnn, ...
  )
  fit_ref    <- res_ref$fit
  data_ref   <- res_ref$data
  splits_ref <- res_ref$splits

  all_datasets <- vector("list", M_total)
  tree_index   <- integer(M_total)
  imputed_mask <- NULL
  se_sum       <- NULL
  pooled_sum   <- NULL
  pooled_n     <- 0L
  trait_cols   <- setdiff(names(traits), species_col)
  idx <- 0L
  # Per-tree Pagel's lambda from each tree's fit_baseline() (NULL only for
  # baselines built before the lambda fields existed). One entry per tree (length T_trees, not M_total).
  lambda_per_trait_by_tree <- vector("list", T_trees)

  for (t in seq_len(T_trees)) {
    if (verbose) cat(sprintf("  Tree %d/%d: baseline only...", t, T_trees))
    t0 <- proc.time()

    # Capture the resolved baseline config from the reference fit
    config <- res_ref$fit$model_config
    baseline_arg <- function(name, default) {
      val <- dots[[name]]
      if (is.null(val) || length(val) == 0L) val <- config[[name]]
      if (is.null(val) || length(val) == 0L) default else val
    }

    # Do not pass graph_ref here. The shared GNN keeps the reference-tree
    # graph inside fit_ref for prediction, but each tree-specific baseline
    # must recompute its own cophenetic distances and phylogenetic covariance.
    #
    # gnn = FALSE (S3, docs/dev-log/arc/2026-09-18-gnn-off-contract.md):
    # fit each tree's baseline tax-free (splits = NULL, every observed cell
    # informs every other cell), matching the production `baseline_full`
    # semantics fit_pigauto()/impute() use elsewhere. gnn = TRUE keeps the
    # held-out splits_ref fit, unchanged.
    # S5b: for gnn = FALSE (splits = NULL below, production-baseline
    # semantics) there are no validation cells for this per-tree refit to
    # re-decide "auto" with, so reuse the reference tree's own per-trait
    # route (fit_ref$model_config$predict_method_by_trait, from its own
    # held-out `baseline`) via predict_route. For gnn = TRUE (splits_ref
    # below, real validation cells), let "auto" re-decide per tree as
    # usual -- no predict_route override needed.
    predict_route_arg <- if (isFALSE(gnn)) {
      baseline_arg("predict_method_by_trait", NULL)
    } else {
      NULL
    }
    baseline_t <- fit_baseline(
      data_ref,
      trees[[t]],
      splits = if (isFALSE(gnn)) NULL else splits_ref,
      graph = NULL,
      # Replay all the config options (joint_solver / predict_method /
      # joint_refine_iter added here -- pre-existing omission, now
      # load-bearing since gnn = FALSE fits go through the joint /
      # threshold-joint baseline dispatch more often than the GNN-on
      # default):
      lambda_mode = baseline_arg("lambda_mode", "estimate"),
      # Fits recorded before discrete_lambda existed used lambda = 1.
      discrete_lambda = baseline_arg("discrete_lambda", "fixed_1"),
      multi_obs_aggregation = baseline_arg("multi_obs_aggregation", "hard"),
      em_iterations = baseline_arg("em_iterations", 0L),
      em_tol = baseline_arg("em_tol", 1e-3),
      em_offdiag = baseline_arg("em_offdiag", FALSE),
      joint_solver = baseline_arg("joint_solver", "inhouse"),
      predict_method = baseline_arg("predict_method", "auto"),
      joint_refine_iter = baseline_arg("joint_refine_iter", 0L),
      predict_route = predict_route_arg,
      seed = if (is.null(seed)) NULL else as.integer(seed) + t
    )
    # NULL only for baselines built before lambda_per_trait existed.
    lambda_per_trait_by_tree[[t]] <- baseline_t$lambda_per_trait
    # draws_method = "conformal": a single deterministic pass per
    # tree, then draw m_per_tree completions from the conformal scores held
    # on fit_ref (calibrated once, on the reference tree) via the same
    # sampling helper multi_impute() uses. draws_method = "mc_dropout"
    # (default): ask predict() for m_per_tree passes. A single pass is a
    # deterministic point prediction conditional on the fitted state;
    # multiple passes add baseline perturbations, with GNN dropout
    # contributing further variation when the GNN is enabled.
    pred_t <- stats::predict(fit_ref, return_se = TRUE,
                              n_imputations = if (use_conformal) 1L else m_per_tree,
                              baseline_override = baseline_t)

    # Build completed data.frames (reuse build_completed).  In the share-GNN
    # path, data_ref / fit_ref are tied to the FIRST tree's preprocess pass;
    # since input row order is independent of the tree (it depends only on
    # the input data.frame's species column), we can reuse data_ref's
    # input_row_order across all trees.  (Multi-obs reordering fix, 2026-04-26.)
    input_row_order <- data_ref$input_row_order

    if (use_conformal) {
      # res_ref$imputed_mask does not depend on the tree (only on which
      # cells are NA in `traits`), so it is safe to reuse across trees --
      # matching the imputed_mask capture below, which is likewise set once.
      completed_list <- .conformal_draws(traits, pred_t, res_ref$imputed_mask,
                                         fit_ref$trait_map, m_per_tree,
                                         species_col = species_col,
                                         seed = if (is.null(seed)) NULL else as.integer(seed) + t - 1L,
                                         input_row_order = input_row_order)
    } else {
      completed_list <- lapply(seq_len(m_per_tree), function(k) {
        one <- if (!is.null(pred_t$imputed_datasets)) pred_t$imputed_datasets[[k]]
               else pred_t$imputed
        build_completed(traits, one, species_col,
                        input_row_order = input_row_order)$completed
      })
    }

    for (k in seq_len(m_per_tree)) {
      idx <- idx + 1L
      completed_i <- completed_list[[k]]
      all_datasets[[idx]] <- completed_i
      tree_index[idx]     <- t
      if (is.null(imputed_mask)) imputed_mask <- res_ref$imputed_mask

      # Pool from the user-facing completed datasets, not raw pred_t$imputed:
      # completed_i is already in input row order, preserves species_col and
      # observed cells, and represents each of the m_per_tree draws.
      if (is.null(pooled_sum)) {
        pooled_sum <- completed_i
        for (nm in trait_cols) {
          if (is.numeric(completed_i[[nm]])) {
            pooled_sum[[nm]] <- as.numeric(completed_i[[nm]])
          }
        }
      } else {
        for (nm in trait_cols) {
          if (is.numeric(completed_i[[nm]])) {
            pooled_sum[[nm]] <- pooled_sum[[nm]] + as.numeric(completed_i[[nm]])
          }
        }
      }
      pooled_n <- pooled_n + 1L
    }

    # Running sum for SE. pred_t$se is in prediction/internal row order, so
    # align it to the user's input rows before returning it beside pooled_point.
    if (!is.null(pred_t$se)) {
      m_se <- as.matrix(pred_t$se)
      if (!is.null(input_row_order)) {
        row_idx <- match(seq_len(nrow(traits)), input_row_order)
        m_se <- m_se[row_idx, , drop = FALSE]
        rownames(m_se) <- rownames(traits)
      }
      se_sum <- if (is.null(se_sum)) m_se else se_sum + m_se
    }

    if (verbose) {
      elapsed <- (proc.time() - t0)[["elapsed"]]
      cat(sprintf(" done (%.1fs)\n", elapsed))
    }
  }

  pooled_point <- pooled_sum
  if (!is.null(pooled_point) && pooled_n > 0L) {
    for (nm in trait_cols) {
      if (is.numeric(pooled_point[[nm]])) {
        pooled_point[[nm]] <- pooled_point[[nm]] / pooled_n
      }
    }
  }
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
    fits         = NULL,
    draws_method = draws_method,
    # List of length T_trees (indexed like tree_index's distinct values),
    # each element the per-tree baseline's lambda_per_trait. NULL entries
    # only for baselines built before this field existed.
    lambda_per_trait_by_tree = lambda_per_trait_by_tree
  )
}


#' @rdname multi_impute_trees
#' @param x A posterior-tree prediction-sensitivity object.
#' @section Printing:
#' Prints tree, row, trait, and targeted-cell counts and the diagnostic-only
#' scope. A targeted cell can remain missing when the prediction is unavailable.
#' The draw-method label does not imply dropout ran with GNN off or with one
#' completion per tree. The `...` argument is currently unused.
#' @return The print method returns `x` invisibly.
#' @export
print.pigauto_mi_trees <- function(x, ...) {
  T_trees    <- x$n_trees
  m_per_tree <- x$m_per_tree
  M_total    <- x$m

  # Count rows and traits from the first dataset
  d1 <- x$datasets[[1]]
  trait_cols <- setdiff(names(d1), x$species_col)
  n_rows <- nrow(d1)

  total_cells <- length(x$imputed_mask)
  n_imp_cells <- sum(x$imputed_mask)
  pct <- if (total_cells > 0) 100 * n_imp_cells / total_cells else 0

  method_label <- switch(x$draws_method %||% "mc_dropout",
    mc_dropout = if (x$m_per_tree == 1L) {
      "mc_dropout method; single deterministic prediction"
    } else {
      "mc_dropout method; baseline perturbations with optional GNN dropout"
    },
    conformal  = "conformal method; trait-specific sampling",
    x$draws_method
  )
  cat("pigauto experimental posterior-tree prediction sensitivity\n")
  cat(sprintf("  Trees     : %d posterior phylogenies\n", T_trees))
  cat(sprintf("  Per tree  : %d dataset(s) (%s)\n",
              m_per_tree, method_label))
  cat(sprintf("  Total     : %d completed datasets\n", M_total))
  cat(sprintf("  Rows      : %d\n", n_rows))
  cat(sprintf("  Traits    : %d -- %s\n", length(trait_cols),
              paste(trait_cols, collapse = ", ")))
  cat(sprintf("  Cells     : %d targeted / %d total (%.1f%%)\n",
              n_imp_cells, total_cells, pct))

  cat("\n  Access datasets:           mi$datasets[[i]]\n")
  cat("  Tree index for dataset i: mi$tree_index[i]\n")
  cat("  Downstream inference:     unsupported for this path\n")
  invisible(x)
}
