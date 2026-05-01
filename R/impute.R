#' Impute missing phylogenetic traits (convenience wrapper)
#'
#' One-call interface to the full pigauto pipeline: preprocessing, baseline
#' fitting, GNN training, and prediction.  For fine-grained control, use the
#' individual functions (\code{\link{preprocess_traits}},
#' \code{\link{fit_baseline}}, \code{\link{fit_pigauto}}, etc.) directly.
#'
#' @param traits data.frame with species as rownames and trait columns,
#'   or (when \code{species_col} is supplied) a data.frame with a species
#'   column that may have multiple rows per species.
#'   Supported column types: \code{numeric} (continuous), \code{integer}
#'   (count), \code{factor} (binary/categorical), \code{ordered} (ordinal),
#'   \code{character} (factor → binary/categorical), \code{logical} (binary).
#'   See the \strong{Trait type auto-detection} section below.
#' @param tree object of class \code{"phylo"}.
#' @param species_col character. Name of the column in \code{traits} that
#'   identifies species.  When supplied, multiple observations per species
#'   are supported.  Default \code{NULL} uses row names (one row per species).
#' @param trait_types named character vector overriding the auto-detected
#'   type for specific trait columns, e.g.
#'   \code{c(Survival = "proportion", Parasites = "zi_count")}. Required for
#'   the two types that cannot be inferred from R class (see
#'   \strong{Trait type auto-detection} below). Default \code{NULL} (auto).
#' @param multi_proportion_groups named list declaring compositional
#'   (\code{multi_proportion}) traits, e.g.
#'   \code{list(colour = c("black", "blue", "red", "yellow"))}. Each list
#'   element names a group and gives the K trait columns that form a
#'   simplex (rows summing to 1). Encoded via CLR + per-component z-score.
#'   Multi_proportion traits \emph{cannot} be declared through
#'   \code{trait_types} — use this argument instead. Default \code{NULL}
#'   (no multi_proportion groups).
#' @param log_transform logical. Auto-log positive continuous columns
#'   (default \code{TRUE}).
#' @param missing_frac numeric. Fraction of observed cells held out for
#'   validation/test evaluation (default \code{0.25}).  Set to \code{0}
#'   to skip splitting (all cells used for training, no evaluation).
#' @param n_imputations integer. Number of MC-dropout imputation sets
#'   (default \code{1}).  Values > 1 enable between-imputation uncertainty.
#' @param covariates data.frame or matrix of environmental covariates
#'   (fully observed — no NAs).  Covariates are conditioners: they inform
#'   imputation but are not themselves imputed.  Numeric/integer columns are
#'   z-scored; factor/ordered columns are one-hot encoded automatically.
#'   If a variable has missing values, include it in \code{traits} instead.
#'   Same number of rows as \code{traits}. Default \code{NULL} (no covariates).
#' @param epochs integer. Maximum GNN training epochs (default \code{2000}).
#' @param verbose logical. Print progress (default \code{TRUE}).
#' @param seed integer. Random seed (default \code{1}).
#' @param multi_obs_aggregation character. How to aggregate multiple
#'   observations per species before the Level-C baseline.  \code{"hard"}
#'   (default) thresholds binary proportions at 0.5 and uses argmax for
#'   categorical.  \code{"soft"} preserves species-level proportions and
#'   dispatches a soft E-step so that intermediate class frequencies
#'   contribute fractional liability evidence.  Passed to
#'   \code{\link{fit_baseline}}.
#' @param em_iterations integer. Phase 6 EM iterations for the
#'   threshold-joint baseline (binary + ordinal + OVR categorical).
#'   Default \code{0L} preserves v0.9.1 behaviour byte-for-byte. When
#'   \code{>= 2L}, the BM rate \eqn{\Sigma} learned by
#'   \code{Rphylopars::phylopars()} at iteration \eqn{k} is fed back as
#'   the per-trait prior SD at iteration \eqn{k+1}, up to
#'   \code{em_iterations} times or until \code{em_tol} convergence.
#'   Passed to \code{\link{fit_baseline}}.
#' @param em_tol numeric. Relative-Frobenius convergence tolerance for
#'   the Phase 6 / 7 EM loop. Default \code{1e-3}.
#' @param em_offdiag logical. Phase 7 opt-in: when \code{TRUE} AND
#'   \code{em_iterations >= 2L}, each liability cell's prior uses the
#'   full conditional-MVN from \eqn{\Sigma}'s off-diagonal entries, so that
#'   observing one discrete trait shifts (not just tightens) the prior
#'   on correlated other traits.  Binary + ordinal only; OVR categorical
#'   stays on Phase 6 diagonal.  Default \code{FALSE}.  Passed to
#'   \code{\link{fit_baseline}}.
#' @param pool_method character. How to pool multiple imputation draws
#'   (\code{n_imputations > 1}) for count, proportion, and zi_count
#'   magnitude traits: \code{"median"} (default) takes the per-cell
#'   median of the \code{M} decoded draws — robust to dropout-noisy
#'   latents amplified by \code{expm1()} / \code{plogis()} decoders.
#'   \code{"mean"} restores the pre-v0.9.2 arithmetic-mean pooling.
#'   \code{"mode"} (Phase H, v0.9.1.9010+) is intended for ordinal
#'   traits: per-cell majority vote across the \code{M} draws,
#'   avoiding the integer-mean-round bias toward middle classes.  For
#'   continuous-family traits, \code{"mode"} falls back to
#'   \code{"median"}.  Binary / categorical / multi_proportion traits
#'   always pool by probability average; unaffected by this argument.  See
#'   \href{https://github.com/itchyshin/pigauto/issues/40}{issue #40}.
#' @param safety_floor logical. When \code{TRUE} (default since
#'   v0.9.1.9002), calibration searches the 3-way simplex
#'   \code{r_BM * BM + r_GNN * GNN + r_MEAN * MEAN} so the grand mean is
#'   always in the candidate set and the calibrated prediction is
#'   guaranteed never to be worse than the grand-mean baseline on
#'   validation.  When \code{FALSE}, the v0.9.1 1-D calibration is used
#'   exactly.  See the Safety floor section below.
#' @param phylo_signal_gate,phylo_signal_threshold,phylo_signal_method
#'   Pass-through to [fit_pigauto()]. See that help page for details.
#' @param clamp_outliers logical.  Phase G (v0.9.1.9011+).  When
#'   \code{TRUE}, post-back-transform predictions for log-transformed
#'   continuous, count, and zi_count magnitude traits are capped at
#'   \code{tm$obs_max * clamp_factor} (and \code{tm$obs_max} is the
#'   observed maximum on the original scale, recorded at preprocess
#'   time).  Targets the AVONET Mass tail-extrapolation mode
#'   documented in
#'   \code{useful/MEMO_2026-05-01_avonet_mass_diag.md} where a
#'   \eqn{+3}-\eqn{4 \sigma} latent overshoot becomes a 50x-100x
#'   value error after \code{expm1()}.  Default \code{FALSE}
#'   preserves v0.9.1 behaviour exactly.
#' @param clamp_factor numeric scalar (>= 1).  Multiplicative factor
#'   on the observed maximum used by \code{clamp_outliers}.  Default
#'   \code{5} (Tukey-style outlier definition: anything >= 5x the
#'   observed max is implausible).  Ignored when
#'   \code{clamp_outliers = FALSE}.
#' @param ... additional arguments passed to \code{\link{fit_pigauto}}.
#' @return An object of class \code{"pigauto_result"} with components:
#'   \describe{
#'     \item{completed}{The input \code{traits} data.frame with observed
#'       values preserved and only missing cells filled in.  This is the
#'       primary output -- typically what users want.}
#'     \item{imputed_mask}{Logical matrix (same shape as \code{completed})
#'       that is \code{TRUE} for cells that were imputed (originally
#'       \code{NA}) and \code{FALSE} for observed cells.}
#'     \item{prediction}{A \code{pigauto_pred} object from
#'       \code{\link{predict.pigauto_fit}} containing raw model
#'       predictions for every cell (observed + missing), standard errors,
#'       class probabilities, and conformal intervals.}
#'     \item{fit}{The trained \code{pigauto_fit} object.}
#'     \item{baseline}{The phylogenetic baseline.}
#'     \item{data}{The preprocessed \code{pigauto_data} object.}
#'     \item{splits}{The val/test splits (or \code{NULL} if
#'       \code{missing_frac = 0}).}
#'     \item{evaluation}{Evaluation metrics on test set (or \code{NULL}).}
#'   }
#' @section Trait type auto-detection:
#' pigauto infers each trait's type from its R class — no \code{trait_types}
#' argument is needed for most data:
#' \tabular{ll}{
#'   \strong{R class} \tab \strong{pigauto type} \cr
#'   \code{numeric} \tab continuous (auto-log if all-positive) \cr
#'   \code{integer} \tab count \cr
#'   \code{factor} with 2 levels \tab binary \cr
#'   \code{factor} (unordered) with >2 levels \tab categorical \cr
#'   \code{ordered} / \code{factor(..., ordered = TRUE)} \tab ordinal \cr
#'   \code{character} \tab → factor → binary or categorical \cr
#'   \code{logical} \tab binary \cr
#' }
#' Two types cannot be inferred from class alone and \emph{must} be declared
#' via \code{trait_types}:
#' \describe{
#'   \item{\code{"proportion"}}{A \code{numeric} bounded 0–1, e.g. survival
#'     rate: \code{trait_types = c(Survival = "proportion")}.}
#'   \item{\code{"zi_count"}}{An \code{integer} with excess zeros, e.g.
#'     parasite count: \code{trait_types = c(Parasites = "zi_count")}.}
#' }
#' Use the \code{trait_types} argument directly (it is an explicit
#' parameter, not a \code{...} pass-through).
#'
#' @section Traits vs covariates:
#' The distinction is \strong{functional, not ontological}: a trait is something
#' you want to impute (NA values allowed in \code{traits}); a covariate is
#' something you use to sharpen imputation accuracy (must be fully observed,
#' passed via \code{covariates}).  The same variable can be either depending on
#' the scientific question.
#'
#' Examples:
#' \itemize{
#'   \item \strong{IUCN status with Data Deficient species} → put it in
#'     \code{traits} as \code{ordered(c("LC","NT","VU","EN","CR"))} so
#'     pigauto predicts the unknown categories.
#'   \item \strong{IUCN status fully known for all species} → pass as a
#'     covariate to inform imputation of other traits (e.g. body mass,
#'     range size).
#'   \item \strong{Realm / biome (factor)} → pass as a covariate; pigauto
#'     one-hot encodes factor columns automatically (v0.6.1+).
#' }
#'
#' Variables that belong in \code{traits}: anything with missing values you
#' care about predicting.  Variables that belong in \code{covariates}: fully
#' observed, exogenous to the trait space (geography, climate, habitat,
#' experimental treatment).
#'
#' @section Safety floor (v0.9.1.9002+):
#'   With \code{safety_floor = TRUE} (the new default), the post-training
#'   calibration grid searches a 3-way convex combination of the
#'   Brownian-motion baseline, the GNN delta, and the per-trait grand
#'   mean.  The simplex is sampled at step 0.05 (231 candidates per latent
#'   column).  Because the corner \code{(0, 0, 1)} — pure grand mean —
#'   is always in the grid, the calibrated validation RMSE is guaranteed
#'   by construction to satisfy \code{calibrated_val_RMSE <=
#'   mean_val_RMSE} on every trait.  The fit object gains four new slots:
#'   \code{r_cal_bm}, \code{r_cal_gnn}, \code{r_cal_mean} (each a named
#'   numeric of length \code{p_latent}), and
#'   \code{mean_baseline_per_col}.
#'
#'   Set \code{safety_floor = FALSE} to reproduce the pre-v0.9.1.9002
#'   1-D calibration bit-identically (no mean term; \code{r_cal_mean = 0};
#'   \code{r_cal_bm = 1 - r_cal_gnn}).  See
#'   \code{specs/2026-04-23-safety-floor-mean-gate-design.md} for the
#'   design rationale and
#'   \code{plans/2026-04-23-safety-floor-mean-gate.md} for the
#'   implementation plan.
#'
#' @examples
#' \dontrun{
#' # Simple case: fill in missing values — type detection is automatic
#' result <- impute(my_traits, my_tree)
#' result$completed              # observed preserved, NAs filled
#' result$imputed_mask           # which cells were imputed
#'
#' # Proportion and zi_count must be declared explicitly
#' result <- impute(my_traits, my_tree,
#'                  trait_types = c(Survival  = "proportion",
#'                                  Parasites = "zi_count"))
#'
#' # Diagnostic: raw predictions for every cell
#' result$prediction$imputed     # model's prediction everywhere
#' result$prediction$se          # per-cell uncertainty
#' result$prediction$probabilities$diet  # class probabilities
#' }
#' @export
impute <- function(traits, tree, species_col = NULL,
                   trait_types = NULL,
                   multi_proportion_groups = NULL,
                   log_transform = TRUE,
                   missing_frac = 0.25, n_imputations = 1L,
                   covariates = NULL,
                   epochs = 2000L, verbose = TRUE, seed = 1L,
                   multi_obs_aggregation = c("hard", "soft"),
                   em_iterations = 0L,
                   em_tol = 1e-3,
                   em_offdiag = FALSE,
                   pool_method = c("median", "mean", "mode"),
                   clamp_outliers = FALSE,
                   clamp_factor = 5,
                   safety_floor = TRUE,
                   phylo_signal_gate = TRUE,
                   phylo_signal_threshold = 0.2,
                   phylo_signal_method = "lambda",
                   ...) {
  multi_obs_aggregation <- match.arg(multi_obs_aggregation)
  pool_method <- match.arg(pool_method)

  # 1. Preprocess
  pd <- preprocess_traits(traits, tree, species_col = species_col,
                          trait_types = trait_types,
                          multi_proportion_groups = multi_proportion_groups,
                          log_transform = log_transform,
                          covariates = covariates)

  # 2. Create val/test splits (unless user opts out).  We give the
  #    validation set a larger share (50%) than the default so that
  #    post-hoc gate calibration has enough cells to detect real
  #    improvements.
  if (missing_frac > 0) {
    splits <- make_missing_splits(pd$X_scaled, missing_frac = missing_frac,
                                  val_frac = 0.5,
                                  seed = seed, trait_map = pd$trait_map)
  } else {
    splits <- NULL
  }

  # 3. Build phylogenetic graph (auto k_eigen scales with tree size).
  #    The graph object carries the cophenetic distance matrix in
  #    graph$D; we pass it to fit_baseline() below so that it is
  #    computed exactly once for the whole pipeline.
  graph <- build_phylo_graph(tree, k_eigen = "auto")

  # 4. Fit phylogenetic baseline (reuses graph$D for label propagation)
  baseline <- fit_baseline(pd, tree, splits = splits, graph = graph,
                           multi_obs_aggregation = multi_obs_aggregation,
                           em_iterations = em_iterations,
                           em_tol = em_tol,
                           em_offdiag = em_offdiag)

  # Free the cached cophenetic distance matrix: fit_pigauto() only
  # needs graph$adj and graph$coords, and at n = 10,000 the ~800 MB
  # D matrix held in R memory during training caused a large
  # per-epoch slowdown. This is safe because no downstream caller
  # reads graph$D -- see the v0.3.1 NEWS entry.
  graph$D     <- NULL
  graph$R_phy <- NULL
  invisible(gc(full = TRUE, verbose = FALSE))

  # 5. Train GNN
  fit <- fit_pigauto(
    data                   = pd,
    tree                   = tree,
    splits                 = splits,
    graph                  = graph,
    baseline               = baseline,
    epochs                 = as.integer(epochs),
    verbose                = verbose,
    seed                   = as.integer(seed),
    safety_floor           = safety_floor,
    phylo_signal_gate      = phylo_signal_gate,
    phylo_signal_threshold = phylo_signal_threshold,
    phylo_signal_method    = phylo_signal_method,
    ...
  )

  # Belt-and-braces GPU memory reclaim before predict().  fit_pigauto()
  # already moves its state_dict to CPU and calls cuda_empty_cache()
  # internally, but doing it again here handles any R-level references
  # to graph/baseline tensors that may linger between calls.  Essential
  # at n >= 5000 on cards with <= 46 GB to avoid OOM on the first
  # predict-stage allocation.
  invisible(gc(full = TRUE, verbose = FALSE))
  if (torch::cuda_is_available()) {
    try(torch::cuda_empty_cache(), silent = TRUE)
  }

  # 6. Predict
  pred <- predict(fit, return_se = TRUE,
                  n_imputations = as.integer(n_imputations),
                  pool_method = pool_method,
                  clamp_outliers = clamp_outliers,
                  clamp_factor   = clamp_factor)

  # 7. Build the completed data.frame: observed values preserved,
  #    missing cells filled with model predictions.  This is the primary
  #    user-facing output.
  #
  # `pd$input_row_order[k]` tells us the original-input row index for
  # internal position k.  build_completed needs this whenever
  # preprocess_traits reordered rows (multi-obs always; single-obs when
  # input row order didn't already match tree-tip order).
  completed_info <- build_completed(traits, pred$imputed, species_col,
                                     input_row_order = pd$input_row_order)
  completed      <- completed_info$completed
  imputed_mask   <- completed_info$imputed_mask

  # 8. Evaluate on test set (if splits exist)
  evaluation <- NULL
  if (!is.null(splits) && length(splits$test_idx) > 0) {
    evaluation <- tryCatch(
      evaluate_imputation(pred, pd$X_scaled, splits),
      error = function(e) NULL
    )
  }

  structure(
    list(
      completed    = completed,
      imputed_mask = imputed_mask,
      prediction   = pred,
      fit          = fit,
      baseline     = baseline,
      data         = pd,
      tree         = tree,
      splits       = splits,
      evaluation   = evaluation
    ),
    class = "pigauto_result"
  )
}


# ---- Internal: merge observed values with model predictions ---------------
# Returns the original traits data.frame with NAs replaced by imputed
# values, plus a boolean mask of which cells were filled.
#
# `input_row_order` (from `pd$input_row_order`) maps internal positions back
# to original input rows: imputed[k, ] corresponds to original[input_row_order[k], ].
# When supplied (post 2026-04-26 row-alignment fix) it is used to recover the
# user's row order; passing NULL falls back to the legacy match-by-rowname /
# 1-to-1 alignment behaviour for backward compat.
build_completed <- function(original, imputed, species_col = NULL,
                             input_row_order = NULL) {
  # Align row order: imputed has rownames = species / obs ids used by the fit
  # Original may have row names OR a species column. We match on the
  # column types.

  # 1. Start from the original so columns are in user order and
  #    non-trait columns (if any, e.g. species_col) are preserved.
  completed <- original
  n_row <- nrow(original)
  trait_cols <- setdiff(names(original), species_col)

  imputed_mask <- matrix(FALSE, nrow = n_row, ncol = length(trait_cols),
                         dimnames = list(rownames(original), trait_cols))

  # 2. Row alignment.
  #
  # Pre-2026-04-26: multi-obs assumed `imputed` rows aligned 1:1 with input
  # rows ("Multi-obs: imputed rows align 1:1 with original rows in input
  # order"), but `preprocess_traits()` reorders rows internally to tree-tip
  # order, breaking that assumption when the input data.frame's species
  # column isn't already tree-tip-sorted.  Result: predictions for masked
  # cells were assigned to the wrong rows in `completed`, producing
  # near-zero correlation with truth on shuffled-input data.
  #
  # Fix: when `input_row_order` is supplied, use it to invert the reorder.
  # input_row_order[k] gives the original input-row index for internal
  # position k; we want the inverse: for each original row i, find the
  # internal position k such that input_row_order[k] == i.  That is
  # `match(seq_len(n_row), input_row_order)`, which is robust to the
  # synthetic NA rows that `preprocess_traits` appends for tree tips that
  # have no input data (those have input_row_order == NA and are simply
  # not matched).
  if (!is.null(input_row_order)) {
    imp_row <- match(seq_len(n_row), input_row_order)
  } else if (is.null(species_col)) {
    # Legacy path: match by rowname (single-obs without input_row_order)
    imp_row <- match(rownames(original), rownames(imputed))
  } else {
    # Legacy path: assume 1:1 (kept for backward compat with any external
    # callers of build_completed; not used by impute() itself anymore).
    imp_row <- seq_len(n_row)
  }

  for (nm in trait_cols) {
    if (!(nm %in% names(imputed))) next
    obs_col <- original[[nm]]
    imp_col <- imputed[[nm]][imp_row]
    missing_i <- is.na(obs_col)
    imputed_mask[, nm] <- missing_i

    if (is.factor(obs_col)) {
      # Preserve factor levels and ordering from the original
      new_col <- obs_col
      if (any(missing_i)) {
        # imp_col is a factor with potentially the same levels
        imp_chars <- as.character(imp_col)
        new_col <- as.character(obs_col)
        new_col[missing_i] <- imp_chars[missing_i]
        new_col <- factor(new_col, levels = levels(obs_col),
                          ordered = is.ordered(obs_col))
      }
      completed[[nm]] <- new_col

    } else if (is.logical(obs_col)) {
      new_col <- obs_col
      if (any(missing_i)) {
        imp_as_logical <- if (is.factor(imp_col)) {
          as.character(imp_col) %in% c("TRUE", "yes", "1")
        } else {
          as.logical(imp_col)
        }
        new_col[missing_i] <- imp_as_logical[missing_i]
      }
      completed[[nm]] <- new_col

    } else if (is.integer(obs_col)) {
      new_col <- obs_col
      if (any(missing_i)) {
        new_col[missing_i] <- as.integer(round(imp_col[missing_i]))
      }
      completed[[nm]] <- new_col

    } else {
      # numeric / double
      new_col <- obs_col
      if (any(missing_i)) {
        new_col[missing_i] <- as.numeric(imp_col[missing_i])
      }
      completed[[nm]] <- new_col
    }
  }

  list(completed = completed, imputed_mask = imputed_mask)
}


#' @export
print.pigauto_result <- function(x, ...) {
  cat("pigauto imputation result\n")
  cat("  Species :", length(x$data$species_names), "\n")
  cat("  Traits  :", length(x$data$trait_map),
      "--", paste(vapply(x$data$trait_map, "[[", character(1), "name"),
                  collapse = ", "), "\n")

  # Trait type counts
  types <- vapply(x$data$trait_map, "[[", character(1), "type")
  type_tab <- table(types)
  cat("  Types   :",
      paste(names(type_tab), type_tab, sep = "=", collapse = ", "), "\n")

  # Imputation summary
  if (!is.null(x$imputed_mask)) {
    per_trait <- colSums(x$imputed_mask)
    total_imp <- sum(per_trait)
    total_cells <- length(x$imputed_mask)
    cat(sprintf("  Cells   : %d imputed / %d total (%.1f%%)\n",
                total_imp, total_cells, 100 * total_imp / total_cells))
    if (total_imp > 0) {
      cat("  Per-trait imputed cell counts:\n")
      for (nm in names(per_trait)) {
        if (per_trait[[nm]] > 0) {
          cat(sprintf("    %-20s %d\n", nm, per_trait[[nm]]))
        }
      }
    }
  }

  # Evaluation metrics (held-out test set)
  if (!is.null(x$evaluation)) {
    test_df <- x$evaluation[x$evaluation$split == "test", ]
    if (nrow(test_df) > 0) {
      cat("\n  Held-out test-set metrics:\n")
      cat(sprintf("    %-20s %-12s %10s %10s\n",
                  "trait", "type", "metric", "value"))
      cat(sprintf("    %-20s %-12s %10s %10s\n",
                  strrep("-", 20), strrep("-", 12),
                  strrep("-", 10), strrep("-", 10)))
      for (i in seq_len(nrow(test_df))) {
        row <- test_df[i, ]
        if (row$type %in% c("continuous", "count")) {
          cat(sprintf("    %-20s %-12s %10s %10.4f\n",
                      row$trait, row$type, "RMSE", row$rmse))
          if (!is.na(row$pearson_r)) {
            cat(sprintf("    %-20s %-12s %10s %10.4f\n",
                        "", "", "r", row$pearson_r))
          }
        } else if (row$type == "ordinal") {
          cat(sprintf("    %-20s %-12s %10s %10.4f\n",
                      row$trait, row$type, "rho", row$spearman_rho))
        } else if (row$type %in% c("binary", "categorical")) {
          cat(sprintf("    %-20s %-12s %10s %10.3f\n",
                      row$trait, row$type, "accuracy", row$accuracy))
        }
      }
    }
  }

  cat("\n  Access the completed data.frame with  result$completed\n")
  invisible(x)
}
