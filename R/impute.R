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
#' Pass \code{trait_types} through \code{...} — it is forwarded to
#' \code{\link{preprocess_traits}}.
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
                   epochs = 2000L, verbose = TRUE, seed = 1L, ...) {

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
  baseline <- fit_baseline(pd, tree, splits = splits, graph = graph)

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
    data     = pd,
    tree     = tree,
    splits   = splits,
    graph    = graph,
    baseline = baseline,
    epochs   = as.integer(epochs),
    verbose  = verbose,
    seed     = as.integer(seed),
    ...
  )

  # 6. Predict
  pred <- predict(fit, return_se = TRUE,
                  n_imputations = as.integer(n_imputations))

  # 7. Build the completed data.frame: observed values preserved,
  #    missing cells filled with model predictions.  This is the primary
  #    user-facing output.
  completed_info <- build_completed(traits, pred$imputed, species_col)
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
      splits       = splits,
      evaluation   = evaluation
    ),
    class = "pigauto_result"
  )
}


# ---- Internal: merge observed values with model predictions ---------------
# Returns the original traits data.frame with NAs replaced by imputed
# values, plus a boolean mask of which cells were filled.
build_completed <- function(original, imputed, species_col = NULL) {
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

  # 2. Row alignment: impute() preserves input order when no species_col
  #    is supplied (row names align). With species_col, multiple rows
  #    per species may exist; imputed has one row per observation.
  if (is.null(species_col)) {
    # Match by rowname
    imp_row <- match(rownames(original), rownames(imputed))
  } else {
    # Multi-obs: imputed rows align 1:1 with original rows in input order
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
