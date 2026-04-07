#' Impute missing phylogenetic traits (convenience wrapper)
#'
#' One-call interface to the full pigauto pipeline: preprocessing, baseline
#' fitting, GNN training, and prediction.  For fine-grained control, use the
#' individual functions (\code{\link{preprocess_traits}},
#' \code{\link{fit_baseline}}, \code{\link{fit_pigauto}}, etc.) directly.
#'
#' @param traits data.frame with species as rownames and trait columns.
#'   Supported column types: numeric (continuous), integer (count),
#'   factor (binary/categorical), ordered (ordinal), logical (binary).
#' @param tree object of class \code{"phylo"}.
#' @param log_transform logical. Auto-log positive continuous columns
#'   (default \code{TRUE}).
#' @param missing_frac numeric. Fraction of observed cells held out for
#'   validation/test evaluation (default \code{0.25}).  Set to \code{0}
#'   to skip splitting (all cells used for training, no evaluation).
#' @param n_imputations integer. Number of MC-dropout imputation sets
#'   (default \code{1}).  Values > 1 enable between-imputation uncertainty.
#' @param epochs integer. Maximum GNN training epochs (default \code{2000}).
#' @param verbose logical. Print progress (default \code{TRUE}).
#' @param seed integer. Random seed (default \code{1}).
#' @param ... additional arguments passed to \code{\link{fit_pigauto}}.
#' @return An object of class \code{"pigauto_result"} with components:
#'   \describe{
#'     \item{prediction}{A \code{pigauto_pred} object from
#'       \code{\link{predict.pigauto_fit}}.}
#'     \item{fit}{The trained \code{pigauto_fit} object.}
#'     \item{baseline}{The phylogenetic baseline.}
#'     \item{data}{The preprocessed \code{pigauto_data} object.}
#'     \item{splits}{The val/test splits (or \code{NULL} if
#'       \code{missing_frac = 0}).}
#'     \item{evaluation}{Evaluation metrics on test set (or \code{NULL}).}
#'   }
#' @examples
#' \dontrun{
#' result <- impute(avonet300, tree300)
#' result$prediction$imputed$Mass       # imputed values
#' result$prediction$se[, "Mass"]       # standard errors
#' }
#' @export
impute <- function(traits, tree, log_transform = TRUE,
                   missing_frac = 0.25, n_imputations = 1L,
                   epochs = 2000L, verbose = TRUE, seed = 1L, ...) {

  # 1. Preprocess
  pd <- preprocess_traits(traits, tree, log_transform = log_transform)

  # 2. Create val/test splits (unless user opts out)
  if (missing_frac > 0) {
    splits <- make_missing_splits(pd$X_scaled, missing_frac = missing_frac,
                                  seed = seed, trait_map = pd$trait_map)
  } else {
    splits <- NULL
  }

  # 3. Build phylogenetic graph
  graph <- build_phylo_graph(tree, k_eigen = 8L)

  # 4. Fit phylogenetic baseline
  baseline <- fit_baseline(pd, tree, splits = splits)

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

  # 7. Evaluate on test set (if splits exist)
  evaluation <- NULL
  if (!is.null(splits) && length(splits$test_idx) > 0) {
    evaluation <- tryCatch(
      evaluate_imputation(pred, pd$X_scaled, splits),
      error = function(e) NULL
    )
  }

  structure(
    list(
      prediction = pred,
      fit        = fit,
      baseline   = baseline,
      data       = pd,
      splits     = splits,
      evaluation = evaluation
    ),
    class = "pigauto_result"
  )
}


#' @export
print.pigauto_result <- function(x, ...) {
  cat("pigauto imputation result\n")
  cat("  Species:", length(x$data$species_names), "\n")
  cat("  Traits: ", length(x$data$trait_map), "\n")
  if (!is.null(x$evaluation)) {
    test_df <- x$evaluation[x$evaluation$split == "test", ]
    if (nrow(test_df) > 0) {
      cat("  Test-set metrics:\n")
      for (i in seq_len(nrow(test_df))) {
        row <- test_df[i, ]
        if (row$type == "continuous") {
          cat(sprintf("    %s: RMSE=%.4f, r=%.4f\n",
                      row$trait, row$rmse, row$pearson_r))
        } else if (row$type == "binary") {
          cat(sprintf("    %s: accuracy=%.3f\n", row$trait, row$accuracy))
        } else if (row$type == "categorical") {
          cat(sprintf("    %s: accuracy=%.3f\n", row$trait, row$accuracy))
        }
      }
    }
  }
  invisible(x)
}
