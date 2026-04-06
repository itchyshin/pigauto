#' Evaluate imputation performance against known values
#'
#' Computes RMSE, Pearson correlation, and 95% interval coverage for the
#' validation and test splits, returning one row per trait per split.
#'
#' @param pred matrix of predicted values (in z-score scale, as returned by
#'   the internal model, not back-transformed).
#' @param truth matrix of true values (in z-score scale).
#' @param splits list (output of \code{\link{make_missing_splits}}).
#' @param pred_se numeric matrix of prediction standard errors (same scale as
#'   \code{pred}). When supplied, 95% interval coverage is computed.
#' @return A \code{data.frame} with columns \code{split}, \code{trait},
#'   \code{n}, \code{rmse}, \code{pearson_r}, and (if \code{pred_se} is
#'   provided) \code{coverage_95}.
#' @examples
#' \dontrun{
#' eval_df <- evaluate_imputation(bl$mu, pd$X_scaled, splits)
#' }
#' @export
evaluate_imputation <- function(pred, truth, splits, pred_se = NULL) {
  if (!is.matrix(pred) || !is.matrix(truth)) {
    stop("'pred' and 'truth' must be matrices.")
  }
  if (!identical(dim(pred), dim(truth))) {
    stop("'pred' and 'truth' must have identical dimensions.")
  }

  trait_names <- colnames(truth)
  if (is.null(trait_names)) trait_names <- paste0("trait", seq_len(ncol(truth)))

  eval_split <- function(idx, label) {
    lapply(seq_len(ncol(truth)), function(j) {
      col_idx <- idx[idx %% nrow(truth) == (j - 1L) %% nrow(truth) |
                       ceiling(idx / nrow(truth)) == j]
      # Simpler: row indices for this column
      rows_j <- idx[((idx - 1L) %% nrow(truth)) + 1L %in%
                      seq_len(nrow(truth)) &
                      ceiling(idx / nrow(truth)) == j]
      # Straightforward: convert linear idx to (row, col) and filter col == j
      row_j <- ((idx - 1L) %% nrow(truth)) + 1L
      col_j <- ceiling(idx / nrow(truth))
      keep  <- which(col_j == j)
      ri    <- row_j[keep]

      if (length(ri) == 0L) return(NULL)

      p_j <- pred[ri, j]
      t_j <- truth[ri, j]
      ok  <- is.finite(p_j) & is.finite(t_j)

      rmse_j <- if (sum(ok) > 0L) rmse_vec(t_j[ok], p_j[ok]) else NA_real_
      r_j    <- if (sum(ok) > 1L) {
        stats::cor(t_j[ok], p_j[ok])
      } else {
        NA_real_
      }

      cov95 <- NA_real_
      if (!is.null(pred_se) && sum(ok) > 0L) {
        se_j  <- pred_se[ri, j]
        lower <- p_j - 1.96 * se_j
        upper <- p_j + 1.96 * se_j
        cov95 <- mean((t_j >= lower & t_j <= upper)[ok], na.rm = TRUE)
      }

      data.frame(
        split      = label,
        trait      = trait_names[j],
        n          = sum(ok),
        rmse       = rmse_j,
        pearson_r  = r_j,
        coverage_95 = cov95,
        stringsAsFactors = FALSE
      )
    }) |>
      (\(lst) do.call(rbind, lst[!sapply(lst, is.null)]))()
  }

  rows <- list()
  if (length(splits$val_idx) > 0L) {
    rows$val <- eval_split(splits$val_idx, "val")
  }
  if (length(splits$test_idx) > 0L) {
    rows$test <- eval_split(splits$test_idx, "test")
  }

  do.call(rbind, rows[!sapply(rows, is.null)])
}
