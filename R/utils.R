# Internal utility functions — not exported

rmse_vec <- function(truth, pred) {
  sqrt(mean((truth - pred)^2, na.rm = TRUE))
}

all_grads_finite <- function(parameters) {
  for (pp in parameters) {
    g <- pp$grad
    if (!is.null(g)) {
      if (!isTRUE(torch::torch_isfinite(g)$all()$item())) return(FALSE)
    }
  }
  TRUE
}

logit <- function(p) {
  p <- pmin(pmax(p, 1e-7), 1 - 1e-7)
  log(p / (1 - p))
}

expit <- function(x) {
  1 / (1 + exp(-x))
}

softmax_rows <- function(mat) {
  e <- exp(mat - apply(mat, 1, max))
  e / rowSums(e)
}

accuracy_vec <- function(truth, pred) {
  ok <- !is.na(truth) & !is.na(pred)
  if (sum(ok) == 0) return(NA_real_)
  mean(truth[ok] == pred[ok])
}

brier_vec <- function(truth_binary, pred_prob) {
  ok <- !is.na(truth_binary) & !is.na(pred_prob)
  if (sum(ok) == 0) return(NA_real_)
  mean((pred_prob[ok] - truth_binary[ok])^2)
}

mae_vec <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (sum(ok) == 0) return(NA_real_)
  mean(abs(truth[ok] - pred[ok]))
}

# Helper: get trait type for a given latent column index
get_trait_for_latent_col <- function(j, trait_map) {
  for (tm in trait_map) {
    if (j %in% tm$latent_cols) return(tm)
  }
  NULL
}

# Helper: get all latent column indices for traits of a given type
get_latent_cols_by_type <- function(trait_map, type) {
  cols <- integer(0)
  for (tm in trait_map) {
    if (tm$type == type) cols <- c(cols, tm$latent_cols)
  }
  cols
}


#' Compute a confusion matrix for categorical or binary predictions
#'
#' @param truth factor or character vector of true classes.
#' @param predicted factor or character vector of predicted classes.
#' @param levels character vector of class levels (default: union of
#'   truth and predicted levels).
#' @return A list with:
#'   \describe{
#'     \item{table}{Confusion matrix (rows = truth, columns = predicted).}
#'     \item{accuracy}{Overall accuracy.}
#'     \item{per_class}{Data.frame with per-class precision, recall, F1.}
#'   }
#' @export
confusion_matrix <- function(truth, predicted, levels = NULL) {
  if (is.null(levels)) {
    levels <- sort(union(unique(as.character(truth)),
                         unique(as.character(predicted))))
  }
  truth     <- factor(truth, levels = levels)
  predicted <- factor(predicted, levels = levels)

  tab <- table(Truth = truth, Predicted = predicted)
  acc <- sum(diag(tab)) / sum(tab)

  per_class <- data.frame(
    class     = levels,
    tp        = integer(length(levels)),
    precision = double(length(levels)),
    recall    = double(length(levels)),
    f1        = double(length(levels)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(levels)) {
    tp <- tab[i, i]
    fp <- sum(tab[, i]) - tp
    fn <- sum(tab[i, ]) - tp
    per_class$tp[i]        <- tp
    per_class$precision[i] <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    per_class$recall[i]    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    p <- per_class$precision[i]
    r <- per_class$recall[i]
    per_class$f1[i] <- if (!is.na(p) && !is.na(r) && (p + r) > 0) {
      2 * p * r / (p + r)
    } else {
      NA_real_
    }
  }

  list(table = tab, accuracy = acc, per_class = per_class)
}


#' Compute calibration data for probability predictions
#'
#' Bins predicted probabilities and computes observed frequencies within
#' each bin.  Useful for calibration plots of binary classifiers.
#'
#' @param truth numeric vector (0/1 or TRUE/FALSE).
#' @param prob numeric vector of predicted probabilities.
#' @param n_bins integer, number of bins (default 10).
#' @return A data.frame with columns \code{bin_mid}, \code{obs_freq},
#'   \code{mean_pred}, \code{n}.
#' @export
calibration_df <- function(truth, prob, n_bins = 10L) {
  ok <- !is.na(truth) & !is.na(prob)
  truth <- as.numeric(truth[ok])
  prob  <- prob[ok]

  breaks <- seq(0, 1, length.out = n_bins + 1L)
  bins   <- cut(prob, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  result <- data.frame(
    bin_mid   = double(n_bins),
    obs_freq  = double(n_bins),
    mean_pred = double(n_bins),
    n         = integer(n_bins)
  )

  for (b in seq_len(n_bins)) {
    idx <- which(bins == b)
    result$bin_mid[b]   <- (breaks[b] + breaks[b + 1]) / 2
    result$n[b]         <- length(idx)
    result$obs_freq[b]  <- if (length(idx) > 0) mean(truth[idx]) else NA_real_
    result$mean_pred[b] <- if (length(idx) > 0) mean(prob[idx]) else NA_real_
  }

  result[result$n > 0, ]
}
