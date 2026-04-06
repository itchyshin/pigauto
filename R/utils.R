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
