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
