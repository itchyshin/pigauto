skip_if_no_libtorch <- function() {
  testthat::skip_if_not_installed("torch")

  backend_ready <- tryCatch(
    isTRUE(torch::torch_is_installed()),
    error = function(...) FALSE
  )
  testthat::skip_if(!backend_ready, "torch backend not installed")

  invisible()
}
