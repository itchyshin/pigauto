skip_if_no_libtorch <- function() {
  testthat::skip_if_not_installed("torch")

  # torch_is_installed() can be TRUE when the R package is present but the
  # Lantern/libtorch runtime is absent (notably on win-builder). Probe one
  # actual tensor operation so tests skip that partial-install state too.
  backend_ready <- tryCatch(
    isTRUE(torch::torch_is_installed()) && {
      probe <- torch::torch_tensor(1)
      isTRUE(as.numeric(probe$item()) == 1)
    },
    error = function(...) FALSE
  )
  testthat::skip_if(!backend_ready, "torch backend not installed")

  invisible()
}
