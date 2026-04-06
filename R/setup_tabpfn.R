#' Install Python dependencies for TabPFN imputation
#'
#' Creates or updates a Python virtual environment with the packages
#' needed for tabular foundation model imputation via
#' \code{\link{fit_baseline_tabpfn}}.
#'
#' @param envname Character. Name of the virtualenv (default
#'   \code{"r-tabpfn"}).
#' @param ... Additional arguments passed to
#'   \code{\link[reticulate]{py_install}}.
#' @return Invisible \code{NULL}. Called for its side effect.
#' @seealso \code{\link{fit_baseline_tabpfn}}
#' @examples
#' \dontrun{
#' setup_tabpfn()
#' }
#' @export
setup_tabpfn <- function(envname = "r-tabpfn", ...) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required. ",
         "Install with install.packages('reticulate').")
  }
  reticulate::py_install(
    packages = c("tabimpute", "numpy"),
    envname  = envname,
    pip      = TRUE,
    ...
  )
  message(
    "TabPFN imputation backend installed in virtualenv '", envname, "'.\n",
    "Activate with: reticulate::use_virtualenv('", envname, "')"
  )
  invisible(NULL)
}
