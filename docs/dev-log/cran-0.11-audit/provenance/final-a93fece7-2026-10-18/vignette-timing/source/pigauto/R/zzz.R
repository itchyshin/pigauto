.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    packageStartupMessage(
      "pigauto requires the torch package.\n",
      "Install it with:\n",
      "  install.packages('torch')\n",
      "  torch::install_torch()"
    )
  }
}
