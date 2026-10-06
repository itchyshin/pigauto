args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L, args[[4L]] %in% c("source", "installed"))
fit_path <- args[[1L]]
backend <- args[[2L]]
location <- args[[3L]]
mode <- args[[4L]]
if (backend %in% loadedNamespaces()) {
  stop("The optional backend namespace was loaded before the saved fit")
}
if (mode == "installed") {
  stopifnot(identical(normalizePath(.libPaths()[[1L]]),
                      normalizePath(location)))
  library(pigauto, lib.loc = location)
  stopifnot(identical(normalizePath(getNamespaceInfo("pigauto", "path")),
                      normalizePath(file.path(location, "pigauto"))))
} else {
  pkgload::load_all(location, quiet = TRUE)
}
if (backend %in% loadedNamespaces()) {
  stop("Loading pigauto loaded the optional backend namespace")
}
fits <- readRDS(fit_path)
if (paste0("package:", backend) %in% search()) {
  stop("Reading the saved fit attached the optional backend package")
}
out <- suppressWarnings(pool_mi(fits))
if (!backend %in% loadedNamespaces()) {
  stop("Pooling did not load the optional backend namespace")
}
stopifnot(inherits(out, "pigauto_pooled"), all(is.finite(out$std.error)))
cat("RELOAD_POOL_OK", backend, "\n")
