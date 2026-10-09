args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L ||
          (length(args) == 3L && identical(args[[2L]], "installed")),
          args[[1L]] %in% c("fixtures", "real", "drm-only", "gllvm-only",
                           "neither", "dependencies"))
mode <- args[[1L]]
installed <- length(args) == 3L
if (installed) {
  stopifnot(mode != "dependencies")
  audit_lib <- normalizePath(args[[3L]], mustWork = TRUE)
  stopifnot(identical(normalizePath(.libPaths()[[1L]]), audit_lib),
            file.exists(file.path(audit_lib, "pigauto", "DESCRIPTION")))
}

if (mode == "dependencies") {
  dcf <- read.dcf("DESCRIPTION")[1L, ]
  dependencies_absent <- function(fields) {
    !any(grepl("(^|[^[:alnum:]_.])(drmTMB|gllvmTMB)([^[:alnum:]_.]|$)",
               fields, perl = TRUE))
  }
  positive_control <- dcf
  positive_control[["Suggests"]] <- paste(dcf[["Suggests"]], "drmTMB", sep = ", ")
  stopifnot(dependencies_absent(dcf), !dependencies_absent(positive_control))
  cat("OPTIONAL_DEPENDENCIES_OK\n")
} else {
  Sys.setenv(NOT_CRAN = "true")
  if (installed) {
    library(pigauto, lib.loc = audit_lib)
    stopifnot(identical(normalizePath(getNamespaceInfo("pigauto", "path")),
                        normalizePath(file.path(audit_lib, "pigauto"))),
              as.character(packageVersion("pigauto")) == "0.11.0")
    Sys.setenv(PIGAUTO_AUDIT_INSTALLED_LIB = audit_lib)
  } else {
    pkgload::load_all(".", quiet = TRUE)
    Sys.unsetenv("PIGAUTO_AUDIT_INSTALLED_LIB")
  }
  files <- if (mode == "fixtures") {
    c("tests/testthat/test-pool-mi-adapter-fixtures.R",
      "tests/testthat/test-mi-provenance.R")
  } else {
    available <- c(drmTMB = requireNamespace("drmTMB", quietly = TRUE),
                   gllvmTMB = requireNamespace("gllvmTMB", quietly = TRUE))
    expected <- switch(mode,
                       real = c(drmTMB = TRUE, gllvmTMB = TRUE),
                       `drm-only` = c(drmTMB = TRUE, gllvmTMB = FALSE),
                       `gllvm-only` = c(drmTMB = FALSE, gllvmTMB = TRUE),
                       neither = c(drmTMB = FALSE, gllvmTMB = FALSE))
    stopifnot(identical(available, expected))
    "script/cran-0.11-integration/test-real-backends.R"
  }
  for (file in files) {
    results <- as.data.frame(testthat::test_file(
      file, reporter = "silent", package = if (installed) "pigauto" else NULL
    ))
    stopifnot(nrow(results) > 0L,
              all(results$failed == 0L),
              !any(results$error),
              all(results$warning == 0L))
    if (mode == "fixtures") {
      stopifnot(!any(results$skipped))
    } else {
      stopifnot(sum(results$skipped) == sum(!expected))
    }
    cat(file, ": ", sum(results$passed), " expectations passed",
        if (installed) " against installed pigauto" else " against source pigauto",
        "\n", sep = "")
  }
  if (identical(mode, "neither")) {
    check_missing_backend_for_saved_class <- function(package, fit) {
      fit_path <- tempfile(fileext = ".rds")
      on.exit(unlink(fit_path), add = TRUE)
      saveRDS(list(fit, fit), fit_path)
      fits <- readRDS(fit_path)
      stopifnot(identical(class(fits[[1L]]), class(fit)),
                !requireNamespace(package, quietly = TRUE))
      error <- tryCatch(pool_mi(fits), error = identity)
      stopifnot(inherits(error, "error"),
                grepl(paste0("automatic ", package,
                             " adapter requires the package"),
                      conditionMessage(error), fixed = TRUE))
    }
    check_missing_backend_for_saved_class(
      "drmTMB",
      structure(list(opt = list(convergence = 0L),
                     sdr = list(pdHess = TRUE),
                     beta = c(x = 1),
                     V = matrix(0.04, 1L, dimnames = list("x", "x"))),
                class = "drmTMB")
    )
    check_missing_backend_for_saved_class(
      "gllvmTMB",
      structure(list(opt = list(convergence = 0L),
                     sd_report = list(pdHess = TRUE)),
                class = "gllvmTMB_multi")
    )
    cat("SAVED_OPTIONAL_CLASS_MISSING_PACKAGE_OK\n")
  }
  cat(switch(mode,
             fixtures = "ADAPTER_FIXTURES_OK\n",
             real = "REAL_ADAPTERS_OK\n",
             `drm-only` = "DRM_ONLY_OK\n",
             `gllvm-only` = "GLLVM_ONLY_OK\n",
             neither = "NEITHER_BACKEND_OK\n"))
}
