args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
root <- args[1]; fixture <- args[2]; receipt <- args[3]; expected_file <- args[4]
threads <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))
if (!all(threads == "1")) stop("child requires one-thread environment")
cat("RESOURCE_THREADS", paste(names(threads), threads, collapse = " "), "\n")
if (file.exists(file.path(root, "Meta", "package.rds"))) {
  library(pigauto, lib.loc = dirname(root))
} else {
  pkgload::load_all(root, quiet = TRUE)
}
fn <- getFromNamespace("exact_conditional_mvn", "pigauto")
source_path <- getSrcFilename(fn, full.names = TRUE)
if (!length(source_path) || !file.exists(source_path)) stop("source fingerprint unavailable")
body_file <- tempfile()
writeLines(deparse(body(fn), width.cutoff = 500L), body_file)
identity <- list(namespace = normalizePath(getNamespaceInfo(asNamespace("pigauto"), "path")),
  source = normalizePath(source_path), source_md5 = unname(tools::md5sum(source_path)),
  body_md5 = unname(tools::md5sum(body_file)))
unlink(body_file)
cat("RESOURCE_IDENTITY", paste(unlist(identity), collapse = " "), "\n")
expected <- readRDS(expected_file)
if (!identical(identity, expected)) stop("resource child source or namespace identity mismatch")
# Probe the limit before allocating the larger fixture; restore it on all paths.
limit_available <- function() {
  old <- mem.maxVSize()
  on.exit({
    restored <- mem.maxVSize(old)
    if (!identical(restored, old)) stop("vector heap limit restoration failed")
  }, add = TRUE)
  target <- ceiling(gc()["Vcells", "gc trigger"] * 8 / 1024^2) + 128
  if (is.finite(old) && old < target) return(FALSE)
  active <- withCallingHandlers(mem.maxVSize(target), warning = function(w) {
    if (grepl("a limit lower than current usage, so ignored", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
  identical(active, as.double(target))
}
if (!limit_available()) {
  saveRDS(list(identity = identity, status = "prerequisite",
    reason = "resource child vector heap limit unavailable"), receipt)
  quit(save = "no", status = 0L)
}
r <- testthat::test_file(fixture, reporter = "summary")
d <- as.data.frame(r)
if (!nrow(d) || any(d$failed) || any(d$error) || any(d$warning))
  stop("child tests did not pass cleanly")
if (any(d$skipped)) {
  skips <- unlist(lapply(r, function(z) Filter(function(x) inherits(x, "expectation_skip"), z$results)), recursive = FALSE)
  if (length(skips) != 1L || !identical(conditionMessage(skips[[1]]),
      "Reason: resource child vector heap limit unavailable")) stop("unexpected child skip")
  saveRDS(list(identity = identity, status = "prerequisite",
    reason = "resource child vector heap limit unavailable"), receipt)
} else {
  if (!sum(d$passed)) stop("child tests had no assertions")
  saveRDS(list(identity = identity, status = "pass", passed = sum(d$passed)), receipt)
}
