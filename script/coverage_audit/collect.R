# Instrument one frozen source export. Run on Totoro, not an interactive laptop.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
expected_sha <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
marker <- file.path(root, "SOURCE_SHA")
stopifnot(file.exists(marker), identical(readLines(marker, warn = FALSE), expected_sha))
threads <- as.integer(Sys.getenv("OMP_NUM_THREADS", "1"))
stopifnot(is.finite(threads), threads >= 1L, threads <= 12L)
source_paths <- c(
  paste0("R/", list.files(file.path(root, "R"), pattern = "\\.R$")),
  paste0("tests/", list.files(file.path(root, "tests"), recursive = TRUE,
                            pattern = "\\.R$")),
  "DESCRIPTION", "NAMESPACE", "codecov.yml", ".github/workflows/test-coverage.yaml")
hashes <- setNames(lapply(source_paths, function(p)
  digest::digest(file = file.path(root, p), algo = "sha256")), source_paths)
metadata <- list(source_sha = expected_sha,
  baseline_source_sha = "033649b93115c6d8191e791f06b1ddf52eb82d19",
  source_hashes = hashes, status = "running", started_at = format(Sys.time(), tz = "UTC"),
  r_version = R.version.string, threads = threads)
write_manifest <- function() jsonlite::write_json(metadata,
  file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
write_manifest()
writeLines(capture.output(sessionInfo()), file.path(out, "session-before.txt"))
installed <- installed.packages()
write.csv(installed[, c("Package", "Version", "LibPath"), drop = FALSE],
          file.path(out, "installed-packages.csv"), row.names = FALSE)
description <- read.dcf(file.path(root, "DESCRIPTION"))
suggests <- trimws(strsplit(gsub("\\([^)]*\\)", "", description[1, "Suggests"]), ",")[[1]])
missing <- suggests[!vapply(suggests, requireNamespace, logical(1), quietly = TRUE)]
metadata$missing_suggests <- unname(missing)
if (length(missing)) {
  metadata$status <- "missing_dependencies"
  write_manifest()
  stop("Missing Suggests: ", paste(missing, collapse = ", "))
}
torch::torch_set_num_threads(threads)
torch::torch_set_num_interop_threads(1L)
stopifnot(isTRUE(torch::torch_is_installed()),
          identical(as.numeric(torch::torch_tensor(1)$item()), 1))
result <- tryCatch({
  cov <- covr::package_coverage(path = root, type = "tests", quiet = FALSE,
    clean = FALSE, install_path = file.path(out, "package"))
  saveRDS(cov, file.path(out, "coverage.rds"))
  covr::to_cobertura(cov, filename = file.path(out, "cobertura.xml"))
  write.csv(as.data.frame(cov), file.path(out, "line-coverage.csv"), row.names = FALSE)
  print(cov)
  metadata$xml_sha256 <- digest::digest(file = file.path(out, "cobertura.xml"), algo = "sha256")
  metadata$status <- "ok"
  TRUE
}, error = function(e) {
  metadata$status <<- "failed"
  metadata$error <<- conditionMessage(e)
  FALSE
})
metadata$finished_at <- format(Sys.time(), tz = "UTC")
write_manifest()
writeLines(capture.output(sessionInfo()), file.path(out, "session-after.txt"))
if (!result) stop("Coverage failed: ", metadata$error)
cat("COVERAGE_COLLECTION_PASS", expected_sha, "\n")
