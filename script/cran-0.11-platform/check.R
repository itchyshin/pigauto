# Check the original compressed archive; do not invoke R CMD build.
dir.create("platform-evidence", showWarnings = FALSE)
Sys.unsetenv("NOT_CRAN")
Sys.setenv(`_R_CHECK_FORCE_SUGGESTS_` = "true",
           `_R_CHECK_TIMINGS_` = "0",
           `_R_CHECK_ALWAYS_LOG_VIGNETTE_OUTPUT_` = "true")
stopifnot(Sys.getenv("NOT_CRAN", unset = "UNSET") == "UNSET")
desc <- read.dcf("frozen-deps/DESCRIPTION")
fields <- intersect(c("Depends", "Imports", "LinkingTo", "Suggests"), colnames(desc))
required <- trimws(unlist(strsplit(paste(desc[1, fields], collapse = ","), ",")))
required <- unique(sub("[[:space:]]*\\(.*", "", required))
required <- setdiff(required, c("R", ""))
available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
write.csv(data.frame(package = required, available = available),
          "platform-evidence/declared-dependencies.csv", row.names = FALSE)
if (!all(available)) stop("Missing declared dependencies: ", paste(required[!available], collapse = ", "))
stopifnot(nzchar(Sys.which("pandoc")), nzchar(Sys.which("pdflatex")))
stopifnot(torch::torch_tensor(1)$item() == 1)
write.csv(as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")]),
          "platform-evidence/installed-packages.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "platform-evidence/session-info.txt")
environment <- Sys.getenv()
keep <- grepl("^(_R_CHECK_|NOT_CRAN$|CI$|GITHUB_(SHA|RUN_ID|RUN_ATTEMPT|HEAD_REF)$|R_(VERSION|LIBS|HOME)|OMP_NUM_THREADS$|OPENBLAS_NUM_THREADS$)", names(environment))
writeLines(c("NOT_CRAN=UNSET", paste(names(environment)[keep], environment[keep], sep = "=")),
           "platform-evidence/check-environment.txt")
archive <- normalizePath("frozen-artifact/pigauto_0.11.0.tar.gz", mustWork = TRUE)
r <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
args <- c("CMD", "check", "--as-cran", "--run-donttest", "--timings", shQuote(archive))
writeLines(c(paste(shQuote(r), paste(args, collapse = " ")),
             "No R CMD build: checking the downloaded gzip archive directly."),
           "platform-evidence/check-command.txt")
dir.create("frozen-check", showWarnings = FALSE)
started <- Sys.time()
console <- normalizePath("platform-evidence", mustWork = TRUE)
console <- file.path(console, "check-console.log")
setwd("frozen-check")
status <- system2(r, args, stdout = console, stderr = console)
setwd("..")
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
output <- readLines(console, warn = FALSE)
cat(paste(output, collapse = "\n"), "\n")
writeLines(c(paste0("exit=", status), paste0("elapsed_seconds=", elapsed)),
           "platform-evidence/check-exit.txt")
check_log <- "frozen-check/pigauto.Rcheck/00check.log"
if (!file.exists(check_log)) stop("No completed R CMD check log")
log <- readLines(check_log, warn = FALSE)
if (!any(grepl("^Status:", log))) stop("R CMD check did not finish")
if (status != 0L) stop("R CMD check failed, exit ", status)
if (!any(grepl("^Status: OK$", log))) {
  stop("Check has NOTE/WARNING findings requiring adjudication; logs retained")
}
cat("DIRECT_FROZEN_ARCHIVE_CHECK_OK elapsed_seconds=", elapsed, "\n", sep = "")
