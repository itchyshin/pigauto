# Move every BACE result whose fit failed into a failure folder (kept as the record) and write retry lists
# "lambda rho seed" per n for rubin_campaign_nibi.sh RETRY_FILE. Pass 1 used results/bace_asshipped_failed/
# (first-pass failures, as shipped). Later passes use results/bace_failed_pass<k>/ and NEVER overwrite: on
# 2026-09-24 pass 2 on fir renamed onto 6 existing pass-1 records (fits that failed again after the fix).
#   Rscript script/rubin_retry_prep.R <root> [dest_subfolder]
a <- commandArgs(trailingOnly = TRUE); root <- a[1]
src <- file.path(root, "results", "bace"); dst <- file.path(root, "results", if (length(a) >= 2) a[2] else "bace_asshipped_failed")
dir.create(dst, showWarnings = FALSE, recursive = TRUE)
fs <- list.files(src, "^rubin_.*\\.rds$", full.names = TRUE); moved <- list()
for (f in fs) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(x) || is.null(x$errors$bace_fit)) next
  to <- file.path(dst, basename(f)); k <- 1L
  while (file.exists(to)) { to <- file.path(dst, sub("\\.rds$", sprintf("_dup%d.rds", k), basename(f))); k <- k + 1L }
  file.rename(f, to)
  moved[[length(moved) + 1L]] <- data.frame(n = x$n, lambda = x$lambda, rho = x$rho, seed = x$seed, err = x$errors$bace_fit)
}
m <- do.call(rbind, moved)
if (is.null(m)) { cat("nothing to retry\n"); quit(save = "no") }
for (n in unique(m$n)) write.table(m[m$n == n, c("lambda", "rho", "seed")], file.path(root, "logs", sprintf("retry_n%d.txt", n)),
                                   row.names = FALSE, col.names = FALSE, quote = FALSE)
print(table(n = m$n, lambda = m$lambda)); print(table(trimws(m$err)))
