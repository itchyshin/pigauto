# Move every BACE result whose fit failed into results/bace_asshipped_failed/ (kept as the record of the
# as-shipped failure rate) and write retry lists "lambda rho seed" per n for rubin_campaign_nibi.sh RETRY_FILE.
#   Rscript script/rubin_retry_prep.R <root>
root <- commandArgs(trailingOnly = TRUE)[1]
src <- file.path(root, "results", "bace"); dst <- file.path(root, "results", "bace_asshipped_failed")
dir.create(dst, showWarnings = FALSE, recursive = TRUE)
fs <- list.files(src, "^rubin_.*\\.rds$", full.names = TRUE); moved <- list()
for (f in fs) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(x) || is.null(x$errors$bace_fit)) next
  file.rename(f, file.path(dst, basename(f)))
  moved[[length(moved) + 1L]] <- data.frame(n = x$n, lambda = x$lambda, rho = x$rho, seed = x$seed, err = x$errors$bace_fit)
}
m <- do.call(rbind, moved)
if (is.null(m)) { cat("nothing to retry\n"); quit(save = "no") }
for (n in unique(m$n)) write.table(m[m$n == n, c("lambda", "rho", "seed")], file.path(root, "logs", sprintf("retry_n%d.txt", n)),
                                   row.names = FALSE, col.names = FALSE, quote = FALSE)
print(table(n = m$n, lambda = m$lambda)); print(table(trimws(m$err)))
