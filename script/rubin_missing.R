# script/rubin_missing.R
#
# Completeness check for the pooled rubin-freq-bace campaign: every expected (n, lambda, rho, seed) per arm set
# against the rds files present (results plus the failure folders). A missing file means the task never wrote one
# (a segmentation fault or timeout kills the task, and with BLOCK > 1 also the seeds after it in that block).
#   Rscript script/rubin_missing.R [POOL_DIR]
# Writes <POOL>/agg/missing_<set>_n<n>.txt as "lambda rho seed" lines (RETRY_FILE format for rubin_campaign_nibi.sh).
a <- commandArgs(trailingOnly = TRUE)
pool <- if (length(a)) a[1] else file.path(Sys.getenv("HOME"), "pigauto_rubin_pool")
dir.create(file.path(pool, "agg"), showWarnings = FALSE)
seeds <- list(bace = list(`100` = 1:200, `300` = 1:200, `1000` = 1:100), freq = list(`100` = 1:200, `300` = 1:200, `1000` = 1:200))
dirs <- list(bace = c("bace", "bace_asshipped_failed", "bace_failed_later"), freq = "freq")
tag <- function(n, l, r, s) sprintf("rubin_types_mixed_BM_l%s_r%s_mcar0.3_n%d_M20_s%d.rds", format(l), format(r), n, s)
for (set in names(seeds)) {
  have <- unique(basename(unlist(lapply(dirs[[set]], function(d) list.files(file.path(pool, d), "\\.rds$", recursive = TRUE)))))
  have <- sub("_dup[0-9]+\\.rds$", ".rds", have)
  if (set == "bace") have <- intersect(have, basename(list.files(file.path(pool, "bace"), "\\.rds$", recursive = TRUE)))  # results of record
  for (n in names(seeds[[set]])) {
    g <- expand.grid(lambda = c(0.3, 0.7, 1), rho = c(0, 0.5), seed = seeds[[set]][[n]])
    miss <- g[!(tag(as.integer(n), g$lambda, g$rho, g$seed) %in% have), ]
    cat(sprintf("%s n=%s: expected %d, present %d, missing %d\n", set, n, nrow(g), nrow(g) - nrow(miss), nrow(miss)))
    if (nrow(miss)) write.table(miss, file.path(pool, "agg", sprintf("missing_%s_n%s.txt", set, n)), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
