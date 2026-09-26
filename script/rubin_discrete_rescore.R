# script/rubin_discrete_rescore.R
#
# Scores the discrete traits of result files written by rubin_cell.R --save_imp, from the saved imputations: rebuilds
# the dataset from (n, lambda, rho, seed), checks it against the stored truth and mask, and writes disc_cells,
# disc_estimands and disc_fill into the rds (the same scores rubin_cell.R --discrete computes in the run).
# Usage: Rscript script/rubin_discrete_rescore.R FILE.rds [FILE.rds ...]
suppressPackageStartupMessages({ library(ape) })
RNGkind("L'Ecuyer-CMRG")
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "campaign_gnn_off_lib.R")); source(file.path(here, "rubin_lib.R"))
source(file.path(here, "rubin_discrete.R"))

for (f in commandArgs(trailingOnly = TRUE)) {
  x <- readRDS(f)
  if (is.null(x$imputations)) { message(f, ": no saved imputations, skipped"); next }
  cell <- make_cell("types_mixed", x$n, x$seed, miss_frac = x$miss_frac, miss = x$miss, lambda = x$lambda,
                    rho = x$rho, thresholds = x$thresholds, driver = x$driver)
  stopifnot(identical(cell$mask, x$mask), isTRUE(all.equal(cell$truth, x$truth, tolerance = 1e-10)))
  truth <- x$truth; mask <- x$mask; eig <- pagel_eigen(cell$tree, rownames(truth))
  dc <- list(); de <- list(); fill <- list()
  for (arm in names(x$imputations)) {
    fd <- fill_degenerate(x$imputations[[arm]], cell$df_miss, disc_traits_of(truth))
    dc[[arm]] <- score_discrete(arm, fd$sets, truth, mask)
    de[[arm]] <- score_discrete_estimand(arm, fd$sets, truth, cell$tree, eig, x$n)
    fill[[arm]] <- fd$n_filled
  }
  x$disc_cells <- do.call(rbind, dc); x$disc_estimands <- do.call(rbind, de); x$disc_fill <- fill
  saveRDS(x, f)
  message(sprintf("%s: scored %s", basename(f), paste(names(x$imputations), collapse = ",")))
}
