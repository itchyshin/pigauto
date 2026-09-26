# script/rubin_discrete_rescore.R
#
# Scores the discrete traits of result files written by rubin_cell.R --save_imp, from the saved imputations: rebuilds
# the dataset from (n, lambda, rho, seed), checks it against the stored truth and mask, and writes disc_cells,
# disc_estimands, disc_fill and disc_est_status into the rds (the same scores rubin_cell.R --discrete computes in the
# run). Stale discrete-scoring errors from the run are removed, and the rescore is recorded in x$rescored.
# A file that fails is reported and left unchanged; the other files are still scored.
# Usage: Rscript script/rubin_discrete_rescore.R FILE.rds [FILE.rds ...]
suppressPackageStartupMessages({ library(ape) })
RNGkind("L'Ecuyer-CMRG")
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "campaign_gnn_off_lib.R")); source(file.path(here, "rubin_lib.R"))
source(file.path(here, "rubin_discrete.R"))
git_hash <- tryCatch(system2("git", c("-C", here, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
                     error = function(e) NA_character_, warning = function(w) NA_character_)
if (is.na(git_hash) && file.exists(file.path(here, "DISC_COMMIT"))) git_hash <- readLines(file.path(here, "DISC_COMMIT"))[1]

rescore_file <- function(f) {
  x <- readRDS(f)
  if (is.null(x$imputations)) { message(f, ": no saved imputations, skipped"); return(invisible(FALSE)) }
  cell <- make_cell("types_mixed", x$n, x$seed, miss_frac = x$miss_frac, miss = x$miss, lambda = x$lambda,
                    rho = x$rho, thresholds = x$thresholds, driver = x$driver)
  stopifnot(identical(cell$mask, x$mask), isTRUE(all.equal(cell$truth, x$truth, tolerance = 1e-10)))
  truth <- x$truth; mask <- x$mask; eig <- pagel_eigen(cell$tree, rownames(truth))
  L5 <- cell$L[rownames(truth), 5]
  dc <- list(); de <- list(); fill <- list(); status <- list()
  for (arm in names(x$imputations)) {
    x$errors[[paste0(arm, "_disc_score")]] <- NULL; x$errors[[paste0(arm, "_disc_est")]] <- NULL
    fd <- fill_degenerate(x$imputations[[arm]], cell$df_miss, disc_traits_of(truth)); fill[[arm]] <- fd$n_filled
    r <- tryCatch(score_discrete(arm, fd$sets, truth, mask), error = function(e) e)
    if (inherits(r, "error")) x$errors[[paste0(arm, "_disc_score")]] <- conditionMessage(r) else dc[[arm]] <- r
    r <- tryCatch(score_discrete_estimand(arm, fd$sets, truth, cell$tree, eig, x$n, L5 = L5, lambda = x$lambda,
                                          rho = x$rho), error = function(e) e)
    if (inherits(r, "error")) {
      x$errors[[paste0(arm, "_disc_est")]] <- conditionMessage(r); status[[arm]] <- "error"
    } else if (is.null(r)) status[[arm]] <- "undefined" else { de[[arm]] <- r; status[[arm]] <- "scored" }
  }
  x$disc_cells <- do.call(rbind, dc); x$disc_estimands <- do.call(rbind, de); x$disc_fill <- fill
  x$disc_est_status <- unlist(status)
  x$rescored <- list(time = Sys.time(), git_hash = git_hash)
  saveRDS(x, f)
  message(sprintf("%s: scored %s", basename(f), paste(names(x$imputations), collapse = ",")))
  invisible(TRUE)
}

files <- commandArgs(trailingOnly = TRUE); n_fail <- 0L
for (f in files) {
  ok <- tryCatch(rescore_file(f), error = function(e) { message(f, ": FAILED, left unchanged: ", conditionMessage(e)); NA })
  if (is.na(ok)) n_fail <- n_fail + 1L
}
message(sprintf("rescored %d of %d files; %d failed", length(files) - n_fail, length(files), n_fail))
