# script/rubin_pig_prior_trial_summary.R
#
# Score the Sigma_E prior trial (script/rubin_pig_prior_trial.R): per variant and dataset, the Rubin-pooled PGLS slope
# error against complete data, the oracle's error on the same dataset (draws from the exact conditional), the pooled SE
# and coverage of rho, and per-cell 95% interval coverage of the missing c1/c2 cells from the 20 draws.
#
#   Rscript script/rubin_pig_prior_trial_summary.R <trial dir>

args <- commandArgs(trailingOnly = TRUE); dir <- args[1]
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
suppressMessages({ source(file.path(here, "campaign_gnn_off_lib.R")); source(file.path(here, "rubin_lib.R")) })
src <- readLines(file.path(here, "rubin_pig_oracle.R"))
eval(parse(text = src[grep("^oracle_sets <- function", src):(grep("^rows <- list\\(\\)", src) - 1)]))   # oracle_sets(), score()
RNGkind("L'Ecuyer-CMRG")

cell_cov <- function(sets, x) {
  mk <- as.matrix(x$mask[, c("c1", "c2")]); tr <- as.matrix(x$truth[, c("c1", "c2")])
  A <- simplify2array(lapply(sets, function(d) as.matrix(d[rownames(x$truth), c("c1", "c2")])))
  lo <- apply(A, c(1, 2), stats::quantile, 0.025); hi <- apply(A, c(1, 2), stats::quantile, 0.975)
  mean((tr >= lo & tr <= hi)[mk])
}

files <- list.files(dir, "^prior_(base|A|B|C|SEP)_n.*[.]rds$", full.names = TRUE)
oracle_cache <- list(); rows <- list()
for (f in files) {
  x <- readRDS(f); key <- sprintf("l%s_s%d", format(x$lambda), x$seed)
  eig <- pagel_eigen(x$tree, rownames(x$truth)); comp <- est_pgls_slope_fast(x$truth, x$tree, eig = eig)
  if (is.null(oracle_cache[[key]])) { set.seed(x$seed + 4242L); oracle_cache[[key]] <- score(oracle_sets(x), x, eig)[["estimate"]] - comp$estimate }
  s <- score(x$pig_post, x, eig)
  rows[[f]] <- data.frame(variant = x$variant, lambda = x$lambda, seed = x$seed,
    err = s[["estimate"]] - comp$estimate, oracle_err = oracle_cache[[key]], se = s[["se"]], covered = s[["covered"]],
    cell_cov = cell_cov(x$pig_post, x), converged = isTRUE(x$pig_diag$converged), wall_s = x$wall_s)
}
d <- do.call(rbind, rows); d$minus_oracle <- d$err - d$oracle_err
options(width = 200)
num <- vapply(d, is.numeric, logical(1)); d2 <- d; d2[num] <- lapply(d2[num], round, 4)
print(d2[order(d2$lambda, d2$seed, d2$variant), ], row.names = FALSE)
cat("\nMeans by variant and lambda:\n")
print(aggregate(cbind(err, oracle_err, minus_oracle, se, covered, cell_cov, converged, wall_s) ~ variant + lambda, d, mean), digits = 3)
saveRDS(d, file.path(dir, "prior_trial_summary.rds"))
