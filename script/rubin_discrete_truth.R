# script/rubin_discrete_truth.R
#
# Truth for the discrete downstream estimand of the Rubin study: the PGLS slope of c1 on bin (bin = "yes" coded 1).
# It is not rho: bin thresholds a liability correlated rho with c1, and the GLS slope on a thresholded predictor
# depends on n and lambda. The truth is the mean complete-data estimate over R fresh datasets per (n, lambda, rho),
# built by the campaign's make_cell() from seeds 1e6 + 1..R (disjoint from the campaign's seeds 1..200).
# At rho = 0 the truth is 0 by symmetry; it is estimated anyway as a check. Datasets whose bin has one class in the
# complete data are left out, as in the campaign scoring.
#
# Usage: Rscript script/rubin_discrete_truth.R --n 100,300 --lambda 0.3,0.7,1 --rho 0,0.5 --R 2000 --cores 4 \
#          --out script/rubin_study/data/discrete/truth_slope_c1_bin.csv
suppressPackageStartupMessages({ library(ape) })
RNGkind("L'Ecuyer-CMRG")
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1L] }
num <- function(x) as.numeric(strsplit(x, ",")[[1]])
ns <- num(get_arg("--n", "100,300")); lams <- num(get_arg("--lambda", "0.3,0.7,1")); rhos <- num(get_arg("--rho", "0,0.5"))
R <- as.integer(get_arg("--R", "2000")); cores <- as.integer(get_arg("--cores", "4"))
out <- get_arg("--out", file.path("script", "rubin_study", "data", "discrete", "truth_slope_c1_bin.csv"))

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "campaign_gnn_off_lib.R")); source(file.path(here, "rubin_lib.R"))
source(file.path(here, "rubin_discrete.R"))

grid <- expand.grid(n = ns, lambda = lams, rho = rhos)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
  p <- grid[g, ]; t0 <- proc.time()[["elapsed"]]
  est <- unlist(parallel::mclapply(seq_len(R), function(i) {
    cell <- make_cell("types_mixed", p$n, 1e6L + i, lambda = p$lambda, rho = p$rho, thresholds = "fixed", driver = TRUE)
    if (nlevels(droplevels(cell$truth$bin)) < 2L) return(NA_real_)
    est_pgls_xy(cell$truth, cell$tree, y = "c1", x = "bin")$estimate
  }, mc.cores = cores))
  ok <- is.finite(est)
  r <- data.frame(p, R = R, R_used = sum(ok), truth = mean(est[ok]), mc_se = stats::sd(est[ok]) / sqrt(sum(ok)),
                  sd_complete = stats::sd(est[ok]), liability_ref = 2 * p$rho * sqrt(2 / pi),
                  wall_s = proc.time()[["elapsed"]] - t0)
  print(r, digits = 4); r
}))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, out, row.names = FALSE)
