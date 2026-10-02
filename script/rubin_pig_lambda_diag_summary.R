# script/rubin_pig_lambda_diag_summary.R
#
# Summarise the lambda diagnostic refits (script/rubin_pig_lambda_diag.R) beside the stored campaign rows for the same
# datasets (pig_post and freqA).
#
#   Rscript script/rubin_pig_lambda_diag_summary.R <diag dir> [pool dir]

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; pool <- if (length(args) >= 2) args[2] else path.expand("~/pigauto_rubin_pool")
tag_of <- function(x) sprintf("rubin_types_mixed_BM_l%s_r%s_mcar0.3_n%d_M20_s%d", format(x$lambda), format(x$rho), x$n, x$seed)
rows <- list()
for (f in list.files(dir, "^lambda_diag_.*rds$", full.names = TRUE)) {
  x <- readRDS(f); bt <- x$block_traits; L <- x$params$lambda; colnames(L) <- bt
  SE <- x$params$Sigma_E; SP <- x$params$Sigma_P
  i1 <- match("c1", bt); i2 <- match("c2", bt)
  se_cor <- SE[i1, i2, ] / sqrt(SE[i1, i1, ] * SE[i2, i2, ]); sp_cor <- SP[i1, i2, ] / sqrt(SP[i1, i1, ] * SP[i2, i2, ])
  tot_cor <- (SP[i1, i2, ] + SE[i1, i2, ]) / sqrt((SP[i1, i1, ] + SE[i1, i1, ]) * (SP[i2, i2, ] + SE[i2, i2, ]))
  pg <- Sys.glob(file.path(pool, "pig", "*", paste0(tag_of(x), ".rds"))); fq <- Sys.glob(file.path(pool, "freq", "*", paste0(tag_of(x), ".rds")))
  est <- function(g, a) { if (!length(g)) return(c(NA, NA)); e <- readRDS(g[1])$estimands; c(e$estimate[e$arm == a & e$estimand == "slope"],
                                                                                         e$estimate[e$arm == "complete" & e$estimand == "slope"]) }
  ep <- est(pg, "pig_post"); ef <- est(fq, "freqA")
  lam_f <- if (length(fq)) mean(readRDS(fq[1])$diag$freqA$lambda_star) else NA
  rows[[f]] <- data.frame(lambda = x$lambda, rho = x$rho, seed = x$seed, converged = x$converged,
    lam_c1 = mean(L[, "c1"]), lam_c2 = mean(L[, "c2"]), lam_min_q025 = min(apply(L, 2, quantile, 0.025)),
    p_lam_gt_099 = mean(L[, c("c1", "c2")] > 0.99), lam_freqA = lam_f,
    SE_c1 = mean(SE[i1, i1, ]), truth_var_c1 = x$truth_cov[i1, i1], SE_cor = mean(se_cor), SP_cor = mean(sp_cor),
    total_cor = mean(tot_cor), truth_cor = x$truth_cor[i1, i2],
    slope_err_pig = ep[1] - ep[2], slope_err_freqA = ef[1] - ef[2])
}
d <- do.call(rbind, rows); num <- vapply(d, is.numeric, logical(1)); d[num] <- lapply(d[num], round, 3)
print(d[order(d$lambda, d$rho, d$seed), ], row.names = FALSE)
cat("\nMeans by cell:\n"); print(aggregate(. ~ lambda + rho, d[, setdiff(names(d), c("seed", "converged"))], mean), digits = 3)
