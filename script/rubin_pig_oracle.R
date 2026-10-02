# script/rubin_pig_oracle.R
#
# Oracle MI on the conditional-diagnostic datasets (script/rubin_pig_cond_diag.R): M = 20 JOINT draws of all missing
# block cells from the exact conditional under the true model, scored with the campaign's own slope estimator and
# Rubin pooling (est_pgls_slope_fast, rubin_pool), beside pig_post and freqA on the same datasets. If the oracle shows
# the same slope error as pig_post, the error is a property of MI here, not of pigauto's sampler.
#
#   Rscript script/rubin_pig_oracle.R <diag dir>

args <- commandArgs(trailingOnly = TRUE); dir <- args[1]
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
suppressMessages({ source(file.path(here, "campaign_gnn_off_lib.R")); source(file.path(here, "rubin_lib.R")) })
RNGkind("L'Ecuyer-CMRG")

oracle_sets <- function(x, M = 20L) {
  bt <- x$block_traits; n <- x$n; K <- length(bt)
  S <- matrix(x$rho, K, K); S[4, ] <- S[, 4] <- 0.35; diag(S) <- 1
  sp <- rownames(x$truth); V <- phylo_corr(x$tree, x$lambda, "BM", 2)[sp, sp]
  C <- kronecker(S, V) + kronecker(diag(c(0, 0, 0.09, 0)), diag(n))
  Y <- as.matrix(x$truth); Y[, "prp"] <- stats::qlogis(Y[, "prp"]); y <- as.vector(Y)
  miss <- as.vector(as.matrix(x$mask)); o <- which(!miss); m <- which(miss)
  W <- C[m, o] %*% solve(C[o, o])
  mu <- as.numeric(W %*% y[o]); Cm <- C[m, m] - W %*% C[o, m]; Cm <- (Cm + t(Cm)) / 2
  R <- chol(Cm + diag(1e-12, length(m)))
  lapply(seq_len(M), function(i) {
    yi <- y; yi[m] <- mu + as.numeric(crossprod(R, stats::rnorm(length(m))))
    d <- as.data.frame(matrix(yi, n, K, dimnames = list(sp, bt))); d$prp <- stats::plogis(d$prp); d
  })
}

score <- function(sets, x, eig) {
  sl <- lapply(sets, function(d) est_pgls_slope_fast(d, x$tree, eig = eig))
  ok <- vapply(sl, function(s) is.finite(s$estimate) && is.finite(s$variance), logical(1))
  ps <- rubin_pool(vapply(sl[ok], `[[`, numeric(1), "estimate"), vapply(sl[ok], `[[`, numeric(1), "variance"), df_com = x$n - 2)
  c(estimate = ps$estimate, se = ps$se, covered = as.numeric(ps$lower <= x$rho & x$rho <= ps$upper),
    lambda_hat = mean(vapply(sl[ok], `[[`, numeric(1), "lambda_hat")))
}

rows <- list()
for (f in list.files(dir, "^cond_diag_.*rds$", full.names = TRUE)) {
  x <- readRDS(f); set.seed(x$seed + 4242L)
  eig <- pagel_eigen(x$tree, rownames(x$truth))
  comp <- est_pgls_slope_fast(x$truth, x$tree, eig = eig)
  for (arm in c("oracle", "pig_post", "freqA")) {
    sets <- if (arm == "oracle") oracle_sets(x) else x[[arm]]
    s <- score(sets, x, eig)
    rows[[length(rows) + 1]] <- data.frame(lambda = x$lambda, seed = x$seed, arm = arm, estimate = s[["estimate"]],
      err_vs_complete = s[["estimate"]] - comp$estimate, se = s[["se"]], se_complete = sqrt(comp$variance),
      covered = s[["covered"]], lambda_hat = s[["lambda_hat"]], lambda_hat_complete = comp$lambda_hat)
  }
}
d <- do.call(rbind, rows); options(width = 200)
num <- vapply(d, is.numeric, logical(1)); d2 <- d; d2[num] <- lapply(d2[num], round, 4)
print(d2[order(d2$lambda, d2$seed, d2$arm), ], row.names = FALSE)
cat("\nMeans:\n"); print(aggregate(cbind(err_vs_complete, se, se_complete, covered, lambda_hat) ~ lambda + arm, d, mean), digits = 3)
saveRDS(d, file.path(dir, "oracle_summary.rds"))
