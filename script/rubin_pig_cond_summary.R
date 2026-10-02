# script/rubin_pig_cond_summary.R
#
# Compare each arm's imputations of the missing c1 and c2 cells (script/rubin_pig_cond_diag.R output) with the EXACT
# conditional under the true data-generating model, given the same observed block cells the arms see.
#
# True model on the block (c1, c2, logit prp, d1) = latent columns 1, 2, 4, 8 of sim_latents():
#   vec(Y) ~ N(0, Sigma_block %x% V_lambda + diag(0, 0, 0.3^2, 0) %x% I_n)
# Sigma_block: rho between c1, c2 and the prp latent, 0.35 between each of them and the driver d1, unit diagonal.
# (prp's clipping at 1e-4 is ignored.)
#
#   Rscript script/rubin_pig_cond_summary.R <diag dir>

args <- commandArgs(trailingOnly = TRUE); dir <- args[1]
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
suppressMessages(source(file.path(here, "campaign_gnn_off_lib.R")))   # phylo_corr()

exact_cond <- function(x) {
  bt <- x$block_traits; n <- x$n; K <- length(bt)
  stopifnot(identical(bt, c("c1", "c2", "prp", "d1")))
  S <- matrix(x$rho, K, K); S[4, ] <- S[, 4] <- 0.35; diag(S) <- 1
  V <- phylo_corr(x$tree, x$lambda, "BM", 2)
  sp <- rownames(x$truth); V <- V[sp, sp]
  C <- kronecker(S, V) + kronecker(diag(c(0, 0, 0.09, 0)), diag(n))
  Y <- as.matrix(x$truth); Y[, "prp"] <- stats::qlogis(Y[, "prp"])
  y <- as.vector(Y); miss <- as.vector(as.matrix(x$mask))
  o <- which(!miss); m <- which(miss)
  W <- C[m, o] %*% solve(C[o, o])                       # m x o
  list(idx = m, mean = as.numeric(W %*% y[o]), sd = sqrt(pmax(diag(C[m, m]) - rowSums(W * C[m, o]), 0)),
       truth = y[m], trait = rep(bt, each = n)[m], row = rep(seq_len(n), K)[m])
}

arm_cells <- function(sets, x, idx) {
  M <- sapply(sets, function(d) { Y <- as.matrix(d[rownames(x$truth), x$block_traits]); Y[, "prp"] <- stats::qlogis(Y[, "prp"]); as.vector(Y)[idx] })
  list(mean = rowMeans(M), sd = apply(M, 1, sd))
}

rows <- list()
for (f in list.files(dir, "^cond_diag_.*rds$", full.names = TRUE)) {
  x <- readRDS(f); ex <- exact_cond(x)
  other_obs <- function(tr) { o <- if (tr == "c1") "c2" else "c1"; !x$mask[ex$row, o] }
  for (arm in c("pig_post", "freqA")) {
    a <- arm_cells(x[[arm]], x, ex$idx)
    for (tr in c("c1", "c2")) for (oo in c(TRUE, FALSE)) {
      k <- ex$trait == tr & other_obs(tr) == oo
      if (sum(k) < 5) next
      rows[[length(rows) + 1]] <- data.frame(lambda = x$lambda, seed = x$seed, arm = arm, trait = tr,
        other_observed = oo, cells = sum(k),
        mean_diff = mean(a$mean[k] - ex$mean[k]),
        shrink = unname(coef(lm(a$mean[k] ~ ex$mean[k]))[2]),
        sd_ratio = sqrt(mean(a$sd[k]^2)) / sqrt(mean(ex$sd[k]^2)),
        rmse_arm = sqrt(mean((a$mean[k] - ex$truth[k])^2)), rmse_exact = sqrt(mean((ex$mean[k] - ex$truth[k])^2)))
    }
  }
}
d <- do.call(rbind, rows)
options(width = 200)
s <- aggregate(cbind(cells, mean_diff, shrink, sd_ratio, rmse_arm, rmse_exact) ~ lambda + arm + trait + other_observed, d, mean)
num <- vapply(s, is.numeric, logical(1)); s[num] <- lapply(s[num], round, 3)
print(s[order(s$lambda, s$trait, s$other_observed, s$arm), ], row.names = FALSE)
saveRDS(d, file.path(dir, "cond_summary.rds"))
