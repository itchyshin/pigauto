# Follow-up diagnostic (2026-09-23): reproduces every number the MI memo quotes
# beyond script 13: true-scale oracle (joint and x-only), per_column vs exact
# centre, oracle proper-MI GLS slopes (10 draws), and independent proper-sized
# noise around the oracle centre. Same tree and seed as the memo.
# Run from the repo root: OPENBLAS_NUM_THREADS=1 Rscript script/mondrian_confirmation/13b_mi_centre_and_scale_diag.R
suppressMessages({devtools::load_all(quiet = TRUE); library(ape)})
set.seed(20260823); n <- 400
tree <- rtree(n); V <- vcv(tree); V <- V / max(V); L <- chol(V + 1e-8 * diag(n))
Sig <- matrix(0.7, 2, 2); diag(Sig) <- 1
Z <- t(L) %*% matrix(rnorm(n * 2), n, 2) %*% chol(Sig)
truth <- data.frame(row.names = tree$tip.label, x = Z[, 1], y = Z[, 2])
miss <- sample(n, 120); df <- truth; df$x[miss] <- NA
C <- kronecker(Sig, V); io <- c(setdiff(1:n, miss), n + 1:n); obs <- c(truth$x[-miss], truth$y)
mu_c <- as.numeric(C[miss, io] %*% solve(C[io, io], obs)); Vc <- C[miss, miss] - C[miss, io] %*% solve(C[io, io], C[io, miss])
Cx <- V; io2 <- setdiff(1:n, miss)
mu_x_only <- as.numeric(Cx[miss, io2] %*% solve(Cx[io2, io2], truth$x[-miss]))
cat("oracle x-only (ignores y) residual SD:", round(sd(mu_x_only - truth$x[miss]), 3), "\n")
cat("oracle (true scale) joint residual SD:", round(sd(mu_c - truth$x[miss]), 3), " mean cond SD:", round(mean(sqrt(diag(Vc))), 3), "\n")
for (pm in c("per_column", "exact")) {
  r <- impute(df, tree, epochs = 150, verbose = FALSE, seed = 1, predict_method = pm)
  xi <- r$completed$x[miss]
  cat(pm, ": point residual SD", round(sd(xi - truth$x[miss]), 3), " cor with oracle joint mean", round(cor(xi, mu_c), 3), " cor with x-only mean", round(cor(xi, mu_x_only), 3), "\n")
}
Lc <- t(chol(Vc + 1e-9 * diag(length(miss))))
gls_slope <- function(d) { d$species <- rownames(d); coef(nlme::gls(y ~ x, correlation = corBrownian(phy = tree, form = ~species), data = d, method = "ML"))[2] }
set.seed(7)
or <- sapply(1:10, function(k) { d <- truth; d$x[miss] <- mu_c + as.numeric(Lc %*% rnorm(length(miss))); gls_slope(d) })
cat("oracle proper MI (true scale) GLS slopes:", round(or, 3), " mean", round(mean(or), 3), "\n")
cat("GLS slope on truth:", round(gls_slope(truth), 3), "\n")
iid <- sapply(1:10, function(k) { d <- truth; d$x[miss] <- mu_c + rnorm(length(miss), 0, sqrt(diag(Vc))); gls_slope(d) })
cat("oracle mean + INDEPENDENT noise of the proper marginal SD:", round(mean(iid), 3), "\n")
