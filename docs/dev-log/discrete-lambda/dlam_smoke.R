# Smoke: discrete-trait accuracy with options(pigauto.discrete_lambda) off vs on.
# DGP as in the four-arm study (simplified): liabilities L = V_lambda^{1/2} Z Sigma^{1/2}, rho among traits,
# bin = 1(L > 0), cat3 = terciles of L, two continuous traits; 30% MCAR per trait.
suppressMessages(devtools::load_all(commandArgs(TRUE)[1], quiet = TRUE))
sim <- function(n, lambda, rho, seed) {
  set.seed(seed); tr <- ape::rcoal(n); tr$edge.length <- tr$edge.length / max(ape::node.depth.edgelength(tr))
  C <- ape::vcv(tr, corr = TRUE); V <- lambda * C + (1 - lambda) * diag(n)
  S <- matrix(rho, 4, 4); diag(S) <- 1
  L <- t(chol(V)) %*% matrix(rnorm(n * 4), n) %*% chol(S); rownames(L) <- tr$tip.label
  d <- data.frame(c1 = L[, 1], c2 = L[, 2],
                  bin = factor(ifelse(L[, 3] > 0, "yes", "no")),
                  cat3 = factor(cut(L[, 4], qnorm(c(0, 1/3, 2/3, 1)), labels = c("a", "b", "c"))),
                  row.names = tr$tip.label)
  truth <- d; for (v in names(d)) d[sample(n, round(0.3 * n)), v] <- NA
  list(tree = tr, d = d, truth = truth)
}
acc <- function(res, s, v) { m <- is.na(s$d[[v]]); mean(as.character(res$completed[m, v]) == as.character(s$truth[m, v])) }
modeacc <- function(s, v) { m <- is.na(s$d[[v]]); md <- names(which.max(table(s$d[[v]]))); mean(as.character(s$truth[m, v]) == md) }
out <- list()
for (lam in c(0.3, 0.7, 1)) for (seed in 1:3) {
  s <- sim(300, lam, 0.5, seed)
  r <- list()
  for (mode in c("fixed_1", "estimate")) {
    options(pigauto.discrete_lambda = mode)
    res <- suppressMessages(suppressWarnings(impute(s$d, s$tree, verbose = FALSE)))
    r[[mode]] <- c(bin = acc(res, s, "bin"), cat3 = acc(res, s, "cat3"))
  }
  out[[length(out) + 1]] <- data.frame(lambda = lam, seed = seed,
    bin_off = r$fixed_1[["bin"]], bin_on = r$estimate[["bin"]], bin_mode = modeacc(s, "bin"),
    cat_off = r$fixed_1[["cat3"]], cat_on = r$estimate[["cat3"]], cat_mode = modeacc(s, "cat3"))
}
o <- do.call(rbind, out); print(o, digits = 3, row.names = FALSE)
print(aggregate(. ~ lambda, o[, -2], mean), digits = 3)
