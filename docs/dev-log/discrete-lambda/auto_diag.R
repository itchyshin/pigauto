suppressMessages(devtools::load_all(".", quiet = TRUE))
sim <- function(n, lambda, seed) { set.seed(seed); tr <- ape::rcoal(n); C <- ape::vcv(tr, corr = TRUE); V <- lambda * C + (1 - lambda) * diag(n)
  L <- t(chol(V)) %*% matrix(rnorm(n * 3), n); rownames(L) <- tr$tip.label
  d <- data.frame(c1 = L[, 1], bin = factor(ifelse(L[, 2] > 0, "yes", "no")), cat3 = factor(cut(L[, 3], qnorm(c(0, 1/3, 2/3, 1)), labels = c("a","b","c"))), row.names = tr$tip.label)
  for (v in names(d)) d[sample(n, round(0.3 * n)), v] <- NA; list(d = d, tree = tr) }
for (lam in c(0.3, 1)) for (seed in 1:3) {
  s <- sim(300, lam, seed); pd <- preprocess_traits(s$d, s$tree); set.seed(seed); sp <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, val_frac = 0.5)
  out <- list()
  for (m in c("fixed_1", "estimate", "auto")) { options(pigauto.discrete_lambda = m); b <- suppressWarnings(suppressMessages(fit_baseline(pd, s$tree, splits = sp)))
    out[[m]] <- b }
  ch <- out$auto$discrete_lambda_chosen
  bcol <- which(colnames(pd$X_scaled) == "bin")
  same_as <- function(a, b) isTRUE(all.equal(a$mu[, bcol], b$mu[, bcol]))
  cat(sprintf("lambda=%s seed=%d chosen: %s | auto bin mu == fixed: %s, == est: %s\n", lam, seed, paste(names(ch), ch, collapse = ", "),
      same_as(out$auto, out$fixed_1), same_as(out$auto, out$estimate)))
}
