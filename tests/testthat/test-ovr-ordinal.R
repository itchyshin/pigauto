# Prototype (feat/discrete-lambda-ordinal): options(pigauto.ordinal_method = "cumulative")
# replaces ordinal baselines with a cumulative decomposition (K-1 binary "class >= k" fits).
# Unset must leave fit_baseline() unchanged.

sim_ordinal_cum <- function(n = 80, seed = 11) {
  set.seed(seed)
  tr <- ape::rcoal(n)
  C <- ape::vcv(tr, corr = TRUE)
  L <- t(chol(0.5 * C + 0.5 * diag(n))) %*% matrix(stats::rnorm(n * 2), n)
  rownames(L) <- tr$tip.label
  d <- data.frame(c1 = L[, 1],
                  ord = factor(cut(L[, 2], stats::qnorm(c(0, .25, .5, .75, 1)),
                                   labels = c("a", "b", "c", "d")), ordered = TRUE),
                  row.names = tr$tip.label)
  for (v in names(d)) d[sample(n, round(0.25 * n)), v] <- NA
  list(tree = tr, d = d)
}

test_that("ordinal_method unset equals 'route' in fit_baseline()", {
  skip_if_not(joint_mvn_available())
  s <- sim_ordinal_cum()
  pd <- preprocess_traits(s$d, s$tree)
  old <- options(pigauto.ordinal_method = NULL); on.exit(options(old), add = TRUE)
  a <- fit_baseline(pd, s$tree)
  options(pigauto.ordinal_method = "route")
  b <- fit_baseline(pd, s$tree)
  expect_identical(a$mu, b$mu)
  expect_identical(a$se, b$se)
})

test_that("decode_cumulative_ordinal gives valid, monotone-safe probabilities", {
  cum <- rbind(c(0.8, 0.5, 0.2),   # monotone
               c(0.3, 0.6, 0.1),   # non-monotone
               c(NA, 0.4, 0.1))    # failed first fit
  dec <- decode_cumulative_ordinal(cum)
  expect_equal(dim(dec$probs), c(3L, 4L))
  expect_equal(rowSums(dec$probs), rep(1, 3), tolerance = 1e-12)
  expect_true(all(dec$probs >= 0 & dec$probs <= 1))
  expect_true(all(dec$class %in% 0:3))
  # non-monotone row: cummin makes P(Y>=2)=P(Y>=1)=0.3, so P(Y=1)=0
  expect_equal(dec$probs[2, 2], 0)
})

test_that("cumulative ordinal path runs through impute() and fills valid levels", {
  skip_if_not(joint_mvn_available())
  s <- sim_ordinal_cum()
  old <- options(pigauto.ordinal_method = "cumulative"); on.exit(options(old), add = TRUE)
  res <- impute(s$d, s$tree, verbose = FALSE)
  miss <- is.na(s$d$ord)
  expect_true(any(miss))
  filled <- res$completed$ord[miss]
  expect_false(anyNA(filled))
  expect_true(all(as.character(filled) %in% levels(s$d$ord)))
})
