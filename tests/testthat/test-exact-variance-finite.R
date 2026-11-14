exact_variance_fixture <- function(height) {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  tree$edge.length <- tree$edge.length * height
  L <- cbind(x = c(1, NA, -1, .5), y = c(NA, .2, .3, .4))
  rownames(L) <- tree$tip.label
  list(L = L, henderson = build_henderson_S_inv(tree), Sigma = diag(1e100, 2L))
}

expect_complete_exact_variance <- function(out, L) {
  expect_named(out, c("mu", "var"))
  expect_identical(dim(out$mu), dim(L))
  expect_identical(dim(out$var), dim(L))
  expect_identical(dimnames(out$mu), dimnames(L))
  expect_identical(dimnames(out$var), dimnames(L))
  expect_true(all(is.finite(out$mu)))
  expect_true(all(is.finite(out$var)))
  observed <- !is.na(L)
  expect_identical(out$mu[observed], L[observed])
  expect_equal(out$var[observed], rep(0, sum(observed)))
  missing <- which(is.na(L), arr.ind = TRUE)
  expect_equal(unname(out$mu[missing]), c(.5, .1), tolerance = 1e-8)
  expect_equal(unname(out$var[missing] / 1e100), c(.75, .75), tolerance = 1e-8)
}

test_that("conditional means and variances stay correct at supported branch scales", {
  for (height in c(1, 1e180)) {
    f <- exact_variance_fixture(height)
    out <- exact_conditional_mvn(f$L, f$Sigma, f$henderson, eps = 0)
    expect_false(is.null(out), info = as.character(height))
    expect_complete_exact_variance(out, f$L)
  }
})

test_that("extreme scaling returns NULL or a complete finite conditional", {
  f <- exact_variance_fixture(1e220)
  out <- exact_conditional_mvn(f$L, f$Sigma, f$henderson, eps = 0)
  expect_true(is.null(out) || is.list(out))
  if (!is.null(out)) expect_complete_exact_variance(out, f$L)
})
