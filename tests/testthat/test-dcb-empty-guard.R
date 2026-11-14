dcb_empty_fixture <- function() {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  traits <- data.frame(x = rep(NA_real_, 4), y = rep(NA_real_, 4), row.names = tree$tip.label)
  list(data = preprocess_traits(traits, tree, log_transform = FALSE), tree = tree)
}

test_that("all-missing cells produce the clear diagnostic with EM covariance", {
  expect_error(draw_conditional_bm(dcb_empty_fixture(), m = 2L, sigma_method = "em"),
               "no observed cells to condition on")
})

test_that("all-missing cells produce the clear diagnostic with inhouse covariance", {
  expect_error(draw_conditional_bm(dcb_empty_fixture(), m = 2L, sigma_method = "inhouse"),
               "no observed cells to condition on")
})

test_that("trait and multi-observation admission precede the empty-cell guard", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  mixed <- data.frame(x = rep(NA_real_, 4),
                      b = factor(rep(NA_character_, 4), levels = c("a", "b")),
                      row.names = tree$tip.label)
  pd <- preprocess_traits(mixed, tree, log_transform = FALSE)
  expect_error(draw_conditional_bm(list(data = pd, tree = tree), m = 2L), "continuous traits only")
  multi <- data.frame(species = rep(tree$tip.label, each = 2L), x = rep(NA_real_, 8))
  pd <- preprocess_traits(multi, tree, species_col = "species", log_transform = FALSE)
  expect_error(draw_conditional_bm(list(data = pd, tree = tree), m = 2L), "single-observation-per-species")
  wrong_tree <- ape::read.tree(text = "(((a:1,b:1):1,c:2):1,(d:1,e:2):1);")
  expect_error(draw_conditional_bm(list(data = dcb_empty_fixture()$data, tree = wrong_tree), m = 2L),
               "tree tip count does not match")
})

test_that("partially observed continuous data retain real conditional draws", {
  tree <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
  traits <- data.frame(x = c(1, NA, 4, 3, NA, 5, 8, 6),
                       y = c(2, 5, 3, 6, 4, 8, 7, 9), row.names = tree$tip.label)
  pd <- preprocess_traits(traits, tree, log_transform = FALSE)
  X <- pd$X_scaled
  R <- cov2cor(ape::vcv(tree))[rownames(X), rownames(X)]
  for (method in c("em", "inhouse")) {
    out <- draw_conditional_bm(list(data = pd, tree = tree), m = 2L, seed = 9L,
                               sigma_method = method)
    C <- kronecker(out$sigma_hat + diag(1e-8, ncol(X)), R)
    v <- as.vector(X); missing <- which(is.na(v)); observed <- which(!is.na(v))
    expected <- C[missing, observed, drop = FALSE] %*%
      solve(C[observed, observed, drop = FALSE], v[observed])
    expect_equal(as.vector(out$mu_cond)[missing], as.vector(expected), tolerance = 1e-6,
                 info = method)
    expect_identical(out$imputed_mask, is.na(X), info = method)
    for (draw in out$datasets) {
      expect_true(all(is.finite(as.matrix(draw))), info = method)
      expect_identical(draw$x[!is.na(traits$x)], traits$x[!is.na(traits$x)], info = method)
      expect_identical(draw$y, traits$y, info = method)
    }
  }
})
