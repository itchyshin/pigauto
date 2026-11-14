test_that("one-missing-cell EM update matches a dense Gaussian Schur complement", {
  set.seed(104)
  tree <- ape::rtree(8)
  L <- cbind(x = c(NA, 2, 4, 3, 7, 5, 8, 6),
             y = c(2, 5, 3, 6, 4, 8, 7, 9))
  rownames(L) <- tree$tip.label
  R <- stats::cov2cor(ape::vcv.phylo(tree))[rownames(L), rownames(L)]
  Rinv <- solve(R)
  zero_filled <- L; zero_filled[is.na(zero_filled)] <- 0
  sigma0 <- crossprod(zero_filled, Rinv %*% zero_filled) / nrow(L)
  C <- kronecker(sigma0, R)
  v <- as.vector(L); missing <- which(is.na(v)); observed <- which(!is.na(v))
  A <- C[missing, observed, drop = FALSE] %*% solve(C[observed, observed, drop = FALSE])
  conditional_mean <- as.numeric(A %*% v[observed])
  conditional_var <- as.numeric(C[missing, missing, drop = FALSE] -
                                  A %*% C[observed, missing, drop = FALSE])
  completed <- L; completed[is.na(completed)] <- conditional_mean
  expected <- crossprod(completed, Rinv %*% completed)
  expected[1, 1] <- expected[1, 1] + Rinv[1, 1] * conditional_var
  expected <- expected / nrow(L)
  got <- tryCatch(.dcb_sigma_em(L, R, max_iter = 1L), error = function(e) NULL)
  expect_true(is.matrix(got))
  if (is.matrix(got)) {
    expect_equal(unname(got), unname(expected), tolerance = 1e-8)
    expect_true(all(is.finite(got)))
    expect_true(isSymmetric(got))
  }
})

test_that("conditional draws preserve observed cells with one missing value", {
  set.seed(104)
  tree <- ape::rtree(8)
  d <- data.frame(x = c(NA, 2, 4, 3, 7, 5, 8, 6),
                  y = c(2, 5, 3, 6, 4, 8, 7, 9), row.names = tree$tip.label)
  data <- preprocess_traits(d, tree, log_transform = FALSE)
  result <- tryCatch(draw_conditional_bm(list(data = data, tree = tree), m = 2, seed = 21),
                     error = function(e) NULL)
  expect_s3_class(result, "pigauto_draws_conditional_bm")
  if (!is.null(result)) {
    expect_true(all(is.finite(result$sigma_hat)))
    expect_true(isSymmetric(result$sigma_hat))
    expect_equal(sum(result$imputed_mask), 1)
    expect_true(all(is.finite(vapply(result$datasets, function(x) x$x[1], numeric(1)))))
    expect_identical(result$datasets[[1]]$x[-1], d$x[-1])
    expect_identical(result$datasets[[1]]$y, d$y)
    X <- data$X_scaled
    R <- stats::cov2cor(ape::vcv.phylo(tree))[rownames(X), rownames(X)]
    C <- kronecker(result$sigma_hat + diag(1e-8, ncol(X)), R)
    v <- as.vector(X); miss <- which(is.na(v)); obs <- which(!is.na(v))
    dense_mean <- C[miss, obs, drop = FALSE] %*% solve(C[obs, obs, drop = FALSE], v[obs])
    expect_equal(as.numeric(result$mu_cond[miss]), as.numeric(dense_mean), tolerance = 1e-6)
  }
})
