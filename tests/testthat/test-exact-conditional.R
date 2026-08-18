# Exact matrix-normal conditional (predict_method = "exact").
# Design + gates: docs/dev-log/2026-08-17-exact-conditional-design.md
#
# The load-bearing test is the first one: the sparse precision-form result
# must equal a dense conditional computed directly from Sigma %x% R. If that
# holds, the method is right; everything else is plumbing.

# Dense reference: conditional of the missing cells of L under
# vec(L) ~ MVN(0, Sigma %x% M). Column-major vec == trait-major cell order.
.dense_ref <- function(L, Sigma, M) {
  V <- kronecker(Sigma, M)
  v <- as.numeric(L)
  oi <- which(!is.na(v)); mi <- which(is.na(v))
  Voo <- V[oi, oi, drop = FALSE]
  list(mu  = as.numeric(V[mi, oi, drop = FALSE] %*% solve(Voo, v[oi])),
       var = diag(V[mi, mi, drop = FALSE] -
                    V[mi, oi, drop = FALSE] %*% solve(Voo, V[oi, mi, drop = FALSE])),
       mi  = mi)
}

.toy <- function(n = 25L, K = 3L, seed = 11L, miss = 0.35) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  A <- ape::vcv(tree); R <- stats::cov2cor(A)
  Sigma <- matrix(c(1, .6, -.3, .6, 1, .2, -.3, .2, 1), K, K)
  L <- t(chol(R)) %*% matrix(stats::rnorm(n * K), n, K) %*% chol(Sigma)
  rownames(L) <- tree$tip.label
  set.seed(seed + 1L)
  L[matrix(stats::runif(n * K) < miss, n, K)] <- NA
  list(tree = tree, A = A, R = R, Sigma = Sigma, L = L)
}

test_that("[exact] matches a dense Sigma-kron-R conditional to machine precision", {
  skip_if_not_installed("Matrix")
  d <- .toy()
  H <- build_henderson_S_inv(d$tree)

  # cor_scale = FALSE works on A = vcv(tree); TRUE works on R = cov2cor(A).
  # Comparing against the wrong scale makes the MEAN still match (means are
  # scale-invariant) while variances differ by the tree height -- so both
  # scales are checked explicitly here.
  for (cs in c(FALSE, TRUE)) {
    ref <- .dense_ref(d$L, d$Sigma, if (cs) d$R else d$A)
    ec <- exact_conditional_mvn(d$L, d$Sigma, H, cor_scale = cs)
    expect_false(is.null(ec))
    expect_equal(as.numeric(ec$mu)[ref$mi],  ref$mu,  tolerance = 1e-6)
    expect_equal(as.numeric(ec$var)[ref$mi], ref$var, tolerance = 1e-6)
  }
})

test_that("[exact] observed cells are echoed with zero variance", {
  skip_if_not_installed("Matrix")
  d <- .toy()
  H <- build_henderson_S_inv(d$tree)
  ec <- exact_conditional_mvn(d$L, d$Sigma, H, cor_scale = TRUE)
  obs <- !is.na(d$L)
  expect_equal(ec$mu[obs], d$L[obs])
  expect_true(all(ec$var[obs] == 0))
  expect_true(all(ec$var[!obs] > 0))
  expect_true(all(is.finite(ec$mu)))
})

test_that("[exact] refuses oversized problems instead of hanging", {
  skip_if_not_installed("Matrix")
  d <- .toy()
  H <- build_henderson_S_inv(d$tree)
  expect_null(exact_conditional_mvn(d$L, d$Sigma, H, max_cells = 5L))
})

test_that("[exact] blocked variance solve is independent of block size", {
  skip_if_not_installed("Matrix")
  d <- .toy()
  H <- build_henderson_S_inv(d$tree)
  a <- exact_conditional_mvn(d$L, d$Sigma, H, cor_scale = TRUE, block = 256L)
  b <- exact_conditional_mvn(d$L, d$Sigma, H, cor_scale = TRUE, block = 3L)
  expect_equal(a$var, b$var)
  expect_equal(a$mu, b$mu)
})

test_that("[exact] fit_baseline default is unchanged and 'exact' differs", {
  skip_if_not_installed("Matrix")
  set.seed(1)
  tree <- ape::rcoal(30)
  df <- data.frame(a = stats::rnorm(30), b = stats::rnorm(30),
                   row.names = tree$tip.label)
  df$a[1:5] <- NA; df$b[6:10] <- NA
  pd <- preprocess_traits(df, tree)
  b_def <- fit_baseline(pd, tree)
  b_pc  <- fit_baseline(pd, tree, predict_method = "per_column")
  b_ex  <- fit_baseline(pd, tree, predict_method = "exact")
  expect_equal(b_def$mu, b_pc$mu)                 # default == per_column
  expect_false(isTRUE(all.equal(b_def$mu, b_ex$mu)))
  expect_true(all(is.finite(b_ex$mu)))
  expect_equal(dim(b_ex$mu), dim(b_def$mu))
})

test_that("[exact] threads through impute() and returns a usable result", {
  skip_if_not_installed("Matrix")
  set.seed(2)
  tree <- ape::rcoal(30)
  df <- data.frame(a = stats::rnorm(30), b = stats::rnorm(30),
                   row.names = tree$tip.label)
  df$a[1:4] <- NA; df$b[5:8] <- NA
  res <- suppressWarnings(
    impute(df, tree, predict_method = "exact", epochs = 5L, verbose = FALSE,
           seed = 2))
  expect_s3_class(res, "pigauto_result")
  expect_false(anyNA(res$completed))
})

test_that("[exact] predict_method reaches the solver on MIXED-TYPE data too", {
  # Regression guard for the void-G4 bug: predict_method reached
  # fit_joint_threshold_baseline() but was never forwarded to
  # fit_joint_solver(), so on any dataset with binary/ordinal traits the
  # option was a SILENT NO-OP. The continuous-only fixture above cannot
  # catch this -- it takes the continuous-joint path. This one takes the
  # threshold-joint path, which is what real mixed-type data uses.
  # Fourth instance of the silent-parameter-drop pattern in this call
  # chain; see docs/dev-log/2026-08-17-exact-conditional-results.md.
  #
  # Asserted on OBSERVED BEHAVIOUR, not on a trace: tracing an internal
  # function is not reliably intercepted after devtools::load_all(), which
  # made an earlier version of this test fail while the feature worked.
  skip_if_not_installed("Matrix")
  set.seed(4)
  n <- 40L
  tree <- ape::rcoal(n)
  df <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    b  = factor(sample(c("no", "yes"), n, TRUE), levels = c("no", "yes")),
    row.names = tree$tip.label
  )
  df$x1[1:6] <- NA; df$x2[7:12] <- NA; df$b[13:18] <- NA
  pd <- preprocess_traits(df, tree)

  b_ex <- fit_baseline(pd, tree, predict_method = "exact")
  b_pc <- fit_baseline(pd, tree, predict_method = "per_column")

  # If predict_method were dropped again, these would be identical.
  expect_false(isTRUE(all.equal(b_pc$mu, b_ex$mu)))
  expect_true(all(is.finite(b_ex$mu)))
  expect_equal(dim(b_ex$mu), dim(b_pc$mu))
})
