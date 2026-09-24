# Internals of the posterior sampler behind
# multi_impute(draws_method = "posterior") (R/mi_posterior.R).
# Items refer to docs/dev-log/mi-posterior/design.md section 5a.

mip_sim <- function(n = 30L, seed = 1L, ultrametric = FALSE,
                    SP = matrix(c(0.6, 0.3, 0.3, 0.5), 2),
                    SE = matrix(c(0.3, 0.1, 0.1, 0.4), 2),
                    mu = c(0.2, -0.1)) {
  set.seed(seed)
  tree <- if (ultrametric) ape::rcoal(n) else ape::rtree(n)
  R <- stats::cov2cor(ape::vcv(tree))[tree$tip.label, tree$tip.label]
  C <- kronecker(SP, R) + kronecker(SE, diag(n))
  y <- as.vector(t(chol(C)) %*% stats::rnorm(2L * n)) + rep(mu, each = n)
  Y <- matrix(y, n, 2L, dimnames = list(tree$tip.label, c("t1", "t2")))
  list(Y = Y, tree = tree, R = R, C = C, SP = SP, SE = SE, mu = mu)
}

dense_cond <- function(C, Y, mu) {
  v <- as.vector(Y)
  mi <- which(is.na(v)); o <- which(!is.na(v))
  mv <- rep(mu, each = nrow(Y))
  list(mean = as.numeric(mv[mi] + C[mi, o] %*% solve(C[o, o], v[o] - mv[o])),
       cov = C[mi, mi] - C[mi, o] %*% solve(C[o, o], C[o, mi]))
}

# ---- Item 1: Qc tip-block Schur complement equals cov2cor(vcv(tree)) -------
test_that("Qc = D Q D has tip-block marginal R on ultrametric and non-ultrametric trees", {
  set.seed(11)
  for (tree in list(ape::rcoal(25), ape::rtree(25))) {
    qc <- .mip_build_Qc(tree)
    Qd <- as.matrix(qc$Qc)
    n <- qc$n
    I <- (n + 1L):qc$N
    S <- Qd[1:n, 1:n] - Qd[1:n, I] %*% solve(Qd[I, I], Qd[I, 1:n])
    R <- stats::cov2cor(ape::vcv(tree))[tree$tip.label, tree$tip.label]
    expect_lt(max(abs(solve(S) - R)), 1e-8)
    expect_lt(max(abs(solve(Qd)[1:n, 1:n] - R)), 1e-8)
    expect_identical(qc$tip_labels, tree$tip.label)
  }
  # Non-ultrametric: tip depths really differ, so the D scaling matters.
  tr <- ape::rtree(25)
  expect_gt(stats::sd(diag(ape::vcv(tr))), 0.01)
})

# ---- Item 2: xi' Qc xi against the dense reference ---------------------------
test_that("the Sigma_W sufficient statistic xi' Qc xi matches the dense form", {
  set.seed(2)
  tree <- ape::rtree(10)
  qc <- .mip_build_Qc(tree)
  xi <- matrix(stats::rnorm(qc$N * 2L), qc$N, 2L)
  sparse <- as.matrix(crossprod(xi, qc$Qc %*% xi))
  A <- solve(as.matrix(qc$Qc))
  dense <- t(xi) %*% solve(A) %*% xi
  expect_equal(sparse, dense, tolerance = 1e-8, ignore_attr = TRUE)
})

# ---- Item 3: (a, mu) precision is positive definite with mu included ---------
test_that("the (a, mu) precision is PD with a flat mu, including all-missing tips", {
  d <- mip_sim(n = 20L, seed = 3L)
  Y <- d$Y
  Y[1:3, ] <- NA                   # tips with nothing observed
  Y[4:8, 1] <- NA
  prob <- .mip_problem(Y, d$tree)
  tpl <- .mip_template(prob, include_mu = TRUE)
  tpl <- .mip_refactor(tpl, d$SP, .mip_lik_blocks(prob, d$SE))
  P <- as.matrix(tpl$P)
  expect_equal(nrow(P), prob$N * 2L + 2L)
  expect_gt(min(eigen(P, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_s4_class(tpl$L, "CHMfactor")
  # The pattern is fixed: a second refactor with other values reuses it.
  tpl2 <- .mip_refactor(tpl, 2 * d$SP, .mip_lik_blocks(prob, 0.5 * d$SE))
  expect_identical(tpl2$P@i, tpl$P@i)
  expect_identical(tpl2$P@p, tpl$P@p)
})

# ---- Item 4: PX mapping, hand-computed K = 2 --------------------------------
test_that("PX mapping Sigma_P = diag(alpha) Sigma_W diag(alpha), a = xi diag(alpha)", {
  alpha <- c(2, -0.5)
  SW <- matrix(c(1, 0.3, 0.3, 2), 2)
  expect_equal(.mip_px_sigma(alpha, SW),
               matrix(c(4, -0.3, -0.3, 0.5), 2))
  xi <- rbind(c(1, 2), c(-3, 4))
  expect_equal(.mip_px_effects(xi, alpha), rbind(c(2, -1), c(-6, -2)))
})

# ---- Item 5: reduction to #187's single-lambda Kronecker model --------------
test_that("Sigma_P = lambda Sigma, Sigma_E = (1 - lambda) Sigma gives Sigma %x% (lambda R + (1 - lambda) I)", {
  set.seed(5)
  tree <- ape::rtree(15)
  qc <- .mip_build_Qc(tree)
  R_tip <- solve(as.matrix(qc$Qc))[1:qc$n, 1:qc$n]
  Sigma <- matrix(c(1.2, 0.4, 0.4, 0.7), 2)
  lam <- 0.6
  implied <- kronecker(lam * Sigma, R_tip) + kronecker((1 - lam) * Sigma,
                                                       diag(qc$n))
  R <- stats::cov2cor(ape::vcv(tree))[tree$tip.label, tree$tip.label]
  target <- kronecker(Sigma, lam * R + (1 - lam) * diag(qc$n))
  expect_lt(max(abs(implied - target)), 1e-8)
  # Per-trait lambda from the pair equals lam.
  SP <- lam * Sigma; SE <- (1 - lam) * Sigma
  expect_equal(diag(SP) / (diag(SP) + diag(SE)), c(lam, lam))
})

# ---- Item 6: only tip rows of a enter residuals and completions -------------
test_that("internal-node rows of a never enter the completed data", {
  d <- mip_sim(n = 15L, seed = 6L)
  Y <- d$Y
  Y[2:5, 1] <- NA; Y[6:9, 2] <- NA; Y[10, ] <- NA
  prob <- .mip_problem(Y, d$tree)
  tpl <- .mip_template(prob)
  lb <- .mip_lik_blocks(prob, d$SE)
  tpl <- .mip_refactor(tpl, d$SP, lb)
  s <- .mip_amu(prob, tpl, .mip_linear(prob, tpl, lb), 0L)
  expect_equal(dim(s$a), c(prob$N, 2L))
  a_tip <- s$a[seq_len(prob$n), , drop = FALSE]
  a_bad <- s$a
  a_bad[(prob$n + 1L):prob$N, ] <- 1e6       # corrupt internal rows only
  y1 <- .mip_fill_ymis(prob, a_tip, s$mu, d$SE, draw = FALSE)
  y2 <- .mip_fill_ymis(prob, a_bad[seq_len(prob$n), , drop = FALSE], s$mu,
                       d$SE, draw = FALSE)
  expect_identical(y1, y2)
  expect_identical(y1[!is.na(Y)], Y[!is.na(Y)])
  # The tip block of the mean is what the dense conditional predicts, so the
  # tips-first index mapping is right.
  expect_equal(.mip_cond_mean(prob, d$SP, d$SE, mu = d$mu),
               dense_cond(d$C, Y, d$mu)$mean, tolerance = 1e-8)
})

# ---- Exactness at fixed parameters (G2 in miniature) ------------------------
test_that("fixed-parameter draws match the dense conditional (mean and covariance)", {
  d <- mip_sim(n = 25L, seed = 7L)
  Y <- d$Y
  set.seed(8)
  Y[sample.int(25L, 6L), 1] <- NA; Y[sample.int(25L, 6L), 2] <- NA
  prob <- .mip_problem(Y, d$tree)
  ref <- dense_cond(d$C, Y, d$mu)
  set.seed(9)
  fx <- .mip_fixed_draws(prob, d$SP, d$SE, 6000L, mu = d$mu)
  mcse <- apply(fx$ymis, 1L, stats::sd) / sqrt(6000)
  expect_lt(max(abs(rowMeans(fx$ymis) - ref$mean) / mcse), 4.5)
  expect_lt(max(abs(stats::cov(t(fx$ymis)) - ref$cov)), 0.06)
  # Flat mu: the conditional mean is the GLS-marginal one.
  n <- 25L
  X1 <- kronecker(diag(2), matrix(1, n, 1))
  v <- as.vector(Y); mi <- which(is.na(v)); o <- which(!is.na(v))
  Ci <- solve(d$C[o, o]); Xo <- X1[o, ]
  bhat <- solve(t(Xo) %*% Ci %*% Xo, t(Xo) %*% Ci %*% v[o])
  m_flat <- X1[mi, ] %*% bhat + d$C[mi, o] %*% Ci %*% (v[o] - Xo %*% bhat)
  expect_equal(.mip_cond_mean(prob, d$SP, d$SE), as.numeric(m_flat),
               tolerance = 1e-8)
})

# ---- Collapsed likelihood ----------------------------------------------------
test_that("the collapsed log-likelihood matches the dense integrated likelihood up to a constant", {
  d <- mip_sim(n = 20L, seed = 10L)
  Y <- d$Y
  Y[1:4, 1] <- NA; Y[3:8, 2] <- NA
  prob <- .mip_problem(Y, d$tree)
  dense <- function(SP, SE) {
    C <- kronecker(SP, d$R) + kronecker(SE, diag(20))
    v <- as.vector(Y); o <- !is.na(v)
    C <- C[o, o]; X <- kronecker(diag(2), matrix(1, 20, 1))[o, ]; vo <- v[o]
    Ci <- solve(C); XtCX <- t(X) %*% Ci %*% X
    as.numeric(-0.5 * determinant(C)$modulus -
                 0.5 * determinant(XtCX)$modulus -
                 0.5 * (t(vo) %*% Ci %*% vo -
                          t(vo) %*% Ci %*% X %*% solve(XtCX, t(X) %*% Ci %*% vo)))
  }
  th <- list(list(d$SP, d$SE),
             list(matrix(c(1.5, -0.4, -0.4, 0.4), 2), diag(c(0.05, 0.8))),
             list(diag(c(0.2, 2)), matrix(c(1, 0.5, 0.5, 1), 2)))
  tpl <- .mip_template(prob)
  sp <- vapply(th, function(t) .mip_collapsed_ll(prob, tpl, t[[1]], t[[2]])$ll,
               numeric(1))
  dn <- vapply(th, function(t) dense(t[[1]], t[[2]]), numeric(1))
  expect_equal(diff(sp), diff(dn), tolerance = 1e-8)
})

# ---- Diagnostics (Vehtari et al. 2021) ---------------------------------------
test_that("bulk ESS and split R-hat behave on known chains", {
  set.seed(12)
  iid <- matrix(stats::rnorm(4000), 1000, 4)
  expect_gt(.mip_ess_bulk(iid), 3000)
  expect_lt(.mip_rhat(iid), 1.01)
  ar <- apply(matrix(stats::rnorm(20000), 5000, 4), 2, function(e) {
    stats::filter(e, 0.9, method = "recursive")
  })
  ess_ar <- .mip_ess_bulk(ar)
  expect_gt(ess_ar, 0.5 * 20000 * 0.1 / 1.9)   # theory 20000 (1 - 0.9)/(1 + 0.9)
  expect_lt(ess_ar, 2 * 20000 * 0.1 / 1.9)
  shifted <- iid + rep(c(0, 0, 0, 3), each = 1000)
  expect_gt(.mip_rhat(shifted), 1.1)
  d <- .mip_diagnostics(list(matrix(stats::rnorm(1000 * 10), 1000),
                             matrix(stats::rnorm(1000 * 10), 1000)), 2L)
  expect_identical(names(d), c("parameter", "rhat", "ess_bulk"))
  expect_true(isTRUE(attr(d, "converged")))
  expect_false(any(grepl("^mu", d$parameter)))
})

# ---- The full sampler: short chains, sanity only -----------------------------
test_that("a short full run returns finite, positive-definite parameter draws", {
  d <- mip_sim(n = 40L, seed = 13L)
  Y <- d$Y
  Y[1:8, 1] <- NA; Y[5:14, 2] <- NA
  ctl <- .mip_resolve_control(list(n_chains = 2L, burnin = 100L, n_iter = 200L,
                                   keep_draws = 40L, seed = 1L), m = 5L)
  f <- .mip_fit(Y, d$tree, ctl)
  expect_equal(dim(f$params$Sigma_P), c(2L, 2L, 40L))
  expect_equal(dim(f$ymis), c(18L, 40L))
  expect_true(all(is.finite(f$ymis)))
  mins <- apply(f$params$Sigma_P, 3L, function(S) min(eigen(S)$values))
  expect_true(all(mins > 0))
  expect_true(all(f$params$lambda > 0 & f$params$lambda < 1))
})

# ---- Item 11: improper mode holds the parameters fixed -----------------------
test_that("param_uncertainty = 'none' fixes Sigma_P and Sigma_E across draws", {
  d <- mip_sim(n = 30L, seed = 14L)
  Y <- d$Y
  Y[1:6, 1] <- NA; Y[4:10, 2] <- NA
  ctl <- .mip_resolve_control(list(n_chains = 2L, burnin = 50L, n_iter = 100L,
                                   keep_draws = 40L, seed = 2L,
                                   param_uncertainty = "none"), m = 5L)
  f <- .mip_fit(Y, d$tree, ctl)
  sp <- f$params$Sigma_P
  se <- f$params$Sigma_E
  expect_true(all(apply(sp, 3L, function(S) isTRUE(all.equal(S, sp[, , 1L])))))
  expect_true(all(apply(se, 3L, function(S) isTRUE(all.equal(S, se[, , 1L])))))
  expect_equal(stats::sd(f$params$lambda[, 1L]), 0)
  expect_gt(stats::sd(f$params$mu[, 1L]), 0)       # mu is still drawn
  expect_gt(min(apply(f$ymis, 1L, stats::sd)), 0)  # cells still vary
})
