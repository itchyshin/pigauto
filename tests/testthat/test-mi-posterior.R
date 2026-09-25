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

# ---- O(n) tip depths: .mip_build_Qc never forms the dense vcv ---------------
test_that("tip depths from node.depth.edgelength() match diag(vcv()) and Qc avoids vcv", {
  set.seed(15)
  for (tree in list(ape::rcoal(30), ape::rtree(30))) {
    n <- length(tree$tip.label)
    h_dense <- build_henderson_S_inv(tree)
    h_sparse <- build_henderson_S_inv(
      tree, tip_depths = ape::node.depth.edgelength(tree)[seq_len(n)])
    expect_lt(max(abs(h_sparse$tip_sqrt_d - h_dense$tip_sqrt_d)), 1e-12)
    expect_identical(names(h_sparse$tip_sqrt_d), tree$tip.label)
    expect_identical(h_sparse$Q, h_dense$Q)
  }
  expect_gt(stats::sd(h_dense$tip_sqrt_d), 0.01)   # the rtree is non-ultrametric
  # The posterior path gets the same values without calling ape::vcv().
  ref <- .mip_build_Qc(tree)
  local_mocked_bindings(vcv = function(...) stop("dense vcv formed"),
                        .package = "ape")
  expect_error(ape::vcv(tree), "dense vcv formed")
  qc <- .mip_build_Qc(tree)
  expect_lt(max(abs(qc$tip_sqrt_d - h_dense$tip_sqrt_d)), 1e-12)
  expect_identical(qc$Qc, ref$Qc)
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
test_that("the completion keeps observed cells and the tips-first mapping is right", {
  d <- mip_sim(n = 15L, seed = 6L)
  Y <- d$Y
  Y[2:5, 1] <- NA; Y[6:9, 2] <- NA; Y[10, ] <- NA
  prob <- .mip_problem(Y, d$tree)
  tpl <- .mip_template(prob)
  lb <- .mip_lik_blocks(prob, d$SE)
  tpl <- .mip_refactor(tpl, d$SP, lb)
  s <- .mip_amu(prob, tpl, .mip_linear(prob, tpl, lb), 0L)
  expect_equal(dim(s$a), c(prob$N, 2L))
  y1 <- .mip_fill_ymis(prob, s$a[seq_len(prob$n), , drop = FALSE], s$mu,
                       d$SE, draw = FALSE)
  expect_identical(y1[!is.na(Y)], Y[!is.na(Y)])
  # The tip block of the mean is what the dense conditional predicts, so the
  # tips-first index mapping is right.
  expect_equal(.mip_cond_mean(prob, d$SP, d$SE, mu = d$mu),
               dense_cond(d$C, Y, d$mu)$mean, tolerance = 1e-8)
})

test_that("inside the sampler, only tip rows of a reach the completion, alpha and Sigma_E steps", {
  # One Gibbs sweep of .mip_run_chain (no burn-in, no Metropolis moves), with
  # probes on the functions that receive rows of a. `internal` overwrites the
  # internal-node rows of the step-1 draw of a before the sweep uses it.
  d <- mip_sim(n = 15L, seed = 6L)
  Y <- d$Y
  Y[2:5, 1] <- NA; Y[6:9, 2] <- NA; Y[10, ] <- NA
  prob <- .mip_problem(Y, d$tree)
  n <- prob$n; N <- prob$N
  hyper <- list(nu_W = 3, S_W = diag(2), V_alpha = 1000, nu_E = 3,
                S_E = diag(0.01 * prob$obs_var, 2))
  start <- list(Sigma_P = d$SP, Sigma_E = d$SE)
  real_amu <- .mip_amu; real_fill <- .mip_fill_ymis
  real_alpha <- .mip_draw_alpha; real_iw <- .mip_riwish
  one_sweep <- function(internal = NULL) {
    cap <- new.env()
    cap$iw <- list()
    local_mocked_bindings(
      .mip_amu = function(prob, tpl, b, n_draws = 1L) {
        s <- real_amu(prob, tpl, b, n_draws)
        if (!is.null(internal)) s$a[(n + 1L):N, ] <- internal
        cap$a <- s$a
        s
      },
      .mip_fill_ymis = function(prob, a_tip, mu, Sigma_E, draw = TRUE) {
        cap$fill_a <- a_tip
        real_fill(prob, a_tip, mu, Sigma_E, draw)
      },
      .mip_draw_alpha = function(xi_tip, Rres, Sigma_E, V_alpha) {
        al <- real_alpha(xi_tip, Rres, Sigma_E, V_alpha)
        cap$xi_tip <- xi_tip; cap$Rres <- Rres; cap$alpha <- al
        al
      },
      .mip_riwish = function(nu, S) {
        cap$iw[[length(cap$iw) + 1L]] <- S
        real_iw(nu, S)
      }
    )
    ch <- .mip_run_chain(prob, start, hyper, burnin = 0L, n_iter = 1L,
                         thin = 1L, seed = 1L, mh = FALSE, gibbs = TRUE)
    list(ch = ch, cap = cap)
  }
  clean <- one_sweep()
  dirty <- one_sweep(matrix(seq(-50, 50, length.out = 2L * (N - n)), N - n, 2L))
  # The corrupted rows did reach the sampler: step 4 (Sigma_W | xi) uses
  # every node, so its scale matrix and Sigma_P change.
  expect_false(isTRUE(all.equal(clean$cap$iw[[1L]], dirty$cap$iw[[1L]])))
  expect_false(isTRUE(all.equal(clean$ch$Sigma_P, dirty$ch$Sigma_P)))
  # The completion (step 1b), alpha (step 3) and Sigma_E (step 5) do not.
  expect_identical(dirty$ch$ymis, clean$ch$ymis)
  expect_identical(dirty$ch$mu, clean$ch$mu)
  expect_identical(dirty$cap$alpha, clean$cap$alpha)
  expect_identical(dirty$ch$Sigma_E, clean$ch$Sigma_E)
  # Exact row checks, which also catch a permutation of the tip rows.
  for (r in list(clean, dirty)) {
    a_tip <- r$cap$a[seq_len(n), , drop = FALSE]
    expect_identical(r$cap$fill_a, a_tip)          # step 1b
    expect_identical(r$cap$xi_tip, a_tip)          # step 3 (alpha = 1 at start)
    Yc <- prob$Y
    Yc[prob$miss] <- r$ch$ymis[, 1L]
    expect_equal(r$cap$Rres, sweep(Yc, 2L, r$ch$mu[1L, ]), ignore_attr = TRUE)
    e <- r$cap$Rres - sweep(a_tip, 2L, r$cap$alpha, "*")   # step 5 residual
    expect_equal(r$cap$iw[[2L]], hyper$S_E + crossprod(e), tolerance = 1e-12,
                 ignore_attr = TRUE)
  }
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
  # Same location, one chain three times wider: only the folded (tail)
  # R-hat detects it (Vehtari et al. 2021), so this fails if .mip_rhat()
  # computed only the bulk version.
  wide <- iid
  wide[, 4] <- wide[, 4] * 3
  expect_lt(.mip_rhat_basic(.mip_rank_normalise(.mip_split_chains(wide))),
            1.01)
  expect_gt(.mip_rhat(wide), 1.05)
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

# ---- Kernel agreement: Metropolis-only vs Gibbs-only (design.md 2.5) ---------
test_that("the Metropolis-only and Gibbs-only kernels target the same posterior", {
  skip_on_cran()
  # Two valid kernels for the same target, run through .mip_run_chain's test
  # hooks: (a) the collapsed Metropolis moves plus step 1 (mh = TRUE,
  # gibbs = FALSE); (b) the parameter-expanded Gibbs steps 1 to 5 (mh =
  # FALSE, gibbs = TRUE). A wrong Jacobian, prior term, degrees of freedom or
  # node precision in either kernel moves its stationary distribution away
  # from the other's. Fixture chosen so that both kernels mix (posterior
  # lambda near 0.5). About 25 s.
  d <- mip_sim(n = 80L, seed = 13L, ultrametric = TRUE,
               SP = matrix(c(0.5, 0.2, 0.2, 0.5), 2),
               SE = matrix(c(0.5, 0.15, 0.15, 0.5), 2))
  Y <- d$Y
  Y[1:8, 1] <- NA; Y[5:14, 2] <- NA
  prob <- .mip_problem(Y, d$tree)
  hyper <- list(nu_W = 3, S_W = diag(2), V_alpha = 1000, nu_E = 3,
                S_E = diag(0.01 * prob$obs_var, 2))
  set.seed(4)
  starts <- .mip_starts(prob, rep(0.5, 2), 4L)
  derive <- function(tr) {
    colnames(tr) <- .mip_param_names(2L)
    cbind(lambda1 = tr[, "lambda[1]"], lambda2 = tr[, "lambda[2]"],
          logSP11 = log(tr[, "Sigma_P[1,1]"]),
          logSP22 = log(tr[, "Sigma_P[2,2]"]),
          logSE11 = log(tr[, "Sigma_E[1,1]"]),
          logSE22 = log(tr[, "Sigma_E[2,2]"]),
          rhoP = tr[, "Sigma_P[1,2]"] /
            sqrt(tr[, "Sigma_P[1,1]"] * tr[, "Sigma_P[2,2]"]),
          rhoE = tr[, "Sigma_E[1,2]"] /
            sqrt(tr[, "Sigma_E[1,1]"] * tr[, "Sigma_E[2,2]"]))
  }
  kernel <- function(mh, gibbs, n_iter, seed0) {
    lapply(seq_len(4L), function(c) {
      derive(.mip_run_chain(prob, starts[[c]], hyper, burnin = 200L,
                            n_iter = n_iter, thin = n_iter, seed = seed0 + c,
                            mh = mh, gibbs = gibbs)$trace)
    })
  }
  # Posterior mean and its MCSE = sd / sqrt(bulk ESS) per quantity.
  summarise <- function(chains) {
    t(vapply(seq_len(ncol(chains[[1L]])), function(j) {
      s <- vapply(chains, function(x) x[, j], numeric(nrow(chains[[1L]])))
      c(mean = mean(s), mcse = stats::sd(as.vector(s)) / sqrt(.mip_ess_bulk(s)))
    }, numeric(2L)))
  }
  mh_only <- kernel(TRUE, FALSE, 1000L, 100L)
  gibbs_only <- kernel(FALSE, TRUE, 4000L, 200L)
  a <- summarise(mh_only)
  b <- summarise(gibbs_only)
  z <- (a[, "mean"] - b[, "mean"]) / sqrt(a[, "mcse"]^2 + b[, "mcse"]^2)
  names(z) <- colnames(mh_only[[1L]])
  expect_true(all(is.finite(z)))
  expect(max(abs(z)) < 3.5,
         paste0("Metropolis-only and Gibbs-only kernels disagree: |z| = ",
                paste(sprintf("%s %.2f", names(z), abs(z)), collapse = ", ")))
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

# ---- Automatic chain extension (design.md 5e) --------------------------------
test_that("posterior_control validates auto_extend and max_extend", {
  def <- .mip_resolve_control(list(), m = 5L)
  expect_identical(def$auto_extend, TRUE)
  expect_identical(def$max_extend, 3L)
  expect_identical(.mip_resolve_control(list(max_extend = 10), m = 5L)$max_extend,
                   10L)
  off <- .mip_resolve_control(list(auto_extend = FALSE, max_extend = 0L), m = 5L)
  expect_identical(off$auto_extend, FALSE)
  expect_identical(off$max_extend, 0L)
  for (bad in list(NA, "yes", c(TRUE, FALSE), 1, logical(0))) {
    expect_error(.mip_resolve_control(list(auto_extend = bad), m = 5L),
                 "`posterior_control\\$auto_extend` must be TRUE or FALSE")
  }
  for (bad in list(-1, 11, 1.5, NA, NA_integer_, "2", c(1, 2), Inf,
                   integer(0))) {
    expect_error(.mip_resolve_control(list(max_extend = bad), m = 5L),
                 "`posterior_control\\$max_extend` must be a whole number from 0 to 10")
  }
})

test_that("a resumed chain is the same chain as one run longer from the start", {
  # Before any subsampling: segment 1 (n_iter = 31) plus a continuation from
  # its saved sampler and RNG state (n_iter = 31) against one run with
  # n_iter = 62 and the same seed and thin. 31 is not a multiple of thin = 3,
  # so the kept-sweep count has to carry across the segment boundary.
  # burnin = 100 so the step sizes adapt (at sweeps 50 and 100) and the
  # resumed segment has to restore the adapted values, not the start values.
  d <- mip_sim(n = 25L, seed = 3L)
  Y <- d$Y
  Y[1:5, 1] <- NA; Y[4:9, 2] <- NA
  prob <- .mip_problem(Y, d$tree)
  hyper <- list(nu_W = 3, S_W = diag(2), V_alpha = 1000, nu_E = 3,
                S_E = diag(0.01 * prob$obs_var, 2))
  start <- list(Sigma_P = d$SP, Sigma_E = d$SE)
  long <- .mip_run_chain(prob, start, hyper, burnin = 100L, n_iter = 62L,
                         thin = 3L, seed = 7L)
  seg1 <- .mip_run_chain(prob, start, hyper, burnin = 100L, n_iter = 31L,
                         thin = 3L, seed = 7L)
  # Adaptation really happened, so restoring the step sizes is exercised.
  expect_false(all(seg1$state$step == 0.5))
  expect_false(all(seg1$state$step_off == 0.2))
  set.seed(99)
  stats::runif(5L)          # the caller's RNG stream moves on in between
  seg2 <- .mip_run_chain(prob, NULL, hyper, burnin = 100L, n_iter = 31L,
                         thin = 3L, resume = seg1$state)
  expect_identical(dim(seg1$ymis), c(nrow(prob$miss), 10L))   # sweeps 3..30
  expect_identical(dim(seg2$ymis), c(nrow(prob$miss), 10L))   # sweeps 33..60
  joined <- .mip_join_segments(seg1, seg2)
  for (k in c("trace", "ymis", "Sigma_P", "Sigma_E", "mu", "mh_step",
              "mh_accept")) {
    expect_identical(joined[[k]], long[[k]], info = k)
  }
  expect_identical(joined$state$j, 62L)
  # The global RNG stream ends where the single long run left it.
  expect_identical(get(".Random.seed", envir = globalenv()), long$state$rng)
})

test_that("an extended fit equals one run longer with the thinning scaled", {
  # After subsampling: n_iter = 30, thin = 3 and one extension against
  # n_iter = 60, thin = 6 and no extension. 2 chains x 60 sweeps cannot
  # reach bulk ESS 400, so exactly one extension is made. "both" also
  # checks the plug-in draws, whose covariances are the posterior means
  # over the whole extended run.
  d <- mip_sim(n = 25L, seed = 3L)
  Y <- d$Y
  Y[1:5, 1] <- NA; Y[4:9, 2] <- NA
  base <- list(n_chains = 2L, burnin = 100L, keep_draws = 20L, seed = 4L,
               param_uncertainty = "both")
  ext <- .mip_fit(Y, d$tree, .mip_resolve_control(
    c(base, list(n_iter = 30L, thin = 3L, max_extend = 1L)), m = 4L))
  one <- .mip_fit(Y, d$tree, .mip_resolve_control(
    c(base, list(n_iter = 60L, thin = 6L, auto_extend = FALSE)), m = 4L))
  expect_identical(attr(ext$diagnostics, "n_extensions"), 1L)
  expect_identical(attr(one$diagnostics, "n_extensions"), 0L)
  expect_identical(attr(ext$diagnostics, "sweeps_per_chain"), 60L)
  expect_identical(attr(one$diagnostics, "sweeps_per_chain"), 60L)
  expect_false(isTRUE(attr(ext$diagnostics, "converged")))
  expect_identical(dim(ext$ymis), c(nrow(ext$miss), 20L))  # 10 per chain
  expect_identical(ext$ymis, one$ymis)
  expect_identical(ext$params, one$params)
  expect_identical(ext$improper, one$improper)
  expect_identical(ext$sweeps, one$sweeps)
  expect_identical(ext$sweeps, 2L * (100L + 60L))
  dg <- ext$diagnostics
  attr(dg, "n_extensions") <- 0L
  expect_identical(dg, one$diagnostics)
  expect_equal(ext$improper$Sigma_P,
               apply(ext$params$Sigma_P, c(1L, 2L), mean), ignore_attr = TRUE)
})

test_that("the first run is pinned to the campaign code (69670d4) on a fixed fixture", {
  # Regression guard for design.md 5e: the re-run campaign reuses 4,756 cells
  # computed at commit 69670d4, so the first-run path of the sampler must not
  # drift. Reference values were generated by this exact fixture at 69670d4
  # and at the commit that added this test (identical); tolerance 1e-8 allows
  # last-bit BLAS differences across platforms.
  skip_on_cran()
  set.seed(11)
  tree <- ape::rtree(40)
  R <- stats::cov2cor(ape::vcv(tree))
  SP <- matrix(c(0.6, 0.3, 0.3, 0.5), 2)
  Y <- t(chol(R + diag(1e-10, 40))) %*% matrix(stats::rnorm(80), 40) %*% chol(SP) +
    matrix(stats::rnorm(80, sd = 0.5), 40)
  dimnames(Y) <- list(tree$tip.label, c("t1", "t2"))
  Y[c(3, 7, 12, 20, 31), 1] <- NA
  Y[c(5, 12, 18, 27), 2] <- NA
  ctl <- .mip_resolve_control(list(n_chains = 3L, burnin = 120L, n_iter = 90L,
                                   keep_draws = 30L, seed = 9L), m = 5L)
  ctl$auto_extend <- FALSE
  f <- .mip_fit(Y, tree, ctl)
  ref_ymis <- matrix(c(-0.421427778495496, -0.0842243460337871, -1.00837690913468,
                       -0.352080602256646, 0.278600117017885, -0.17153508451454,
                       -0.700673544028394, 0.359473689403797, -0.991244807893159), 3L)
  ref_lambda <- matrix(c(0.940063997170256, 0.984335306033873, 0.999429438633398,
                         0.997665520927872, 0.994380395228885, 0.995454170036249), 3L)
  expect_equal(unname(f$ymis[1:3, 1:3]), ref_ymis, tolerance = 1e-8)
  expect_equal(unname(f$params$lambda[1:3, ]), ref_lambda, tolerance = 1e-8)
  expect_false(isTRUE(attr(f$diagnostics, "converged")))
})
