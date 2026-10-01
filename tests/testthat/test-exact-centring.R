# Mean-model consistency for predict_method = "exact" (S1,
# docs/dev-log/exact-default/S1-centring-report.md).
#
# exact_conditional_mvn() assumes vec(L) ~ MVN(0, Sigma %x% R(lambda_block)) --
# a ZERO column mean -- the same zero-root assumption henderson_bm_predict()
# makes for K = 1 (R/joint_mvn_solver.R, .mvn_init_per_column()'s own
# centring comment). `exact_centre` (default TRUE in fit_mvn_bm_inhouse())
# centres each column at its own GLS phylogenetic mean at lambda_block before
# the exact solve and adds the mean back afterwards; `exact_centre = FALSE`
# reproduces the pre-S1 (implicit zero-mean) behaviour exactly.

# Shared DGP helper for (b) and (c): draws L ~ MVN(0, Sigma %x% R(lambda)) + a
# constant column shift, from a random coalescent tree, then MCARs a
# fraction `miss` of cells (with a floor of >= 5 observed cells per column so
# the GLS mean / exact solve stay well-posed).
.centring_dgp <- function(seed, n = 200L, K = 3L, lambda_true = 0.3,
                           shift = 0, miss = 0.3) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv(tree))
  R_lambda <- lambda_true * R
  diag(R_lambda) <- lambda_true * diag(R) + (1 - lambda_true)
  Sigma <- matrix(c(1, .5, -.2, .5, 1, .3, -.2, .3, 1), K, K)
  Z <- matrix(stats::rnorm(n * K), n, K)
  L_true <- t(chol(R_lambda)) %*% Z %*% chol(Sigma) + shift
  rownames(L_true) <- tree$tip.label
  colnames(L_true) <- c("a", "b", "c")

  set.seed(seed + 10000L)
  miss_mask <- matrix(stats::runif(n * K) < miss, n, K)
  for (j in seq_len(K)) {
    if (sum(!miss_mask[, j]) < 5L) {
      obs_idx <- sample(seq_len(n), 5L)
      miss_mask[obs_idx, j] <- FALSE
    }
  }
  L_obs <- L_true
  L_obs[miss_mask] <- NA
  list(tree = tree, L_true = L_true, L_obs = L_obs, miss_mask = miss_mask,
       lambda_true = lambda_true)
}

test_that("[exact-centre] exact_centre = FALSE reproduces the pre-S1 exact output to 1e-12", {
  skip_if_not_installed("Matrix")
  # Fixture and reference values were captured by running the UNMODIFIED
  # (pre-S1) fit_mvn_bm_inhouse(predict_method = "exact") on this exact
  # seed/tree/L, i.e. before exact_centre existed (equivalent to today's
  # exact_centre = FALSE, since that flag reproduces the old code path
  # verbatim -- no centring, no add-back).
  set.seed(42)
  n <- 20L; K <- 3L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv(tree))
  Sigma_true <- matrix(c(1, .5, -.2, .5, 1, .3, -.2, .3, 1), K, K)
  L <- t(chol(R)) %*% matrix(stats::rnorm(n * K), n, K) %*% chol(Sigma_true)
  rownames(L) <- tree$tip.label
  colnames(L) <- c("a", "b", "c")
  set.seed(43)
  L[matrix(stats::runif(n * K) < 0.3, n, K)] <- NA

  fit <- fit_mvn_bm_inhouse(L, tree = tree, predict_method = "exact",
                             exact_centre = FALSE)

  expect_equal(fit$anc_recon["t1", "a"],  -0.781828470664128, tolerance = 1e-12)
  expect_equal(fit$anc_recon["t13", "b"], -2.16498808350981,  tolerance = 1e-12)
  expect_equal(fit$anc_recon["t19", "c"],  0.571073617764007, tolerance = 1e-12)
  expect_equal(fit$anc_recon["t20", "a"], -0.124932501524891, tolerance = 1e-12)
  expect_equal(fit$anc_recon["t11", "b"], -1.06300238214169,  tolerance = 1e-12)
  expect_equal(fit$anc_recon["t9",  "c"],  0.184162871117541, tolerance = 1e-12)
})

test_that("[exact-centre] centring reduces missing-cell RMSE under a mean-shifted DGP", {
  skip_if_not_installed("Matrix")
  # lambda = 0.3 is the exact value at which Rose's review measured the
  # per-column path's own pre-centring bias (up to 0.25 in mu on
  # non-centered data; docs/dev-log/lambda-default/S2-solver-report.md).
  # miss = 0.85 (rather than a lower fraction) is deliberate: at low
  # missingness the uncentred exact conditional partially self-corrects
  # the mean by averaging over many weakly-correlated observed cells (a
  # law-of-large-numbers effect measured while designing this test --
  # at miss = 0.3 the centred/uncentred RMSE gap was within noise, 10/20
  # seed-wins each way). High missingness leaves too few observed cells
  # for that self-correction, which is exactly the regime this fix
  # targets.
  n_seeds <- 20L
  rmse_c  <- numeric(n_seeds)
  rmse_uc <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    d <- .centring_dgp(seed = s, shift = 1.5, miss = 0.85)
    fit_c  <- fit_mvn_bm_inhouse(d$L_obs, tree = d$tree, predict_method = "exact",
                                  lambda = d$lambda_true, exact_centre = TRUE)
    fit_uc <- fit_mvn_bm_inhouse(d$L_obs, tree = d$tree, predict_method = "exact",
                                  lambda = d$lambda_true, exact_centre = FALSE)
    rmse_c[s]  <- sqrt(mean((fit_c$anc_recon[d$miss_mask]  - d$L_true[d$miss_mask])^2))
    rmse_uc[s] <- sqrt(mean((fit_uc$anc_recon[d$miss_mask] - d$L_true[d$miss_mask])^2))
  }
  # Measured (2026-09-25): mean centred RMSE 0.893, mean uncentred RMSE
  # 0.929, centred wins 16/20 seeds.
  expect_lt(mean(rmse_c), mean(rmse_uc))
})

test_that("[exact-centre] centred and uncentred nearly agree under a zero-mean DGP", {
  skip_if_not_installed("Matrix")
  # No shift: the model's own zero-mean assumption is (approximately) true,
  # so exact_centre should change little -- the GLS mean estimate itself is
  # close to zero, so subtracting and adding it back is close to a no-op.
  # "Close" is bounded by SAMPLING noise in a single realisation's GLS mean
  # (finite n_obs even under a true-zero DGP), not by an exact-zero
  # guarantee -- measured max |centred - uncentred| across 20 seeds was
  # 0.283 (mean 0.087); the assertion below uses a tolerance with margin
  # above that measurement.
  n_seeds <- 20L
  max_diff <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    d <- .centring_dgp(seed = s, shift = 0, miss = 0.3)
    fit_c  <- fit_mvn_bm_inhouse(d$L_obs, tree = d$tree, predict_method = "exact",
                                  lambda = d$lambda_true, exact_centre = TRUE)
    fit_uc <- fit_mvn_bm_inhouse(d$L_obs, tree = d$tree, predict_method = "exact",
                                  lambda = d$lambda_true, exact_centre = FALSE)
    max_diff[s] <- max(abs(fit_c$anc_recon - fit_uc$anc_recon))
  }
  expect_lt(max(max_diff), 0.35)
})

test_that("[exact-centre] fallback to per-column still works when the exact solve is unusable", {
  skip_if_not_installed("Matrix")
  # fit_mvn_bm_inhouse() doesn't expose exact_conditional_mvn()'s own
  # max_cells argument (that function's "refuses oversized problems" case
  # is already covered directly in test-exact-conditional.R); what this
  # test protects is the SURROUNDING fallback logic in fit_mvn_bm_inhouse
  # -- ec <- exact_conditional_mvn(...); if NULL, warn and fall through to
  # the per-column path -- still working correctly now that a centring step
  # runs before that call. Mocking exact_conditional_mvn() to return NULL
  # reproduces exactly the oversized-problem return value without needing
  # a genuinely huge problem.
  local_mocked_bindings(exact_conditional_mvn = function(...) NULL)

  set.seed(7)
  n <- 25L; K <- 3L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv(tree))
  Sigma <- matrix(c(1, .4, -.1, .4, 1, .2, -.1, .2, 1), K, K)
  L <- t(chol(R)) %*% matrix(stats::rnorm(n * K), n, K) %*% chol(Sigma)
  rownames(L) <- tree$tip.label
  L[matrix(stats::runif(n * K) < 0.3, n, K)] <- NA

  expect_warning(
    fit_ex <- fit_mvn_bm_inhouse(L, tree = tree, predict_method = "exact"),
    "falling back to the per-column path"
  )
  fit_pc <- fit_mvn_bm_inhouse(L, tree = tree, predict_method = "per_column")

  expect_equal(fit_ex$anc_recon, fit_pc$anc_recon)
  expect_equal(fit_ex$anc_var,   fit_pc$anc_var)
  expect_true(all(is.finite(fit_ex$anc_recon)))
})
