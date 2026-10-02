#!/usr/bin/env Rscript
# script/mi_gls/gate_exactness.R
#
# G2 (docs/dev-log/mi-posterior/design.md section 5b): with the parameters
# fixed, step 1 of the posterior sampler is an exact draw of y_mis | y_obs.
#
# G2a  n = 40, K = 2, Sigma_P and Sigma_E positive definite and fixed,
#      Sigma_E not small, mu fixed (internal hook: .mip_fixed_draws(mu = )).
#      20,000 draws of y_mis compared with the dense conditional of
#      y_mis | y_obs under N(1 mu', Sigma_P %x% R + Sigma_E %x% I_n):
#        - every cell mean within 4 Monte Carlo SE;
#        - every entry of the y_mis covariance within 0.03.
# G2b  Sigma_E = 1e-6 I, mu = 0, Sigma_P = the prototype's Sigma estimate:
#      the exact conditional mean of y_mis (latent scale) matches
#      draw_conditional_bm()$mu_cond within 1e-3 at every missing cell.
#
# Prints EXACTNESS_OK only if both pass; otherwise the failing numbers.
# Usage (cwd = worktree root): Rscript script/mi_gls/gate_exactness.R

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
suppressMessages(devtools::load_all(quiet = TRUE))

ok_all <- TRUE

# ---- G2a --------------------------------------------------------------------
set.seed(20260924)
n <- 40L
tree <- ape::rtree(n)
R <- stats::cov2cor(ape::vcv(tree))[tree$tip.label, tree$tip.label]
SP <- matrix(c(0.6, 0.3, 0.3, 0.5), 2)
SE <- matrix(c(0.3, 0.1, 0.1, 0.4), 2)
mu <- c(0.2, -0.1)
C <- kronecker(SP, R) + kronecker(SE, diag(n))
y <- as.vector(t(chol(C)) %*% stats::rnorm(2L * n)) + rep(mu, each = n)
Y <- matrix(y, n, 2L, dimnames = list(tree$tip.label, c("t1", "t2")))
Y[sample.int(n, 10L), 1L] <- NA
Y[sample.int(n, 10L), 2L] <- NA
Y[1:3, ] <- NA                      # tips with nothing observed
prob <- .mip_problem(Y, tree)
v <- as.vector(Y)
mi <- which(is.na(v)); o <- which(!is.na(v))
stopifnot(identical(mi, as.integer((prob$miss[, 2L] - 1L) * n + prob$miss[, 1L])))
m_dense <- as.numeric(rep(mu, each = n)[mi] +
  C[mi, o] %*% solve(C[o, o], v[o] - rep(mu, each = n)[o]))
V_dense <- C[mi, mi] - C[mi, o] %*% solve(C[o, o], C[o, mi])

n_draws <- 20000L
t0 <- proc.time()[["elapsed"]]
fx <- .mip_fixed_draws(prob, SP, SE, n_draws, mu = mu)
wall <- proc.time()[["elapsed"]] - t0
mcse <- apply(fx$ymis, 1L, stats::sd) / sqrt(n_draws)
z <- (rowMeans(fx$ymis) - m_dense) / mcse
cov_err <- abs(stats::cov(t(fx$ymis)) - V_dense)
g2a_mean <- max(abs(z)) < 4
g2a_cov <- max(cov_err) < 0.03
cat(sprintf(paste0("G2a: %d missing cells, %d draws (%.1fs). max |mean - ",
                   "dense| / MCSE = %.3f (< 4: %s); max |cov - dense| = %.4f ",
                   "(< 0.03: %s)\n"),
            length(mi), n_draws, wall, max(abs(z)), g2a_mean,
            max(cov_err), g2a_cov))
if (!g2a_mean) {
  bad <- which(abs(z) >= 4)
  cat("G2a FAILED mean cells:\n")
  print(data.frame(cell = mi[bad], dense = m_dense[bad],
                   sampled = rowMeans(fx$ymis)[bad], z = z[bad]))
}
if (!g2a_cov) {
  w <- which(cov_err >= 0.03, arr.ind = TRUE)
  cat("G2a FAILED covariance entries (first 10):\n")
  print(utils::head(data.frame(i = w[, 1L], j = w[, 2L],
                               dense = V_dense[w],
                               sampled = stats::cov(t(fx$ymis))[w]), 10L))
}
ok_all <- ok_all && g2a_mean && g2a_cov

# ---- G2b --------------------------------------------------------------------
set.seed(7)
n2 <- 60L
tree2 <- ape::rtree(n2)
R2 <- stats::cov2cor(ape::vcv(tree2))
S0 <- matrix(c(1, 0.6, 0.6, 1), 2)
Y2 <- t(chol(R2 + diag(1e-10, n2))) %*% matrix(stats::rnorm(2L * n2), n2) %*%
  chol(S0)
df <- data.frame(t1 = Y2[, 1L], t2 = Y2[, 2L], row.names = tree2$tip.label)
df$t1[sample.int(n2, 15L)] <- NA
df$t2[sample.int(n2, 15L)] <- NA
pd <- preprocess_traits(df, tree2, log_transform = FALSE)
dcb <- draw_conditional_bm(list(data = pd, tree = tree2), m = 2L, seed = 1L)
X <- pd$X_scaled
rownames(X) <- pd$species_names
X <- X[tree2$tip.label, , drop = FALSE]
prob2 <- .mip_problem(X, tree2)
cm <- .mip_cond_mean(prob2, dcb$sigma_hat, diag(1e-6, 2L), mu = c(0, 0))
proto <- dcb$mu_cond[tree2$tip.label, , drop = FALSE][prob2$miss]
g2b_err <- max(abs(cm - proto))
g2b <- g2b_err < 1e-3
cat(sprintf(paste0("G2b: %d missing cells, Sigma_E = 1e-6 I: max |posterior ",
                   "conditional mean - draw_conditional_bm()| = %.2e ",
                   "(< 1e-3: %s)\n"),
            nrow(prob2$miss), g2b_err, g2b))
if (!g2b) {
  bad <- which(abs(cm - proto) >= 1e-3)
  cat("G2b FAILED cells:\n")
  print(data.frame(row = prob2$miss[bad, 1L], col = prob2$miss[bad, 2L],
                   posterior = cm[bad], prototype = proto[bad]))
}
ok_all <- ok_all && g2b

if (ok_all) {
  cat("EXACTNESS_OK\n")
} else {
  cat("G2 FAILED\n")
  quit(save = "no", status = 1L)
}
