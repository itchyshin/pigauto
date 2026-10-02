# Calibration diagnostic for the posterior sampler (not a gate).
# Same DGP as script/mi_gls/gate_recovery.R; records 95% credible-interval
# coverage of the true lambda_k and phylogenetic correlation, plus REML lambda.
# Usage (cwd = code dir): Rscript g3_calibration.R <settings "B,C,D"> <seed_from> <seed_to> <cores> <out.rds>
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressMessages(devtools::load_all(quiet = TRUE))
a <- commandArgs(trailingOnly = TRUE)
sel <- strsplit(a[1], ",")[[1]]; seeds <- as.integer(a[2]):as.integer(a[3])
cores <- as.integer(a[4]); out <- a[5]
n <- 1000L; miss_frac <- 0.25
settings <- data.frame(id = c("A", "B", "C", "D"), lambda1 = c(1, 0.5, 0.3, 0.05),
                       lambda2 = c(1, 0.5, 0.9, 0.95), rho = c(0.7, 0.7, 0, 0.7))
settings <- settings[settings$id %in% sel, ]
sim <- function(lambda, rho, seed) {
  set.seed(seed)
  tree <- ape::rtree(n)
  R <- stats::cov2cor(ape::vcv(tree))
  Cr <- matrix(c(1, rho, rho, 1), 2)
  SP <- diag(sqrt(lambda)) %*% Cr %*% diag(sqrt(lambda))
  Lr <- t(chol(R + diag(1e-10, n)))
  Y <- Lr %*% matrix(stats::rnorm(2 * n), n) %*% chol(SP)
  Y <- Y + sweep(matrix(stats::rnorm(2 * n), n), 2L, sqrt(1 - lambda), "*")
  dimnames(Y) <- list(tree$tip.label, c("t1", "t2"))
  for (k in 1:2) Y[sample.int(n, round(miss_frac * n)), k] <- NA
  list(Y = Y, tree = tree)
}
ctl <- .mip_resolve_control(list(), m = 20L)
jobs <- expand.grid(s = seq_len(nrow(settings)), seed = seeds)
one <- function(j) {
  st <- settings[jobs$s[j], ]; seed <- jobs$seed[j]
  lam <- c(st$lambda1, st$lambda2)
  d <- sim(lam, st$rho, seed)
  c1 <- ctl; c1$seed <- seed
  f <- .mip_fit(d$Y, d$tree, c1)
  L <- f$params$lambda
  rho <- apply(f$params$Sigma_P, 3L, function(S) stats::cov2cor(S)[1L, 2L])
  q <- function(x) stats::quantile(x, c(0.025, 0.975), names = FALSE)
  data.frame(setting = st$id, seed = seed,
    lam1_true = lam[1], lam1_mean = mean(L[, 1]), lam1_sd = sd(L[, 1]),
    lam1_lo = q(L[, 1])[1], lam1_hi = q(L[, 1])[2], lam1_reml = f$start$lambda_reml[1],
    lam2_true = lam[2], lam2_mean = mean(L[, 2]), lam2_sd = sd(L[, 2]),
    lam2_lo = q(L[, 2])[1], lam2_hi = q(L[, 2])[2], lam2_reml = f$start$lambda_reml[2],
    rho_true = st$rho, rho_mean = mean(rho), rho_sd = sd(rho), rho_lo = q(rho)[1], rho_hi = q(rho)[2],
    converged = isTRUE(attr(f$diagnostics, "converged")),
    min_ess = min(f$diagnostics$ess_bulk))
}
rows <- parallel::mclapply(seq_len(nrow(jobs)), function(j) tryCatch(one(j), error = function(e) NULL),
                           mc.cores = cores)
res <- do.call(rbind, rows)
saveRDS(res, out)
cat(sprintf("fits: %d of %d\n", nrow(res), nrow(jobs)))
for (s in unique(res$setting)) {
  r <- res[res$setting == s, ]
  cv <- function(t, lo, hi) mean(t >= lo & t <= hi)
  cat(sprintf("%s (n=%d): lam1 cover %.3f  sd(mean) %.3f  mean(post sd) %.3f | lam2 cover %.3f  sd(mean) %.3f  mean(post sd) %.3f | rho cover %.3f  sd(mean) %.3f  mean(post sd) %.3f | REML sd %.3f %.3f | converged %d\n",
    s, nrow(r), cv(r$lam1_true, r$lam1_lo, r$lam1_hi), sd(r$lam1_mean), mean(r$lam1_sd),
    cv(r$lam2_true, r$lam2_lo, r$lam2_hi), sd(r$lam2_mean), mean(r$lam2_sd),
    cv(r$rho_true, r$rho_lo, r$rho_hi), sd(r$rho_mean), mean(r$rho_sd),
    sd(r$lam1_reml), sd(r$lam2_reml), sum(r$converged)))
}
