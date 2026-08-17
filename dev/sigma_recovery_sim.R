# Known-Sigma recovery simulation (Tier 2.2) — gates the Fisher-ML port.
# Design (pre-registered): docs/dev-log/2026-08-17-sigma-recovery-design.md
# Estimators: E0 = fit_mvn_bm_inhouse(sigma_method="single_pass") [current]
#             E1 = fit_mvn_bm_inhouse(sigma_method="fisher_ml")   [the port]
#             E2 = Rphylopars::phylopars REML                     [yardstick]
# Run AFTER the port lands; if sigma_method is absent, E1 is skipped and
# reported as such (lets E0/E2 arms run early).

suppressPackageStartupMessages({ library(ape) })
devtools::load_all(".", quiet = TRUE)
t0 <- proc.time()[["elapsed"]]
say <- function(...) { cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
                           ..., "\n", sep = ""); flush.console() }

K <- 4L; MISS <- 0.30; N_REP <- 50L
GRID <- expand.grid(n = c(100L, 300L), lambda = c(1.0, 0.2),
                    sigma_design = c("exch07", "exch03", "hetero"),
                    stringsAsFactors = FALSE)
has_fisher <- "sigma_method" %in% names(formals(fit_mvn_bm_inhouse))
say("fisher_ml available: ", has_fisher)

make_sigma <- function(design) {
  if (design == "exch07") { S <- matrix(0.7, K, K); diag(S) <- 1; return(S) }
  if (design == "exch03") { S <- matrix(0.3, K, K); diag(S) <- 1; return(S) }
  v <- c(0.5, 1, 1.5, 2)
  Rho <- matrix(c(1, .5, -.3, .2,  .5, 1, .2, -.15,
                  -.3, .2, 1, .4,  .2, -.15, .4, 1), K, K)
  S <- diag(sqrt(v)) %*% Rho %*% diag(sqrt(v))
  stopifnot(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) > 0.05)
  S
}

fit_arm <- function(arm, L, tree) {
  if (arm == "E0") return(fit_mvn_bm_inhouse(L = L, tree = tree))
  if (arm == "E1") return(fit_mvn_bm_inhouse(L = L, tree = tree,
                                             sigma_method = "fisher_ml"))
  if (arm == "E1i1") return(fit_mvn_bm_inhouse(L = L, tree = tree,
                              sigma_method = "fisher_ml", max_iter = 1L))
  if (arm == "E1i3") return(fit_mvn_bm_inhouse(L = L, tree = tree,
                              sigma_method = "fisher_ml", max_iter = 3L))
  if (arm == "E0i1") return(fit_mvn_bm_inhouse(L = L, tree = tree,
                              max_iter = 1L))
  Rphylopars::phylopars(data.frame(species = rownames(L), L,
                                   stringsAsFactors = FALSE),
                        tree = tree, model = "BM",
                        phylo_correlated = TRUE, pheno_correlated = TRUE,
                        REML = TRUE)
}

rows <- list()
for (g in seq_len(nrow(GRID))) {
  n <- GRID$n[g]; lam <- GRID$lambda[g]; sd_name <- GRID$sigma_design[g]
  Sigma <- make_sigma(sd_name)
  say(sprintf("cell %d/%d: n=%d lambda=%.1f %s", g, nrow(GRID), n, lam, sd_name))
  for (r in seq_len(N_REP)) {
    seed <- 20260817L + g * 1000L + r
    set.seed(seed)
    tree <- ape::rtree(n)
    V <- ape::vcv(tree); V <- V / max(V)
    R <- lam * stats::cov2cor(V) + (1 - lam) * diag(n)
    cR <- chol(R + 1e-8 * diag(n)); cS <- chol(Sigma)
    L_true <- t(cR) %*% matrix(stats::rnorm(n * K), n, K) %*% cS
    rownames(L_true) <- tree$tip.label
    colnames(L_true) <- paste0("t", 1:K)
    mask <- matrix(stats::runif(n * K) < MISS, n, K)
    L <- L_true; L[mask] <- NA

    # lam-transformed tree for the fit so all arms see the SAME (correct) R:
    # scale internal structure by building a tree whose vcv is R -- simplest:
    # keep the original tree at lam=1; at lam=0.2, transform edge lengths via
    # pigauto's own transform (per-column path uses R directly; the solvers
    # take a tree). Use transform_tree_pagel when available, else skip
    # lam != 1 for tree-taking arms and note it.
    tr_fit <- if (lam == 1) tree else
      tryCatch(pigauto:::transform_tree_pagel(tree, lam),
               error = function(e) NULL)
    if (is.null(tr_fit)) next
    # height-normalize so E2's phylocov is on the simulation's unit scale
    # (E0/E1 cov2cor internally -- unaffected; fixes the E2 frob artifact)
    tr_fit$edge.length <- tr_fit$edge.length / max(ape::node.depth.edgelength(tr_fit))

    for (arm in c("E0", "E0i1", if (has_fisher) c("E1", "E1i1", "E1i3"), "E2")) {
      t1 <- proc.time()[["elapsed"]]
      fit <- tryCatch(suppressWarnings(fit_arm(arm, L, tr_fit)),
                      error = function(e) e)
      wall <- proc.time()[["elapsed"]] - t1
      if (inherits(fit, "error")) {
        rows[[length(rows)+1L]] <- data.frame(cell = g, n = n, lambda = lam,
          sigma_design = sd_name, rep = r, arm = arm, failed = TRUE,
          frob = NA, sign_acc = NA, rmse = NA, cover = NA, wall = wall)
        next
      }
      S_hat <- tryCatch(as.matrix(fit$pars$phylocov), error = function(e) NULL)
      frob <- if (!is.null(S_hat) && all(dim(S_hat) == K))
        norm(S_hat - Sigma, "F") / norm(Sigma, "F") else NA_real_
      off <- upper.tri(Sigma)
      sign_acc <- if (!is.null(S_hat) && all(dim(S_hat) == K))
        mean(sign(S_hat[off]) == sign(Sigma[off])) else NA_real_
      tip <- match(rownames(L), rownames(fit$anc_recon))
      mu <- fit$anc_recon[tip, , drop = FALSE]
      va <- fit$anc_var[tip, , drop = FALSE]
      err <- (mu - L_true)[mask]
      rmse <- sqrt(mean(err^2))
      se <- sqrt(pmax(va[mask], 0))
      cover <- mean(abs(err) <= 1.96 * se)
      rows[[length(rows)+1L]] <- data.frame(cell = g, n = n, lambda = lam,
        sigma_design = sd_name, rep = r, arm = arm, failed = FALSE,
        frob = frob, sign_acc = sign_acc, rmse = rmse, cover = cover,
        wall = wall)
    }
  }
}
d <- do.call(rbind, rows)
agg <- do.call(rbind, lapply(split(d[!d$failed, ],
    d[!d$failed, c("n","lambda","sigma_design","arm")], drop = TRUE),
  function(gg) data.frame(gg[1, c("n","lambda","sigma_design","arm")],
    n_rep = nrow(gg), frob = mean(gg$frob), frob_mcse = sd(gg$frob)/sqrt(nrow(gg)),
    sign_acc = mean(gg$sign_acc), rmse = mean(gg$rmse),
    rmse_mcse = sd(gg$rmse)/sqrt(nrow(gg)), cover = mean(gg$cover),
    cover_mcse = sd(gg$cover)/sqrt(nrow(gg)), wall = mean(gg$wall))))
rownames(agg) <- NULL
agg <- agg[order(agg$lambda, agg$sigma_design, agg$n, agg$arm), ]
say("failures: ", sum(d$failed))
print(agg, digits = 3)
saveRDS(list(raw = d, agg = agg), "dev/sigma_recovery_sim.rds")
say("wrote dev/sigma_recovery_sim.rds")
