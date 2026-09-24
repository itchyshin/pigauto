# script/mi_gls/dgp_v2.R
#
# Shared data-generating process for script/mi_gls/01_cell_v2.R and
# script/mi_gls/gate_convergence.R (arc/mi-posterior). Regimes 1-16 use the
# ORIGINAL tree-Pagel-transform DGP from script/mi_gls/01_cell.R verbatim
# (same seed formula -> byte-identical simulated data, so v1 and v2 outputs
# join by regime_id + rep). Regimes 17-24 use the Sigma_P / Sigma_E
# Kronecker DGP of docs/dev-log/mi-posterior/design.md section 1, over the
# ORIGINAL tree (no Pagel tree-transform -- per-trait lambda is built into
# Sigma_P / Sigma_E directly).
#
# Caller must already have run devtools::load_all() (for
# transform_tree_pagel()) and library(ape) before sourcing this file.

source(file.path("script", "mi_gls", "regimes.R"))

.mi_gls_v2_seed <- function(regime_id, rep_i) 20260923L + regime_id * 10000L + rep_i

# ---- missingness masks (verbatim from script/mi_gls/01_cell.R) --------------
mar_phylo_mask <- function(tree, m_miss = 0.3) {
  n_tip <- ape::Ntip(tree)
  nodes <- (n_tip + 2L):(n_tip + tree$Nnode)
  sizes <- vapply(nodes, function(nd)
    length(ape::extract.clade(tree, nd)$tip.label), integer(1))
  cand <- nodes[sizes >= floor(0.15 * n_tip) & sizes <= ceiling(0.35 * n_tip)]
  picked <- if (length(cand) >= 2L) sample(cand, 2L) else
    nodes[order(abs(sizes - 0.25 * n_tip))[1:2]]
  in_clade <- rep(FALSE, n_tip); names(in_clade) <- tree$tip.label
  for (nd in picked) in_clade[ape::extract.clade(tree, nd)$tip.label] <- TRUE
  pvec <- ifelse(in_clade[tree$tip.label], 7, 1)
  pvec <- pvec * (m_miss * n_tip) / sum(pvec)
  pvec <- pmin(pvec, 0.95)
  mask <- stats::runif(n_tip) < pvec
  if (sum(!mask) < 20L) {
    keep <- sample(which(mask), sum(mask) - (n_tip - 20L))
    mask[keep] <- FALSE
  }
  mask
}
mcar_mask <- function(n_tip, m_miss = 0.3) stats::runif(n_tip) < m_miss

# ---- simulate one (regime_id, rep_i) cell ------------------------------------
# Returns: regime_id, rep, seed, regime (1-row data.frame from `regimes`),
# tree, truth (fully observed 2-col data.frame, rownames = tip labels), df
# (truth with regime's missingness applied), mask_x, mask_y (logical,
# length n), true_beta (population regression coefficient of y on x; see
# script/mi_gls/regimes.R for the 17-24 caveat).
simulate_regime_cell <- function(regime_id, rep_i) {
  reg <- regimes[regimes$regime_id == regime_id, ]
  if (nrow(reg) != 1L) stop("unknown regime_id: ", regime_id, call. = FALSE)
  seed <- .mi_gls_v2_seed(regime_id, rep_i)
  set.seed(seed)

  n <- reg$n
  tree <- ape::rtree(n)

  if (regime_id <= 16L) {
    # ---- verbatim from script/mi_gls/01_cell.R's "---- DGP ----" block ----
    sim_tree <- if (reg$lambda == 1) tree else transform_tree_pagel(tree, reg$lambda)
    V_sim <- ape::vcv(sim_tree); V_sim <- V_sim / max(V_sim)
    Lc  <- chol(V_sim + 1e-8 * diag(n))
    Sig <- matrix(reg$rho, 2, 2); diag(Sig) <- 1
    Z <- t(Lc) %*% matrix(stats::rnorm(n * 2), n, 2) %*% chol(Sig)
    true_beta <- reg$rho
  } else {
    # ---- Sigma_P / Sigma_E Kronecker DGP (design.md section 1) ----
    sg <- build_sigma_pair(reg$lambda_x, reg$lambda_y, reg$corr_phylo, reg$corr_resid)
    R  <- stats::cov2cor(ape::vcv(tree))                 # R = cov2cor(vcv(tree)) per design.md
    Lr <- t(chol(R + 1e-8 * diag(n)))
    A  <- Lr %*% matrix(stats::rnorm(n * 2), n, 2) %*% chol(sg$Sigma_P)   # phylogenetic effect
    E  <- matrix(stats::rnorm(n * 2), n, 2) %*% chol(sg$Sigma_E)          # i.i.d. residual
    Z  <- A + E
    true_beta <- reg$true_beta_pop
  }
  truth <- data.frame(row.names = tree$tip.label, x = Z[, 1], y = Z[, 2])

  make_mask <- function() {
    if (identical(reg$mechanism, "MCAR")) mcar_mask(n) else mar_phylo_mask(tree)
  }
  mask_x <- make_mask()
  mask_y <- if (identical(reg$missing, "both")) make_mask() else rep(FALSE, n)

  df <- truth
  df$x[mask_x] <- NA
  df$y[mask_y] <- NA

  list(regime_id = regime_id, rep = rep_i, seed = seed, regime = reg,
      tree = tree, truth = truth, df = df, mask_x = mask_x, mask_y = mask_y,
      true_beta = true_beta)
}
