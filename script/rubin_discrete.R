# script/rubin_discrete.R
#
# Discrete-trait arms and scores for the freq-vs-BACE Rubin study (script/rubin_study/study.qmd, section
# "Discrete traits"). Assumes campaign_gnn_off_lib.R and rubin_lib.R are sourced (rubin_cell.R does this).
#
# Frequentist discrete arms mirror the continuous freq A / freq B:
#   freq B analogue: one castor Mk fit (castor::hsp_mk_model; ER for binary and categorical, SUEDE for
#     ordinal, as simulation v1's run_freq()), then M joint conditional draws of the missing tips at the
#     fitted rates.
#   freq A analogue: for each m, a parametric bootstrap of the rates (simulate the trait on the tree at the
#     fitted Q, keep the observed tips, refit), then one joint conditional draw at the refitted Q.
# The joint draw (mk_joint_draw) is forward filtering, backward sampling on the tree: exact draws from the
# joint distribution of all missing tips given the observed tips, under the same root prior castor uses
# ("empirical": observed class frequencies). castor's hsp gives per-tip marginals only; the per-tip
# marginals of the joint draws equal castor's (test-discrete.R).
#
# castor models one trait at a time and has no covariates, so these arms ignore c1, c2 and the driver d1.
# A trait with fewer than two observed classes cannot be fitted (castor stops; BACE leaves it out); every arm,
# BACE included, then imputes the single observed class (fill_degenerate()), and the fill is counted.

disc_traits_of <- function(truth) names(truth)[vapply(truth, is.factor, logical(1))]

.edge_P <- function(Q, t) {
  e <- eigen(Q)
  if (is.complex(e$values) || any(abs(Im(e$vectors)) > 0)) return(lapply(t, function(s) as.matrix(Matrix::expm(Q * s))))
  V <- e$vectors; Vi <- solve(V); d <- e$values
  lapply(t, function(s) { P <- V %*% (exp(d * s) * Vi); P[P < 0] <- 0; P / rowSums(P) })
}

#' Joint conditional draws of the missing tip states of one Mk trait.
#' @param tree ape phylo; @param tip integer states 1..K in tree$tip.label order, NA = missing
#' @param Q K x K transition-rate matrix; @param root_prior length-K probabilities; @param M number of draws
#' @return M x Ntip integer matrix of complete tip states (observed tips kept)
mk_joint_draw <- function(tree, tip, Q, root_prior, M) {
  K <- nrow(Q); nt <- length(tree$tip.label); nn <- tree$Nnode
  tr <- ape::reorder.phylo(tree, "postorder")
  P <- .edge_P(Q, tr$edge.length)
  lik <- matrix(1, nt + nn, K)
  obs <- !is.na(tip); lik[which(obs), ] <- 0; lik[cbind(which(obs), tip[obs])] <- 1
  for (i in seq_len(nrow(tr$edge))) {
    p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
    lik[p, ] <- lik[p, ] * as.numeric(P[[i]] %*% lik[c, ])
    lik[p, ] <- lik[p, ] / max(lik[p, ])   # rescale against underflow; sampling ignores the scale
  }
  root <- nt + 1L
  draw_cat <- function(pr) {                 # pr: M x K unnormalised probabilities, one row per draw
    cp <- t(apply(pr / rowSums(pr), 1, cumsum))
    1L + rowSums(cp < stats::runif(nrow(pr)))
  }
  st <- matrix(NA_integer_, M, nt + nn)
  st[, root] <- draw_cat(matrix(root_prior * lik[root, ], M, K, byrow = TRUE))
  for (i in rev(seq_len(nrow(tr$edge)))) {   # reverse postorder visits every parent before its children
    p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
    st[, c] <- draw_cat(P[[i]][st[, p], , drop = FALSE] * matrix(lik[c, ], M, K, byrow = TRUE))
  }
  st[, seq_len(nt), drop = FALSE]
}

.empirical_prior <- function(tip, K) { f <- tabulate(tip[!is.na(tip)], K); f / sum(f) }

.castor_fit <- function(tree, tip, K, rate_model) {
  h <- castor::hsp_mk_model(tree, tip_states = tip, Nstates = K, rate_model = rate_model, Ntrials = 3,
                            root_prior = "empirical", Nthreads = 1)
  if (!isTRUE(h$success)) stop("castor hsp_mk_model failed: ", h$error %||% "unknown")
  h$transition_matrix
}
`%||%` <- function(x, y) if (is.null(x)) y else x

#' castor arms for every discrete trait of a cell.
#' @param cell a make_cell() result; @param M draws; @param proper TRUE = freq A analogue (bootstrap rates)
#' @param base_sets optional list of M data.frames (e.g. the freq A continuous draws) whose discrete columns are
#'   filled; default copies of df_miss
#' @return list(datasets, diag = per-trait data.frame of fitted rate, degenerate flag, bootstrap retries)
mi_castor <- function(cell, M, proper, base_sets = NULL) {
  tree <- cell$tree; df <- cell$df_miss
  sets <- base_sets %||% rep(list(df), M)
  dg <- list()
  for (v in disc_traits_of(cell$truth)) {
    f <- df[[v]]; lev <- levels(f); K <- length(lev)
    tip <- as.integer(f)[match(tree$tip.label, rownames(df))]
    rm <- if (is.ordered(f)) "SUEDE" else "ER"
    miss_sp <- tree$tip.label[is.na(tip)]
    n_cls <- length(unique(tip[!is.na(tip)]))
    if (n_cls < 2L) {                        # castor cannot fit; leave NA for fill_degenerate()
      dg[[v]] <- data.frame(trait = v, degenerate = TRUE, q_hat = NA, q_star_med = NA, n_resim = 0L)
      next
    }
    Q <- .castor_fit(tree, tip, K, rm); prior <- .empirical_prior(tip, K)
    if (!proper) {
      draws <- mk_joint_draw(tree, tip, Q, prior, M); q_star <- rep(max(-diag(Q)), M); n_resim <- 0L
    } else {
      draws <- matrix(NA_integer_, M, length(tip)); q_star <- numeric(M); n_resim <- 0L
      obs <- !is.na(tip)
      for (m in seq_len(M)) {
        repeat {                               # a bootstrap sample needs >= 2 observed classes to refit
          sim <- castor::simulate_mk_model(tree, Q, root_probabilities = prior, include_nodes = FALSE)$tip_states
          tip_b <- ifelse(obs, sim, NA_integer_)
          if (length(unique(tip_b[obs])) >= 2L || n_resim > 50L * M) break
          n_resim <- n_resim + 1L
        }
        Qm <- tryCatch(.castor_fit(tree, tip_b, K, rm), error = function(e) Q)
        draws[m, ] <- mk_joint_draw(tree, tip, Qm, prior, 1L); q_star[m] <- max(-diag(Qm))
      }
    }
    colnames(draws) <- tree$tip.label
    for (m in seq_len(M)) {
      x <- as.character(sets[[m]][[v]]); x[match(miss_sp, rownames(df))] <- lev[draws[m, miss_sp]]
      sets[[m]][[v]] <- factor(x, levels = lev, ordered = is.ordered(f))
    }
    dg[[v]] <- data.frame(trait = v, degenerate = FALSE, q_hat = max(-diag(Q)), q_star_med = stats::median(q_star),
                            n_resim = n_resim)
  }
  list(datasets = sets, diag = do.call(rbind, dg))
}

#' A trait with a single observed class: impute that class in every missing cell, for every arm alike.
#' @return list(sets, n_filled = number of cells filled per trait)
fill_degenerate <- function(sets, df_miss, traits) {
  n_filled <- setNames(integer(length(traits)), traits)
  for (v in traits) {
    cls <- unique(as.character(df_miss[[v]][!is.na(df_miss[[v]])]))
    if (length(cls) != 1L) next
    for (m in seq_along(sets)) {
      x <- as.character(sets[[m]][[v]]); na <- is.na(x); n_filled[v] <- n_filled[v] + sum(na)
      x[na] <- cls; sets[[m]][[v]] <- factor(x, levels = levels(df_miss[[v]]), ordered = is.ordered(df_miss[[v]]))
    }
  }
  list(sets = sets, n_filled = n_filled)
}

#' Per-value scores of one arm on the discrete traits, from its M draws only (every arm on the same basis).
#' accuracy: of the modal class, ties counted fractionally (the expected accuracy of a random tie-break);
#' brier: multi-class Brier score of the draw frequencies (sum over classes, mean over cells);
#' ece: top-label expected calibration error, 10 equal-width bins of the modal frequency;
#' set_coverage / set_size: the smallest set of classes holding >= 95% of the draws;
#' mae_class: mean absolute class distance of the modal class (ordinal only); frac_unanimous: all M draws agree.
#' A trait with one class in the complete data (the DGP keeps only realised levels) has nothing to predict and is
#' skipped; its absence from the output is counted at aggregation.
score_discrete <- function(arm, sets, truth, mask, traits = disc_traits_of(truth)) {
  rows <- list()
  for (v in traits) {
    idx <- which(mask[, v]); if (!length(idx)) next
    lev <- levels(truth[[v]]); K <- length(lev)
    if (K < 2L) next
    D <- vapply(sets, function(s) as.integer(factor(as.character(s[[v]][idx]), levels = lev)), integer(length(idx)))
    D <- matrix(D, nrow = length(idx))                 # cells x M
    n_na <- sum(is.na(D))
    if (n_na == length(D)) next
    Pm <- t(apply(D, 1, function(d) tabulate(d[!is.na(d)], K) / sum(!is.na(d))))
    tr <- as.integer(truth[[v]][idx]); Y <- diag(K)[tr, , drop = FALSE]
    top <- apply(Pm, 1, max); is_top <- Pm == top
    acc <- rowSums(is_top & Y == 1) / rowSums(is_top)
    ord_p <- t(apply(Pm, 1, order, decreasing = TRUE))
    in_set <- vapply(seq_along(idx), function(i) {
      o <- ord_p[i, ]; k <- which(cumsum(Pm[i, o]) >= 0.95 - 1e-9)[1]; c(tr[i] %in% o[seq_len(k)], k)
    }, numeric(2))
    bins <- cut(top, seq(0, 1, 0.1), include.lowest = TRUE)
    ece <- sum(tapply(seq_along(idx), bins, function(j) length(j) * abs(mean(acc[j]) - mean(top[j]))), na.rm = TRUE) /
      length(idx)
    mode_cls <- apply(Pm, 1, which.max)
    rows[[v]] <- data.frame(arm = arm, trait = v, n_cells = length(idx), accuracy = mean(acc),
                            brier = mean(rowSums((Pm - Y)^2)), ece = ece,
                            set_coverage = mean(in_set[1, ]), set_size = mean(in_set[2, ]),
                            mae_class = if (is.ordered(truth[[v]])) mean(abs(mode_cls - tr)) else NA_real_,
                            frac_unanimous = mean(top == 1), n_na = n_na)
  }
  do.call(rbind, rows)
}

#' PGLS slope of y on x (Pagel's lambda by REML, the eigenbasis estimator of est_pgls_slope_fast) for any two
#' columns; x may be a two-level factor, coded 0/1 (second level = 1).
est_pgls_xy <- function(df, tree, y, x, eig = NULL) {
  d2 <- data.frame(row.names = rownames(df), c1 = df[[x]], c2 = df[[y]])
  if (is.factor(d2$c1)) d2$c1 <- as.numeric(d2$c1 == levels(d2$c1)[2L])
  est_pgls_slope_fast(d2, tree, eig = eig)
}

#' Rubin-pooled downstream estimand on the discrete trait: PGLS slope of c1 on bin (bin = "yes" coded 1).
#' The truth is not rho; it is filled at aggregation from a Monte Carlo of the complete-data estimate (0 at rho = 0).
#' Undefined (NULL) when bin has one class in the complete data.
score_discrete_estimand <- function(arm, sets, truth, tree, eig, n) {
  if (nlevels(droplevels(truth$bin)) < 2L) return(NULL)
  sl <- lapply(sets, est_pgls_xy, tree = tree, y = "c1", x = "bin", eig = eig)
  ok <- vapply(sl, function(s) is.finite(s$estimate) && is.finite(s$variance), logical(1))
  if (sum(ok) < 2L) stop(sprintf("%s: fewer than 2 analysable imputations for c1 ~ bin", arm))
  ps <- rubin_pool(vapply(sl[ok], `[[`, numeric(1), "estimate"), vapply(sl[ok], `[[`, numeric(1), "variance"),
                   df_com = n - 2)
  ref <- est_pgls_xy(truth, tree, y = "c1", x = "bin", eig = eig)
  data.frame(arm = arm, estimand = "slope_c1_bin", estimate = ps$estimate, se = ps$se, lower = ps$lower,
             upper = ps$upper, df = ps$df, fmi = ps$fmi, m_ok = sum(ok), complete_data = ref$estimate)
}
