# script/rubin_discrete.R
#
# Discrete-trait arms and scores for the freq-vs-BACE Rubin study (script/rubin_study/study.qmd, section
# "Discrete traits"). Assumes campaign_gnn_off_lib.R and rubin_lib.R are sourced (rubin_cell.R does this).
#
# Frequentist discrete arms mirror the continuous freq A / freq B:
#   freq B analogue: one castor Mk fit (castor::hsp_mk_model), then M joint conditional draws of the missing tips at
#     the fitted rates.
#   freq A analogue: for each m, a parametric bootstrap of the rates (simulate the trait on the tree at the fitted Q,
#     keep the observed tips, refit), then one joint conditional draw at the refitted Q.
# Rate models (review of 2026-09-26): the default is the flexible set CASTOR_FLEX (all rates different, ARD, for binary
# and categorical traits; stepwise with every rate different, SRD, for the ordinal trait). Its stationary distribution
# can match the observed class frequencies, so at low phylogenetic signal it falls back to them, as the continuous
# arm's Pagel lambda falls back to the sample mean. Equal rates (ER) cannot: at low lambda its rate runs to the bound
# and it imputes uniform classes, below the observed-frequency floor. Simulation v1's set (ER, ER, SUEDE) is kept as
# CASTOR_V1 for a comparison arm.
# The joint draw (mk_joint_draw) is forward filtering, backward sampling on the tree: exact draws from the joint
# distribution of all missing tips given the observed tips, with the root prior castor's fit uses ("empirical":
# observed class frequencies; at very high fitted rates the draws follow the model's stationary distribution instead).
# Transition probabilities use the matrix exponential (Matrix::expm), which also handles the defective rate matrices
# that boundary fits produce (a stepwise rate at exactly 0). castor's own hsp marginals are exact only for reversible
# models, so the tests compare the joint draws with brute-force marginals.
#
# castor models one trait at a time and has no covariates, so these arms ignore c1, c2 and the driver d1.
# A trait with fewer than two observed classes cannot be fitted (castor stops; BACE leaves it out); every arm, BACE
# included, then imputes the single observed class (fill_degenerate()), and the filled cells are counted. A castor
# failure on one trait leaves that trait missing for the arm (recorded in the arm's diag) and never touches the others.

disc_traits_of <- function(truth) names(truth)[vapply(truth, is.factor, logical(1))]

CASTOR_FLEX <- c(binary = "ARD", categorical = "ARD", ordinal = "SRD")
CASTOR_V1   <- c(binary = "ER",  categorical = "ER",  ordinal = "SUEDE")
.trait_kind <- function(f) if (is.ordered(f)) "ordinal" else if (nlevels(f) == 2L) "binary" else "categorical"

.edge_P <- function(Q, t) {
  lapply(t, function(s) { P <- as.matrix(Matrix::expm(Q * s)); P[P < 0] <- 0; P / rowSums(P) })
}

#' Joint conditional draws of the missing tip states of one Mk trait.
#' @param tree ape phylo; @param tip integer states 1..K in tree$tip.label order, NA = missing
#' @param Q K x K transition-rate matrix; @param root_prior length-K probabilities; @param M number of draws
#' @return M x Ntip integer matrix of complete tip states (observed tips kept). Stops if the observed tips have zero
#'   likelihood under Q (possible for a bootstrap refit with rates at exactly 0), instead of returning NA draws.
mk_joint_draw <- function(tree, tip, Q, root_prior, M) {
  K <- nrow(Q); nt <- length(tree$tip.label); nn <- tree$Nnode
  tr <- ape::reorder.phylo(tree, "postorder")
  P <- .edge_P(Q, tr$edge.length)
  lik <- matrix(1, nt + nn, K)
  obs <- !is.na(tip); lik[which(obs), ] <- 0; lik[cbind(which(obs), tip[obs])] <- 1
  for (i in seq_len(nrow(tr$edge))) {
    p <- tr$edge[i, 1]; c <- tr$edge[i, 2]
    lik[p, ] <- lik[p, ] * as.numeric(P[[i]] %*% lik[c, ])
    mx <- max(lik[p, ])
    if (!is.finite(mx) || mx <= 0) stop("the observed tips have zero likelihood under Q")
    lik[p, ] <- lik[p, ] / mx              # rescale against underflow; sampling ignores the scale
  }
  root <- nt + 1L
  draw_cat <- function(pr) {                 # pr: M x K unnormalised probabilities, one row per draw
    rs <- rowSums(pr)
    if (any(!is.finite(rs) | rs <= 0)) stop("the observed tips have zero likelihood under Q")
    cp <- t(apply(pr / rs, 1, cumsum))
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
                            root_prior = "empirical", include_likelihoods = FALSE, Nthreads = 1)
  if (!isTRUE(h$success)) stop("castor hsp_mk_model failed: ", h$error %||% "unknown")
  h$transition_matrix
}
`%||%` <- function(x, y) if (is.null(x)) y else x

#' castor arms for every discrete trait of a cell.
#' @param cell a make_cell() result; @param M draws; @param proper TRUE = freq A analogue (bootstrap rates)
#' @param base_sets optional list of M data.frames (e.g. the freq A continuous draws) whose discrete columns are
#'   filled; default copies of df_miss
#' @param models named rate models for binary, categorical and ordinal traits (CASTOR_FLEX or CASTOR_V1)
#' @return list(datasets, diag = per-trait data.frame: rate model, degenerate flag, largest exit rate of the fit and
#'   median over the refits used, bootstrap re-simulations (< 2 observed classes), failed refits, rejected refits
#'   (observed tips impossible under Q*), draws that fell back to the fitted Q, error)
mi_castor <- function(cell, M, proper, base_sets = NULL, models = CASTOR_FLEX) {
  tree <- cell$tree; df <- cell$df_miss
  sets <- base_sets %||% rep(list(df), M)
  dg <- list()
  for (v in disc_traits_of(cell$truth)) {
    f <- df[[v]]; lev <- levels(f); K <- length(lev)
    tip <- as.integer(f)[match(tree$tip.label, rownames(df))]
    rm <- models[[.trait_kind(f)]]
    miss_sp <- tree$tip.label[is.na(tip)]
    row <- data.frame(trait = v, rate_model = rm, degenerate = FALSE, q_hat = NA_real_, q_star_med = NA_real_,
                      n_resim = 0L, n_refit_fail = 0L, n_reject = 0L, n_fallback = 0L, error = NA_character_)
    if (length(unique(tip[!is.na(tip)])) < 2L) {   # castor cannot fit; leave NA for fill_degenerate()
      row$degenerate <- TRUE; dg[[v]] <- row; next
    }
    res <- tryCatch({
      Q <- .castor_fit(tree, tip, K, rm); prior <- .empirical_prior(tip, K)
      n_resim <- 0L; n_fail <- 0L; n_reject <- 0L; n_fallback <- 0L
      if (!proper) {
        draws <- mk_joint_draw(tree, tip, Q, prior, M); q_star <- rep(max(-diag(Q)), M)
      } else {
        draws <- matrix(NA_integer_, M, length(tip)); q_star <- numeric(M)
        obs <- !is.na(tip)
        for (m in seq_len(M)) {
          dm <- NULL; tries <- 0L
          # A bootstrap Q* is kept only if the observed tips are possible under it (a refit can put rates at exactly
          # 0); otherwise it is redrawn. After 20 failed tries the fitted Q is used. All three events are counted.
          while (is.null(dm) && tries < 20L) {
            tries <- tries + 1L
            repeat {                           # a bootstrap sample needs >= 2 observed classes to refit
              sim <- castor::simulate_mk_model(tree, Q, root_probabilities = prior, include_nodes = FALSE)$tip_states
              tip_b <- ifelse(obs, sim, NA_integer_)
              if (length(unique(tip_b[obs])) >= 2L || n_resim > 50L * M) break
              n_resim <- n_resim + 1L
            }
            Qm <- tryCatch(.castor_fit(tree, tip_b, K, rm), error = function(e) NULL)
            if (is.null(Qm)) { n_fail <- n_fail + 1L; next }
            dm <- tryCatch(mk_joint_draw(tree, tip, Qm, prior, 1L), error = function(e) NULL)
            if (is.null(dm)) n_reject <- n_reject + 1L
          }
          if (is.null(dm)) { Qm <- Q; dm <- mk_joint_draw(tree, tip, Q, prior, 1L); n_fallback <- n_fallback + 1L }
          draws[m, ] <- dm; q_star[m] <- max(-diag(Qm))
        }
      }
      list(draws = draws, q_hat = max(-diag(Q)), q_star = q_star, n_resim = n_resim, n_fail = n_fail,
           n_reject = n_reject, n_fallback = n_fallback)
    }, error = function(e) e)
    if (inherits(res, "error")) { row$error <- conditionMessage(res); dg[[v]] <- row; next }
    draws <- res$draws; colnames(draws) <- tree$tip.label
    for (m in seq_len(M)) {
      x <- as.character(sets[[m]][[v]]); x[match(miss_sp, rownames(df))] <- lev[draws[m, miss_sp]]
      sets[[m]][[v]] <- factor(x, levels = lev, ordered = is.ordered(f))
    }
    row$q_hat <- res$q_hat; row$q_star_med <- stats::median(res$q_star)
    row$n_resim <- res$n_resim; row$n_refit_fail <- res$n_fail; row$n_reject <- res$n_reject
    row$n_fallback <- res$n_fallback
    dg[[v]] <- row
  }
  list(datasets = sets, diag = do.call(rbind, dg))
}

#' A trait with a single observed class: impute that class in every missing cell, for every arm alike.
#' @return list(sets, n_filled = number of missing cells filled per trait, per dataset)
fill_degenerate <- function(sets, df_miss, traits) {
  n_filled <- setNames(integer(length(traits)), traits)
  for (v in traits) {
    cls <- unique(as.character(df_miss[[v]][!is.na(df_miss[[v]])]))
    if (length(cls) != 1L) next
    n_filled[v] <- sum(is.na(df_miss[[v]]))
    for (m in seq_along(sets)) {
      x <- as.character(sets[[m]][[v]]); x[is.na(x)] <- cls
      sets[[m]][[v]] <- factor(x, levels = levels(df_miss[[v]]), ordered = is.ordered(df_miss[[v]]))
    }
  }
  list(sets = sets, n_filled = n_filled)
}

#' Per-value scores of one arm on the discrete traits, from its M draws only (every arm on the same basis).
#' accuracy: of the modal class; brier: multi-class Brier score of the draw frequencies (sum over classes, mean over
#' cells); ece: top-label expected calibration error, 10 equal-width bins of the modal frequency;
#' set_coverage / set_size: the smallest set of classes holding >= 95% of the draws; mae_class: mean absolute class
#' distance of the modal class (ordinal only); frac_unanimous: all M draws agree. Ties (for the modal class and at the
#' edge of the 95% set) are counted fractionally: the expected value under a random tie-break.
#' A trait with one realised class in the complete data has nothing to predict and is skipped (its absence is counted
#' at aggregation). NA draws are dropped cell by cell; cells whose draws are all NA are dropped and counted.
score_discrete <- function(arm, sets, truth, mask, traits = disc_traits_of(truth)) {
  tol <- 1e-12; rows <- list()
  for (v in traits) {
    idx <- which(mask[, v]); if (!length(idx)) next
    if (nlevels(droplevels(truth[[v]])) < 2L) next
    lev <- levels(truth[[v]]); K <- length(lev)
    D <- vapply(sets, function(s) as.integer(factor(as.character(s[[v]][idx]), levels = lev)), integer(length(idx)))
    D <- matrix(D, nrow = length(idx))                 # cells x M
    n_na <- sum(is.na(D))
    keep <- rowSums(!is.na(D)) > 0L; n_all_na <- sum(!keep)
    if (!any(keep)) next
    idx <- idx[keep]; D <- D[keep, , drop = FALSE]
    Pm <- t(apply(D, 1, function(d) tabulate(d[!is.na(d)], K) / sum(!is.na(d))))
    tr <- as.integer(truth[[v]][idx]); Y <- diag(K)[tr, , drop = FALSE]
    top <- apply(Pm, 1, max); is_top <- abs(Pm - top) <= tol
    acc <- rowSums(is_top & Y == 1) / rowSums(is_top)
    in_set <- vapply(seq_along(idx), function(i) {
      p <- Pm[i, ]; ps <- sort(p, decreasing = TRUE)
      k <- which(cumsum(ps) >= 0.95 - 1e-9)[1]; thr <- ps[k]
      n_gt <- sum(p > thr + tol); n_eq <- sum(abs(p - thr) <= tol)
      pt <- p[tr[i]]
      c(if (pt > thr + tol) 1 else if (abs(pt - thr) <= tol) (k - n_gt) / n_eq else 0, k)
    }, numeric(2))
    bins <- cut(top, seq(0, 1, 0.1), include.lowest = TRUE)
    ece <- sum(tapply(seq_along(idx), bins, function(j) length(j) * abs(mean(acc[j]) - mean(top[j]))), na.rm = TRUE) /
      length(idx)
    mae <- if (is.ordered(truth[[v]])) mean(vapply(seq_along(idx), function(i) mean(abs(which(is_top[i, ]) - tr[i])), numeric(1))) else NA_real_
    rows[[v]] <- data.frame(arm = arm, trait = v, n_cells = length(idx), accuracy = mean(acc),
                            brier = mean(rowSums((Pm - Y)^2)), ece = ece,
                            set_coverage = mean(in_set[1, ]), set_size = mean(in_set[2, ]),
                            mae_class = mae, frac_unanimous = mean(top == 1), n_na = n_na, n_cells_all_na = n_all_na)
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

#' GLS slope of y on x at a fixed Pagel lambda (eigenbasis of pagel_eigen()).
.gls_slope_at <- function(y, x, eig, lam) {
  X <- cbind(1, x); Xs <- crossprod(eig$U, X); ys <- as.numeric(crossprod(eig$U, y))
  w <- 1 / (lam * eig$d + (1 - lam))
  solve(crossprod(Xs, w * Xs), crossprod(Xs, w * ys))[2, 1]
}

#' Rubin-pooled downstream estimand on the discrete trait: PGLS slope of c1 on bin (bin = "yes" coded 1).
#' Two targets are stored for coverage. The population target (filled at aggregation from the Monte Carlo of the
#' complete-data estimate) is not what any per-dataset interval estimates at lambda = 1, rho = 0.5, where the
#' dataset-specific value varies more than the standard error (review of 2026-09-26: complete-data coverage 0.31 to
#' 0.56). The per-dataset target, target_cond = rho times the GLS slope (at the true lambda) of the bin liability L5
#' on bin, is E[slope | tree, L5] with lambda known; complete-data coverage of it is about 0.96 in every cell.
#' Undefined (NULL) when bin has one class in the complete data, or is constant in every imputed dataset (one
#' observed class, filled by fill_degenerate()).
#' @param L5 the bin liability (make_cell()$L[, 5]) in rownames(truth) order; @param lambda,rho the cell's values
score_discrete_estimand <- function(arm, sets, truth, tree, eig, n, L5 = NULL, lambda = NA_real_, rho = NA_real_) {
  if (nlevels(droplevels(truth$bin)) < 2L) return(NULL)
  if (all(vapply(sets, function(s) nlevels(droplevels(s$bin)) < 2L && !anyNA(s$bin), logical(1)))) return(NULL)
  if (is.null(eig)) eig <- pagel_eigen(tree, rownames(truth))
  sl <- lapply(sets, est_pgls_xy, tree = tree, y = "c1", x = "bin", eig = eig)
  ok <- vapply(sl, function(s) is.finite(s$estimate) && is.finite(s$variance), logical(1))
  if (sum(ok) < 2L) stop(sprintf("%s: fewer than 2 analysable imputations for c1 ~ bin", arm))
  ps <- rubin_pool(vapply(sl[ok], `[[`, numeric(1), "estimate"), vapply(sl[ok], `[[`, numeric(1), "variance"),
                   df_com = n - 2)
  ref <- est_pgls_xy(truth, tree, y = "c1", x = "bin", eig = eig)
  # the complete-data interval (t with n - 2 df), the coverage baseline, as the continuous "complete" row
  ref_se <- sqrt(ref$variance); tq <- stats::qt(0.975, n - 2)
  target_cond <- if (!is.null(L5) && is.finite(lambda) && is.finite(rho))
    rho * .gls_slope_at(L5, as.numeric(truth$bin == levels(truth$bin)[2L]), eig, lambda) else NA_real_
  data.frame(arm = arm, estimand = "slope_c1_bin", estimate = ps$estimate, se = ps$se, lower = ps$lower,
             upper = ps$upper, df = ps$df, fmi = ps$fmi, m_ok = sum(ok), complete_data = ref$estimate,
             complete_se = ref_se, complete_lower = ref$estimate - tq * ref_se,
             complete_upper = ref$estimate + tq * ref_se, target_cond = target_cond)
}
