# Discrete arms (script/rubin_discrete.R): the joint Mk draw must reproduce castor's per-tip marginals and keep
# observed tips; the castor arms fill only missing cells; the scores and the c1 ~ bin estimand behave on known inputs.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages({
  source(file.path(root, "script", "campaign_gnn_off_lib.R")); source(file.path(root, "script", "rubin_lib.R"))
  source(file.path(root, "script", "rubin_discrete.R"))
})
old_kind <- RNGkind()[1]; RNGkind("Mersenne-Twister"); on.exit(RNGkind(old_kind), add = TRUE)

set.seed(11)
tree <- ape::rcoal(40); tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
# exact per-tip marginals by brute force: the full pruning likelihood with the tip fixed to each state
brute_marginal <- function(tip, Q, prior, i) {
  K <- nrow(Q); tr <- ape::reorder.phylo(tree, "postorder"); P <- .edge_P(Q, tr$edge.length)
  vapply(seq_len(K), function(k) {
    tp <- tip; tp[i] <- k; lik <- matrix(1, 79, K); o <- !is.na(tp); lik[which(o), ] <- 0; lik[cbind(which(o), tp[o])] <- 1
    for (e in seq_len(nrow(tr$edge))) lik[tr$edge[e, 1], ] <- lik[tr$edge[e, 1], ] * as.numeric(P[[e]] %*% lik[tr$edge[e, 2], ])
    sum(prior * lik[41, ])
  }, numeric(1)) |> (\(b) b / sum(b))()
}
check_marginals <- function(K, rate_model, q) {
  Q <- castor::get_random_mk_transition_matrix(K, rate_model = rate_model, max_rate = q)
  tip <- castor::simulate_mk_model(tree, Q)$tip_states
  tip[sample(40, 12)] <- NA
  prior <- .empirical_prior(tip, K)
  d <- mk_joint_draw(tree, tip, Q, prior, 6000L)
  emp <- t(vapply(seq_len(40), function(i) tabulate(d[, i], K) / nrow(d), numeric(K)))
  miss <- which(is.na(tip))
  exact <- t(vapply(miss, function(i) brute_marginal(tip, Q, prior, i), numeric(K)))
  h <- castor::hsp_mk_model(tree, tip, Nstates = K, transition_matrix = Q, root_prior = "empirical")
  list(maxdiff = max(abs(emp[miss, , drop = FALSE] - exact)),
       castor_diff = max(abs(emp - h$likelihoods[1:40, , drop = FALSE])),
       kept = all(d[, !is.na(tip)] == rep(tip[!is.na(tip)], each = nrow(d))))
}

# castor's hsp marginals are exact only for reversible models under its own root handling (rerooting method);
# they agree with the joint draws for ER but not for SUEDE, so the reference here is the brute-force marginal.
testthat::test_that("joint Mk draws have the exact per-tip marginals and keep observed tips (ER 2, ER 3, SUEDE 4)", {
  for (cfg in list(list(2, "ER", 2), list(3, "ER", 3), list(4, "SUEDE", 4))) {
    r <- check_marginals(cfg[[1]], cfg[[2]], cfg[[3]])
    message(sprintf("[test-discrete] K=%d %s max |marginal diff| exact %.3f, castor %.3f", cfg[[1]], cfg[[2]],
                    r$maxdiff, r$castor_diff))
    testthat::expect_lt(r$maxdiff, 0.03)   # 6000 draws: Monte Carlo SE <= 0.0065
    if (cfg[[2]] == "ER") testthat::expect_lt(r$castor_diff, 0.03)
    testthat::expect_true(r$kept)
  }
})

cell <- make_cell("types_mixed", 60L, 7L, lambda = 0.7, rho = 0.5, thresholds = "fixed", driver = TRUE)
cell$df_miss <- cell$truth; cell$df_miss[cell$mask] <- NA
for (v in names(cell$truth)) { x <- cell$truth[[v]]; x[cell$mask[, v]] <- NA; cell$df_miss[[v]] <- x }

testthat::test_that("castor arms fill every missing discrete cell and leave observed cells and continuous traits alone", {
  for (cfg in list(list(FALSE, CASTOR_FLEX), list(TRUE, CASTOR_FLEX), list(TRUE, CASTOR_V1))) {
    proper <- cfg[[1]]
    a <- mi_castor(cell, 4L, proper = proper, models = cfg[[2]])
    testthat::expect_true(all(is.na(a$diag$error)))
    testthat::expect_length(a$datasets, 4L)
    for (s in a$datasets) for (v in disc_traits_of(cell$truth)) {
      testthat::expect_false(anyNA(s[[v]]))
      testthat::expect_identical(as.character(s[[v]][!cell$mask[, v]]), as.character(cell$truth[[v]][!cell$mask[, v]]))
      testthat::expect_identical(levels(s[[v]]), levels(cell$truth[[v]]))
    }
    testthat::expect_identical(a$datasets[[1]]$c1, cell$df_miss$c1)
  }
})

testthat::test_that("fill_degenerate imputes a single observed class and counts the cells", {
  df <- data.frame(b = factor(c("yes", "yes", NA, NA), levels = c("no", "yes")))
  f <- fill_degenerate(list(df, df), df, "b")
  testthat::expect_true(all(f$sets[[2]]$b == "yes")); testthat::expect_equal(f$n_filled[["b"]], 2L)   # cells, not cell-draws
})

testthat::test_that("score_discrete: perfect draws score 1/0/1/1; a 50:50 binary tie counts half", {
  s <- score_discrete("x", rep(list(cell$truth), 20), cell$truth, cell$mask)
  testthat::expect_equal(s$accuracy, rep(1, 3)); testthat::expect_equal(s$brier, rep(0, 3))
  testthat::expect_equal(s$set_coverage, rep(1, 3)); testthat::expect_equal(s$set_size, rep(1, 3))
  flip <- cell$truth; flip$bin <- factor(ifelse(flip$bin == "yes", "no", "yes"), levels = levels(flip$bin))
  s2 <- score_discrete("x", c(rep(list(cell$truth), 10), rep(list(flip), 10)), cell$truth, cell$mask, "bin")
  testthat::expect_equal(s2$accuracy, 0.5); testthat::expect_equal(s2$brier, 0.5); testthat::expect_equal(s2$set_size, 2)
})

testthat::test_that("est_pgls_xy codes a two-level factor as 0/1 and matches est_pgls_slope_fast", {
  d <- data.frame(row.names = rownames(cell$truth), c1 = as.numeric(cell$truth$bin == "yes"), c2 = cell$truth$c1)
  testthat::expect_equal(est_pgls_xy(cell$truth, cell$tree, y = "c1", x = "bin")$estimate,
                         est_pgls_slope_fast(d, cell$tree)$estimate)
})

testthat::test_that("a trait with one class in the complete data is skipped by both scores", {
  tr <- cell$truth; tr$bin <- factor(rep("no", nrow(tr)))
  s <- score_discrete("x", rep(list(tr), 3), tr, cell$mask)
  testthat::expect_false("bin" %in% s$trait)
  testthat::expect_null(score_discrete_estimand("x", rep(list(tr), 3), tr, cell$tree, NULL, nrow(tr)))
})

testthat::test_that("the c1 ~ bin row carries the complete-data interval (t, n - 2 df)", {
  e <- score_discrete_estimand("x", rep(list(cell$truth), 3), cell$truth, cell$tree, NULL, nrow(cell$truth))
  testthat::expect_equal(e$complete_upper - e$complete_data, stats::qt(0.975, nrow(cell$truth) - 2) * e$complete_se)
  testthat::expect_equal(e$estimate, e$complete_data)
})


testthat::test_that("transition probabilities are right for a defective (boundary) rate matrix, and the joint draw runs", {
  a <- 1.3; Q <- matrix(c(-a, a, 0, 0), 2, byrow = TRUE)     # one rate exactly 0: a Jordan block, not diagonalisable
  P <- .edge_P(Q, c(0, 0.4, 2))
  testthat::expect_equal(P[[1]], diag(2))
  testthat::expect_equal(P[[2]], matrix(c(exp(-a * 0.4), 1 - exp(-a * 0.4), 0, 1), 2, byrow = TRUE), tolerance = 1e-12)
  sued <- function(up, down, K = 4) { Q <- matrix(0, K, K); for (i in 1:(K - 1)) { Q[i, i + 1] <- up; Q[i + 1, i] <- down }; diag(Q) <- -rowSums(Q); Q }
  for (Qb in list(sued(2, 0), sued(0, 2))) {
    Pb <- .edge_P(Qb, c(0.1, 0.7))
    testthat::expect_equal(Pb[[1]] %*% Pb[[1]], .edge_P(Qb, 0.2)[[1]], tolerance = 1e-10)   # semigroup
    testthat::expect_equal(rowSums(Pb[[2]]), rep(1, 4))
  }
  Qb <- sued(2, 0); tip <- castor::simulate_mk_model(tree, Qb, root_probabilities = c(.7, .1, .1, .1))$tip_states
  tip[sample(40, 12)] <- NA; K <- 4; prior <- .empirical_prior(tip, K)
  d <- mk_joint_draw(tree, tip, Qb, prior, 4000L)
  miss <- which(is.na(tip))
  emp <- t(vapply(miss, function(i) tabulate(d[, i], K) / nrow(d), numeric(K)))
  exact <- t(vapply(miss, function(i) brute_marginal(tip, Qb, prior, i), numeric(K)))
  testthat::expect_lt(max(abs(emp - exact)), 0.035)
})

testthat::test_that("a castor failure on one trait leaves that trait missing and fills the others", {
  orig <- .castor_fit
  .castor_fit <<- function(tree, tip, K, rate_model) if (rate_model == "SRD") stop("boom") else orig(tree, tip, K, rate_model)
  on.exit(.castor_fit <<- orig, add = TRUE)
  a <- mi_castor(cell, 3L, proper = TRUE)
  testthat::expect_equal(a$diag$error[a$diag$trait == "ord"], "boom")
  testthat::expect_true(all(is.na(a$datasets[[1]]$ord[cell$mask[, "ord"]])))
  testthat::expect_false(anyNA(a$datasets[[1]]$bin)); testthat::expect_false(anyNA(a$datasets[[1]]$cat3))
  s <- score_discrete("x", a$datasets, cell$truth, cell$mask)
  testthat::expect_setequal(s$trait, c("bin", "cat3"))      # all-NA cells dropped: ord has no scorable cell
})

testthat::test_that("score_discrete: fractional ties at the edge of the 95% set, one-class traits skipped, all-NA cells counted", {
  tr <- data.frame(row.names = paste0("s", 1:3), k = factor(c("A", "B", "A"), levels = c("A", "B", "C")),
                   o = factor(c("L2", "L2", "L2"), levels = c("L1", "L2", "L3"), ordered = TRUE))
  mk <- matrix(c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE), 3, dimnames = list(rownames(tr), c("k", "o")))
  draw <- function(cls2, cls3) { d <- tr; d$k <- factor(c("A", cls2, cls3), levels = levels(tr$k)); d }
  sets <- c(rep(list(draw("A", NA)), 18), list(draw("B", NA)), list(draw("C", NA)))   # cell 2: 18 A, 1 B, 1 C; cell 3: all NA
  s <- score_discrete("x", sets, tr, mk)
  testthat::expect_equal(nrow(s), 1L)                        # o has one realised class in the complete data: skipped
  testthat::expect_equal(s$n_cells, 1L); testthat::expect_equal(s$n_cells_all_na, 1L)
  testthat::expect_equal(s$set_size, 2); testthat::expect_equal(s$set_coverage, 0.5)   # B tied with C for the second place
  testthat::expect_equal(s$accuracy, 0)
})

testthat::test_that("c1 ~ bin: undefined when bin is constant in every imputed dataset; per-dataset target is rho x GLS slope of L5", {
  one <- cell$truth; one$bin[] <- levels(one$bin)[2L]
  testthat::expect_null(score_discrete_estimand("x", rep(list(one), 3), cell$truth, cell$tree, NULL, nrow(one)))
  eig <- pagel_eigen(cell$tree, rownames(cell$truth)); L5 <- cell$L[rownames(cell$truth), 5]
  e <- score_discrete_estimand("x", rep(list(cell$truth), 3), cell$truth, cell$tree, eig, nrow(one), L5 = L5,
                               lambda = 0.7, rho = 0.5)
  x01 <- as.numeric(cell$truth$bin == "yes"); V <- 0.7 * cov2cor(ape::vcv(cell$tree))[rownames(cell$truth), rownames(cell$truth)] + 0.3 * diag(nrow(one))
  X <- cbind(1, x01); b <- solve(t(X) %*% solve(V, X), t(X) %*% solve(V, L5))[2]
  testthat::expect_equal(e$target_cond, 0.5 * b, tolerance = 1e-8)
  e0 <- score_discrete_estimand("x", rep(list(cell$truth), 3), cell$truth, cell$tree, eig, nrow(one), L5 = L5, lambda = 0.7, rho = 0)
  testthat::expect_equal(e0$target_cond, 0)
})

testthat::test_that("the joint draw stops (never returns NA) when the observed tips are impossible under Q", {
  tip2 <- rep(c(1L, 2L), 20); tip2[1:5] <- NA
  testthat::expect_error(mk_joint_draw(tree, tip2, matrix(0, 2, 2), c(.5, .5), 3L), "zero likelihood")
})

testthat::test_that("freq A rejects bootstrap rates under which the observed tips are impossible, and never leaves NA", {
  orig <- .castor_fit; calls <- 0L
  .castor_fit <<- function(tree, tip, K, rate_model) {             # every refit after the first fit returns Q = 0
    calls <<- calls + 1L; if (calls == 1L) orig(tree, tip, K, rate_model) else matrix(0, K, K)
  }
  on.exit(.castor_fit <<- orig, add = TRUE)
  c1 <- cell; c1$truth <- cell$truth[, c("c1", "bin")]; c1$df_miss <- cell$df_miss[, c("c1", "bin")]
  a <- mi_castor(c1, 2L, proper = TRUE)
  testthat::expect_equal(a$diag$n_fallback, 2L); testthat::expect_equal(a$diag$n_reject, 40L)
  testthat::expect_false(anyNA(a$datasets[[1]]$bin)); testthat::expect_false(anyNA(a$datasets[[2]]$bin))
})
