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
  for (proper in c(FALSE, TRUE)) {
    a <- mi_castor(cell, 4L, proper = proper)
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
  testthat::expect_true(all(f$sets[[2]]$b == "yes")); testthat::expect_equal(f$n_filled[["b"]], 4L)
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
