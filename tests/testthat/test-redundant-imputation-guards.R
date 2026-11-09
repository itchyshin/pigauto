test_that("baseline-only fits retain explicit safe weights with and without validation", {
  set.seed(71)
  tree <- ape::rtree(30)
  traits <- data.frame(x = sin(seq_len(30) / 3), row.names = tree$tip.label)
  traits$x[c(3, 8)] <- NA_real_
  no_split <- impute(traits, tree, trait_types = c(x = "continuous"),
    log_transform = FALSE, gnn = FALSE, missing_frac = 0,
    verbose = FALSE, seed = 41, predict_method = "per_column")
  expect_equal(unname(no_split$fit$r_cal_bm), 1)
  expect_equal(unname(no_split$fit$r_cal_gnn), 0)
  expect_equal(unname(no_split$fit$r_cal_mean), 0)
  expect_true(all(is.finite(no_split$prediction$imputed_latent)))
  expect_equal(unname(no_split$prediction$imputed_latent),
               unname(no_split$fit$baseline_full$mu), tolerance = 1e-8)
  expect_warning(validated <- impute(traits, tree, trait_types = c(x = "continuous"),
    log_transform = FALSE, gnn = FALSE, missing_frac = .2,
    verbose = FALSE, seed = 41, predict_method = "per_column"), "Small validation set")
  expect_equal(unname(validated$fit$r_cal_gnn), 0)
  expect_equal(unname(validated$fit$r_cal_bm + validated$fit$r_cal_gnn +
                      validated$fit$r_cal_mean), 1)
  expect_true(all(is.finite(validated$prediction$imputed_latent)))
  expect_equal(unname(validated$prediction$imputed_latent),
               unname(validated$fit$baseline_full$mu), tolerance = 1e-8)
})

test_that("public active rankings retain whole-row composition and ZI fallbacks", {
  set.seed(88)
  tree <- ape::rtree(12)
  i <- seq_len(12)
  p <- seq(.15, .8, length.out = 12)
  traits <- data.frame(x = sin(i / 3) + i / 10,
    b = factor(rep(c("no", "yes"), 6)), cat = factor(rep(c("a", "b", "c"), 4)),
    z = as.integer(c(2, 0, 0, 3, rep(0, 8))), p = p, q = 1 - p,
    row.names = tree$tip.label)
  traits[c(10, 11), c("x", "b", "cat", "z")] <- NA
  traits$p[c(10, 11)] <- NA_real_
  traits$q[11] <- NA_real_
  result <- impute(traits, tree, trait_types = c(x = "continuous", z = "zi_count"),
    multi_proportion_groups = list(comp = c("p", "q")), gnn = FALSE,
    missing_frac = 0, verbose = FALSE, seed = 12, predict_method = "per_column")
  X <- result$data$X_scaled
  map <- result$data$trait_map
  comp_cols <- map$comp$latent_cols
  expect_true(!is.na(traits$q[10]))
  expect_true(all(is.na(X[c(10, 11), comp_cols, drop = FALSE])))
  expect_true(all(is.finite(X[-c(10, 11), comp_cols, drop = FALSE])))
  expected_species <- tree$tip.label[c(10, 11)]
  rankings <- lapply(c("continuous", "binary", "categorical", "zi_count", "multi_proportion"),
    function(type) suggest_next_observation(result, types = type, top_n = 100))
  names(rankings) <- c("continuous", "binary", "categorical", "zi_count", "multi_proportion")
  for (type in names(rankings)) {
    ranking <- rankings[[type]]
    expect_s3_class(ranking, "pigauto_active")
    expect_equal(nrow(ranking), 2L, info = type)
    expect_setequal(ranking$species, expected_species)
    expect_true(all(is.finite(ranking$delta)) && all(ranking$delta >= 0), info = type)
  }
  zi <- rankings$zi_count
  mag_col <- map$z$latent_cols[2]
  expect_equal(sum(is.finite(X[, mag_col])), 2L)
  expect_equal(zi$delta_var_total, c(0, 0))
  expect_true(all(is.finite(zi$delta_entropy_total)) && all(zi$delta_entropy_total > 0))
  expect_true(all(is.na(rankings$continuous$delta_entropy_total)))
  expect_true(all(is.na(rankings$binary$delta_var_total)))
  expect_true(all(is.na(rankings$categorical$delta_var_total)))

  # Independently recompute each continuous BM variance reduction by
  # conditioning on each candidate while holding the original rate fixed.
  bm_reference <- function(y) {
    obs <- which(!is.na(y)); miss <- which(is.na(y)); R <- phylo_cor_matrix(tree)
    nugget <- 1e-6
    Roo <- R[obs, obs, drop = FALSE] + diag(nugget, length(obs))
    inv1 <- solve(Roo, rep(1, length(obs)))
    invy <- solve(Roo, y[obs])
    mu <- sum(invy) / sum(inv1)
    e <- y[obs] - mu
    sigma2 <- as.numeric(crossprod(e, solve(Roo, e))) / (length(obs) - 1L)
    total_var <- function(observed, missing) {
      if (!length(missing)) return(0)
      A <- R[observed, observed, drop = FALSE] + diag(nugget, length(observed))
      B <- R[observed, missing, drop = FALSE]
      sum(diag(R[missing, missing, drop = FALSE] - t(B) %*% solve(A, B)))
    }
    before <- total_var(obs, miss)
    vapply(miss, function(k) sigma2 * (before - total_var(c(obs, k), setdiff(miss, k))),
           numeric(1))
  }
  ref_x <- bm_reference(X[, map$x$latent_cols])
  ranked_x <- rankings$continuous
  expect_equal(unname(ranked_x$delta[match(expected_species, ranked_x$species)]),
               unname(ref_x), tolerance = 1e-5)
  ref_comp <- Reduce(`+`, lapply(comp_cols, function(j) bm_reference(X[, j])))
  ranked_comp <- rankings$multi_proportion
  expect_equal(unname(ranked_comp$delta[match(expected_species, ranked_comp$species)]),
               unname(ref_comp), tolerance = 1e-5)
  # Older saved results can lack multi_obs and the cached correlation matrix.
  # Their valid whole-row composition mask is retained while the tree rebuilds R.
  legacy <- result
  legacy$data$multi_obs <- NULL
  legacy$fit$graph$R_phy <- NULL
  legacy_comp <- suggest_next_observation(legacy, types = "multi_proportion", top_n = 100)
  expect_equal(legacy_comp, ranked_comp, tolerance = 1e-8)
})
