# P1-9: observed zi_count ZEROS were excluded from val/test splits.
#
# zi_count uses two latent columns: gate (0/1) and log1p-z magnitude, where
# the magnitude is NA when the observed count is zero. make_missing_splits()
# required ALL latent columns non-NA to call a cell observed, so every
# observed zero was treated as unobserved and could never be drawn into
# val/test — silently restricting zi_count calibration and evaluation to the
# non-zero subset.

make_zi_fixture <- function(n = 60L, seed = 7L) {
  set.seed(seed)
  tr <- ape::rtree(n)
  sp <- tr$tip.label
  # ~40% structural zeros, the rest small counts
  cnt <- ifelse(runif(n) < 0.4, 0L, rpois(n, 4) + 1L)
  df <- data.frame(row.names = sp, z = as.integer(cnt),
                   x = as.numeric(ape::rTraitCont(tr, model = "BM")))
  list(tree = tr, df = df, zeros = sum(cnt == 0L))
}

test_that("[P1-9] zi_count zero cells are eligible for val/test splits", {
  fx <- make_zi_fixture()
  d <- preprocess_traits(fx$df, fx$tree, trait_types = c(z = "zi_count"))
  tmz <- Filter(function(t) identical(t$type, "zi_count"), d$trait_map)[[1]]
  gate_lc <- tmz$latent_cols[1L]; mag_lc <- tmz$latent_cols[2L]

  # Sanity: the fixture really does encode zeros as gate-observed,
  # magnitude-NA — i.e. the situation the bug turned on.
  zero_rows <- which(!is.na(d$X_scaled[, gate_lc]) & is.na(d$X_scaled[, mag_lc]))
  expect_gt(length(zero_rows), 0L)

  sp <- make_missing_splits(d$X_scaled, missing_frac = 0.5, val_frac = 0.5,
                            seed = 3L, trait_map = d$trait_map)
  # Which (row, trait) cells were held out for the zi_count trait?
  zi_j <- which(vapply(d$trait_map, function(t) identical(t$type, "zi_count"),
                       logical(1)))
  n_rows <- nrow(d$X_scaled)
  held_trait <- c(sp$val_idx_trait, sp$test_idx_trait)
  held_rows_zi <- ((held_trait - 1L) %% n_rows) + 1L
  held_j <- ((held_trait - 1L) %/% n_rows) + 1L
  held_rows_zi <- held_rows_zi[held_j == zi_j]

  # The defect: this intersection was ALWAYS empty.
  expect_gt(length(intersect(held_rows_zi, zero_rows)), 0L)
})

test_that("[P1-9] non-zi traits keep the all-columns-observed rule", {
  set.seed(9L)
  n <- 50L; tr <- ape::rtree(n); sp <- tr$tip.label
  df <- data.frame(row.names = sp,
                   a = as.numeric(ape::rTraitCont(tr, model = "BM")),
                   b = factor(sample(c("p", "q", "r"), n, TRUE)))
  df$a[1:5] <- NA_real_
  d <- preprocess_traits(df, tr)
  sp_out <- make_missing_splits(d$X_scaled, missing_frac = 0.4, val_frac = 0.5,
                                seed = 2L, trait_map = d$trait_map)
  # A genuinely missing continuous cell must never be selected as held-out.
  n_rows <- nrow(d$X_scaled)
  held <- c(sp_out$val_idx_trait, sp_out$test_idx_trait)
  held_rows <- ((held - 1L) %% n_rows) + 1L
  held_j <- ((held - 1L) %/% n_rows) + 1L
  a_j <- which(vapply(d$trait_map, function(t) identical(t$name, "a"), logical(1)))
  expect_length(intersect(held_rows[held_j == a_j], 1:5), 0L)
})

test_that("[P1-9] a zi_count fit still runs end to end with zeros held out", {
  skip_if_no_libtorch()
  fx <- make_zi_fixture(n = 60L)
  res <- impute(fx$df, fx$tree, trait_types = c(z = "zi_count"),
                epochs = 20L, verbose = FALSE, seed = 4L)
  expect_s3_class(res, "pigauto_result")
  expect_true(all(!is.na(res$completed$z)))
  expect_true(all(res$completed$z >= 0))
})
