# Tests for multi_impute(), with_imputations(), and pool_mi().
# Uses tiny synthetic data (small tree, short training) to keep runtime low.

# ---- Helpers ----------------------------------------------------------------

make_mi_test_data <- function(n = 40, p = 3, miss_frac = 0.15, seed = 123) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label
  df <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5,
    tr3 = abs(stats::rnorm(n)) + 0.5
  )
  # Punch holes in the matrix so we have actual missing cells to impute.
  n_cells   <- n * p
  n_missing <- max(1L, round(n_cells * miss_frac))
  idx <- sample.int(n_cells, n_missing)
  m   <- as.matrix(df)
  m[idx] <- NA
  df <- as.data.frame(m)
  rownames(df) <- sp
  list(tree = tree, df = df, n_missing = n_missing)
}

quick_mi <- function(m = 3L, seed = 123, draws_method = "conformal") {
  td <- make_mi_test_data(seed = seed)
  mi <- multi_impute(
    traits        = td$df,
    tree          = td$tree,
    m             = m,
    draws_method  = draws_method,
    epochs        = 20L,
    missing_frac  = 0.25,
    verbose       = FALSE,
    seed          = seed,
    eval_every    = 10L,
    patience      = 5L
  )
  # Include traits_with_na so tests can identify originally-missing cells
  list(mi = mi, td = td, traits_with_na = td$df)
}


# ---- 1. multi_impute() structural shape -------------------------------------

test_that("multi_impute returns pigauto_mi with M datasets of the right shape", {
  setup <- quick_mi(m = 3L, seed = 11)
  mi <- setup$mi

  expect_s3_class(mi, "pigauto_mi")
  expect_equal(mi$m, 3L)
  expect_true(is.list(mi$datasets))
  expect_equal(length(mi$datasets), 3L)

  # Each dataset should have the same shape as the input.
  for (d in mi$datasets) {
    expect_s3_class(d, "data.frame")
    expect_equal(nrow(d), nrow(setup$td$df))
    expect_equal(ncol(d), ncol(setup$td$df))
    expect_equal(names(d), names(setup$td$df))
    expect_equal(rownames(d), rownames(setup$td$df))
  }

  # Slots expected by downstream helpers.
  expect_true("imputed_mask" %in% names(mi))
  expect_true("fit"          %in% names(mi))
  expect_true("tree"         %in% names(mi))
})


# ---- 2. Dataset consistency across imputations ------------------------------

test_that("all M datasets have matching rownames, colnames, and column classes", {
  setup <- quick_mi(m = 3L, seed = 12)
  mi <- setup$mi

  ref <- mi$datasets[[1]]
  ref_rn <- rownames(ref)
  ref_cn <- colnames(ref)
  ref_cls <- vapply(ref, function(x) class(x)[1], character(1))

  for (d in mi$datasets) {
    expect_equal(rownames(d), ref_rn)
    expect_equal(colnames(d), ref_cn)
    expect_equal(vapply(d, function(x) class(x)[1], character(1)), ref_cls)
  }
})


# ---- 3. Observed cells are identical across imputations ---------------------

test_that("cells that were observed in the input are identical across datasets", {
  setup <- quick_mi(m = 3L, seed = 13)
  mi <- setup$mi
  orig <- setup$td$df

  # For every originally observed cell, all M datasets should return
  # exactly the input value.
  obs_mask <- !is.na(orig)
  for (d in mi$datasets) {
    for (j in seq_len(ncol(orig))) {
      col_obs <- obs_mask[, j]
      expect_equal(
        d[col_obs, j],
        orig[col_obs, j],
        tolerance = 0,  # exact match expected for observed cells
        info = paste("column", colnames(orig)[j])
      )
    }
  }
})


# ---- 4. Imputed cells vary across draws -------------------------------------

test_that("cells that were originally missing vary across at least some datasets", {
  setup <- quick_mi(m = 3L, seed = 14)
  mi <- setup$mi
  miss <- is.na(setup$td$df)

  # Collect the values at originally-missing positions for each draw.
  if (sum(miss) == 0L) skip("no missing cells in test data")

  mat_stack <- do.call(rbind, lapply(mi$datasets, function(d) {
    as.numeric(as.matrix(d)[miss])
  }))
  # Should have M rows, n_missing cols
  expect_equal(nrow(mat_stack), 3L)

  # At least one missing cell should have non-zero variance across M draws.
  # (BM posterior draws guarantee this even when calibrated gates are 0.)
  col_sds <- apply(mat_stack, 2, stats::sd, na.rm = TRUE)
  expect_true(any(col_sds > 1e-10),
              info = "draws should produce varying imputations")
})


# ---- 5. with_imputations() returns a list of length M -----------------------

test_that("with_imputations() returns a pigauto_mi_fits list of length M", {
  setup <- quick_mi(m = 3L, seed = 15)
  mi <- setup$mi

  fits <- with_imputations(mi, function(d) {
    stats::lm(tr1 ~ tr2 + tr3, data = d)
  }, .progress = FALSE)

  expect_s3_class(fits, "pigauto_mi_fits")
  expect_equal(length(fits), 3L)
  expect_true(all(vapply(fits, inherits, logical(1), "lm")))
  expect_equal(attr(fits, "n_failed"), 0L)
})


# ---- 6. with_imputations() captures errors per-imputation -------------------

test_that("with_imputations() captures errors and continues when .on_error = 'continue'", {
  setup <- quick_mi(m = 3L, seed = 16)
  mi <- setup$mi

  counter <- 0L
  make_partial_fail <- function(d) {
    counter <<- counter + 1L
    if (counter == 2L) stop("synthetic failure on imputation 2")
    stats::lm(tr1 ~ tr2, data = d)
  }

  fits <- suppressWarnings(
    with_imputations(mi, make_partial_fail, .progress = FALSE,
                     .on_error = "continue")
  )

  expect_equal(length(fits), 3L)
  expect_equal(attr(fits, "n_failed"), 1L)
  expect_equal(attr(fits, "failed"), 2L)
  expect_s3_class(fits[[2]], "pigauto_mi_error")
  expect_s3_class(fits[[1]], "lm")
  expect_s3_class(fits[[3]], "lm")
})


# ---- 7. pool_mi() on lm fits returns a tidy data.frame ----------------------

test_that("pool_mi() on lm fits returns a data.frame with Rubin's rules columns", {
  setup <- quick_mi(m = 3L, seed = 17)
  mi <- setup$mi

  fits <- with_imputations(mi, function(d) {
    stats::lm(tr1 ~ tr2 + tr3, data = d)
  }, .progress = FALSE)

  pooled <- pool_mi(fits)

  expect_s3_class(pooled, "pigauto_pooled")
  expect_s3_class(pooled, "data.frame")
  expected_cols <- c("term", "estimate", "std.error", "df", "statistic",
                     "p.value", "conf.low", "conf.high", "fmi", "riv")
  expect_equal(names(pooled), expected_cols)

  # Three coefficients: (Intercept), tr2, tr3
  expect_equal(nrow(pooled), 3L)
  expect_true(all(c("(Intercept)", "tr2", "tr3") %in% pooled$term))
  expect_true(all(is.finite(pooled$estimate)))
  expect_true(all(pooled$std.error > 0))
  expect_true(all(pooled$df > 0))
  expect_equal(attr(pooled, "m"), 3L)
})


# ---- 8. Numerical correctness of Rubin's rules against hand calculation ----

test_that("pool_mi() matches a hand-computed Rubin's rules reference (M=3, p=2)", {
  # Three synthetic fits with known coef and vcov.
  make_fake <- function(beta, var_diag) {
    V <- diag(var_diag)
    dimnames(V) <- list(names(beta), names(beta))
    list(beta = beta, V = V)
  }
  fits <- list(
    make_fake(c(a = 1.00, b = 2.00), c(0.010, 0.040)),
    make_fake(c(a = 1.20, b = 2.10), c(0.020, 0.050)),
    make_fake(c(a = 0.90, b = 1.80), c(0.015, 0.030))
  )

  pooled <- pool_mi(
    fits,
    coef_fun = function(f) f$beta,
    vcov_fun = function(f) f$V
  )

  M  <- 3L
  beta_a <- c(1.00, 1.20, 0.90)
  beta_b <- c(2.00, 2.10, 1.80)
  theta_a <- mean(beta_a)
  theta_b <- mean(beta_b)
  W_a <- mean(c(0.010, 0.020, 0.015))
  W_b <- mean(c(0.040, 0.050, 0.030))
  B_a <- stats::var(beta_a)
  B_b <- stats::var(beta_b)
  T_a <- W_a + (1 + 1 / M) * B_a
  T_b <- W_b + (1 + 1 / M) * B_b
  r_a <- (1 + 1 / M) * B_a / W_a
  r_b <- (1 + 1 / M) * B_b / W_b
  df_a <- (M - 1) * (1 + 1 / r_a)^2
  df_b <- (M - 1) * (1 + 1 / r_b)^2

  # Order of rows in pooled should match input coef names.
  row_a <- which(pooled$term == "a")
  row_b <- which(pooled$term == "b")

  expect_equal(pooled$estimate[row_a],  theta_a,    tolerance = 1e-10)
  expect_equal(pooled$estimate[row_b],  theta_b,    tolerance = 1e-10)
  expect_equal(pooled$std.error[row_a], sqrt(T_a),  tolerance = 1e-10)
  expect_equal(pooled$std.error[row_b], sqrt(T_b),  tolerance = 1e-10)
  expect_equal(pooled$riv[row_a],       r_a,        tolerance = 1e-10)
  expect_equal(pooled$riv[row_b],       r_b,        tolerance = 1e-10)
  expect_equal(pooled$df[row_a],        df_a,       tolerance = 1e-10)
  expect_equal(pooled$df[row_b],        df_b,       tolerance = 1e-10)

  # fmi bounded in [0, 1] and positive because B > 0.
  expect_true(all(pooled$fmi > 0 & pooled$fmi < 1))
})


# ---- 9. End-to-end with glmmTMB + propto() ----------------------------------

test_that("pool_mi() works end-to-end with glmmTMB + propto()", {
  skip_if_not_installed("glmmTMB")
  # `propto()` is an internal glmmTMB helper resolved by glmmTMB's formula
  # parser; it is not exported, so we attach glmmTMB to put it on the
  # search path instead of writing glmmTMB::propto.
  # Note: we deliberately do NOT detach/unload glmmTMB here. Unloading its
  # DLL breaks TMB finalizers used by other packages (e.g. Rphylopars) and
  # causes `FreeADFunObject` errors in downstream test files. Leaving it
  # attached for the rest of the test run is harmless.
  suppressPackageStartupMessages(library(glmmTMB))

  # propto() needs more species than predictors to be identifiable; the
  # 40-species default quick_mi() is too small, so build an 80-species
  # tree here.
  td <- make_mi_test_data(n = 80, seed = 18)
  mi <- multi_impute(
    traits       = td$df,
    tree         = td$tree,
    m            = 3L,
    epochs       = 20L,
    missing_frac = 0.25,
    verbose      = FALSE,
    seed         = 18,
    eval_every   = 10L,
    patience     = 5L
  )
  tree <- mi$tree
  # Vphy must be a CORRELATION matrix, not a covariance: propto() in
  # glmmTMB estimates sigma^2 freely, so the raw vcv(tree) (diagonal
  # = tree height) would double-count the variance scale.
  Vphy <- cov2cor(ape::vcv(tree))

  fits <- suppressWarnings(with_imputations(mi, function(d) {
    d$species <- factor(rownames(d), levels = rownames(Vphy))
    d$dummy   <- factor(1)
    glmmTMB(
      tr1 ~ tr2 + propto(0 + species | dummy, Vphy),
      data = d
    )
  }, .progress = FALSE, .on_error = "continue"))

  # Some fits may still fail on synthetic data; require at least 2 successful.
  ok <- vapply(fits, function(f) !inherits(f, "pigauto_mi_error"),
               logical(1))
  skip_if(sum(ok) < 2L,
          "Fewer than 2 glmmTMB fits succeeded on synthetic test data")

  pooled <- suppressWarnings(pool_mi(
    fits[ok],
    coef_fun = function(f) fixef(f)$cond,
    vcov_fun = function(f) stats::vcov(f)$cond
  ))

  expect_s3_class(pooled, "pigauto_pooled")
  expect_true(nrow(pooled) >= 2L)
  expect_true(all(is.finite(pooled$estimate)))
  expect_true(all(pooled$std.error >= 0))
})


# ---- 10. pool_mi() rejects MCMCglmm fits ------------------------------------

test_that("pool_mi() errors cleanly when passed MCMCglmm fits", {
  # Fake MCMCglmm-classed objects -- no need to actually run MCMCglmm.
  fake <- structure(list(), class = "MCMCglmm")
  expect_error(
    pool_mi(list(fake, fake)),
    regexp = "MCMCglmm|posterior"
  )
})


# ---- 11. pool_mi() errors on inconsistent coefficient names -----------------

test_that("pool_mi() errors when fits have inconsistent coefficient names", {
  make_fake <- function(beta, var_diag) {
    V <- diag(var_diag)
    dimnames(V) <- list(names(beta), names(beta))
    list(beta = beta, V = V)
  }
  fits <- list(
    make_fake(c(a = 1.0, b = 2.0), c(0.01, 0.02)),
    make_fake(c(a = 1.1, c = 2.1), c(0.01, 0.02))   # 'c' instead of 'b'
  )
  expect_error(
    pool_mi(
      fits,
      coef_fun = function(f) f$beta,
      vcov_fun = function(f) f$V
    ),
    regexp = "names differ|Rubin"
  )
})


# ---- 12. multi_impute_trees() structural shape --------------------------------

test_that("multi_impute_trees returns pigauto_mi_trees with T*m datasets", {
  # Build 3 small random trees sharing the same tip labels
  n <- 40
  set.seed(200)
  tree1 <- ape::rtree(n)
  sp <- tree1$tip.label

  # Make 2 more trees by randomly perturbing edge lengths
  tree2 <- tree1; tree2$edge.length <- tree1$edge.length * stats::runif(length(tree1$edge.length), 0.5, 1.5)
  tree3 <- tree1; tree3$edge.length <- tree1$edge.length * stats::runif(length(tree1$edge.length), 0.5, 1.5)
  trees <- list(tree1, tree2, tree3)
  class(trees) <- "multiPhylo"

  df <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5,
    tr3 = abs(stats::rnorm(n)) + 0.5
  )
  # Punch some holes
  m <- as.matrix(df)
  idx <- sample.int(n * 3, 15)
  m[idx] <- NA
  df <- as.data.frame(m); rownames(df) <- sp

  # share_gnn = FALSE exercises the legacy per-tree-fit path, which this
  # test asserts on via $fits. suppressWarnings() silences the expected
  # low-M runtime warning triggered by the tiny T*m used for test speed.
  mi_t <- suppressWarnings(multi_impute_trees(
    traits       = df,
    trees        = trees,
    m_per_tree   = 2L,
    epochs       = 20L,
    missing_frac = 0.25,
    verbose      = FALSE,
    seed         = 200L,
    share_gnn    = FALSE,
    eval_every   = 10L,
    patience     = 5L
  ))

  # Class hierarchy

  expect_s3_class(mi_t, "pigauto_mi_trees")
  expect_s3_class(mi_t, "pigauto_mi")

  # Structure
  expect_equal(mi_t$m, 6L)           # 3 trees * 2 imputations
  expect_equal(mi_t$n_trees, 3L)
  expect_equal(mi_t$m_per_tree, 2L)
  expect_equal(length(mi_t$datasets), 6L)
  expect_equal(length(mi_t$tree_index), 6L)
  expect_equal(mi_t$tree_index, c(1L, 1L, 2L, 2L, 3L, 3L))

  # Each dataset has correct shape
  for (d in mi_t$datasets) {
    expect_s3_class(d, "data.frame")
    expect_equal(nrow(d), n)
    expect_equal(ncol(d), 3L)
    expect_equal(names(d), c("tr1", "tr2", "tr3"))
  }

  # Fits list has one per tree (legacy behaviour, share_gnn = FALSE)
  expect_equal(length(mi_t$fits), 3L)
  for (f in mi_t$fits) {
    expect_s3_class(f, "pigauto_fit")
  }
})


# ---- 13. multi_impute_trees() compatible with with_imputations() + pool_mi() --

test_that("multi_impute_trees result works with with_imputations() and pool_mi()", {
  n <- 40
  set.seed(210)
  tree1 <- ape::rtree(n)
  sp <- tree1$tip.label
  tree2 <- tree1; tree2$edge.length <- tree1$edge.length * stats::runif(length(tree1$edge.length), 0.5, 1.5)
  trees <- list(tree1, tree2)

  df <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5
  )
  m <- as.matrix(df); idx <- sample.int(n * 2, 10); m[idx] <- NA
  df <- as.data.frame(m); rownames(df) <- sp

  # suppressWarnings() silences the expected low-M warning (T*m = 4 < 10
  # in this speed-tuned test). Default share_gnn = TRUE is fine here
  # because this test only exercises the downstream with_imputations() +
  # pool_mi() compatibility, not the $fits / $fit distinction.
  mi_t <- suppressWarnings(multi_impute_trees(
    traits = df, trees = trees, m_per_tree = 2L,
    epochs = 20L, missing_frac = 0.25, verbose = FALSE,
    seed = 210L, eval_every = 10L, patience = 5L
  ))

  # with_imputations should accept it via pigauto_mi inheritance
  fits <- with_imputations(mi_t, function(d) {
    stats::lm(tr1 ~ tr2, data = d)
  }, .progress = FALSE)

  expect_s3_class(fits, "pigauto_mi_fits")
  expect_equal(length(fits), 4L)  # 2 trees * 2 imputations

  pooled <- pool_mi(fits)
  expect_s3_class(pooled, "pigauto_pooled")
  expect_equal(nrow(pooled), 2L)  # (Intercept) + tr2
  expect_true(all(is.finite(pooled$estimate)))
  expect_true(all(pooled$std.error > 0))
  expect_equal(attr(pooled, "m"), 4L)
})


# ---- 14. draws_method = "conformal" explicitly ----------------------------

test_that("draws_method='conformal' stores method and produces non-zero between-dataset variance", {
  setup <- quick_mi(m = 5L, seed = 77L, draws_method = "conformal")
  mi    <- setup$mi

  expect_equal(mi$draws_method, "conformal")

  # Collect values at originally-missing positions
  miss <- is.na(setup$traits_with_na[["tr1"]])
  if (sum(miss) == 0L) skip("no missing cells in test data")

  mat_stack <- do.call(rbind, lapply(mi$datasets, function(d) {
    as.numeric(d[["tr1"]][miss])
  }))
  col_sds <- apply(mat_stack, 2, stats::sd, na.rm = TRUE)
  expect_true(any(col_sds > 1e-10),
              info = "conformal draws should produce non-zero between-dataset variance")

  # Observed cells must be identical across all M datasets
  obs_vals <- lapply(mi$datasets, function(d) d[["tr1"]][!miss])
  for (i in seq_along(obs_vals)[-1])
    expect_identical(obs_vals[[i]], obs_vals[[1]],
                     info = paste("dataset", i, "should not alter observed cells"))
})

test_that("multi_impute() multi-obs aligns datasets when input is shuffled", {
  # Regression test for the multi_impute path of the row-alignment bug
  # found 2026-04-26 (commit a3e6d39).  preprocess_traits() reorders rows
  # internally to tree-tip order, and build_completed() needs the inverse
  # mapping to write predictions back to the user's input row order.  This
  # test exercises both draws_method paths (default conformal + mc_dropout)
  # to make sure both code paths in multi_impute.R were patched.
  set.seed(20260426)
  tree <- ape::rtree(15)
  # Two obs per species, species column SHUFFLED.
  sp_shuffled <- sample(rep(tree$tip.label, each = 2L))
  sp_truth <- setNames(rnorm(15, mean = 5, sd = 3), tree$tip.label)
  trait_full <- sp_truth[sp_shuffled] + rnorm(length(sp_shuffled), sd = 0.05)
  df <- data.frame(species = sp_shuffled, trait = trait_full,
                    stringsAsFactors = FALSE)

  # Mask one of each species's two observations.
  mask_idx <- vapply(unique(df$species),
                      function(sp) sample(which(df$species == sp), 1L),
                      integer(1))
  df_obs <- df
  df_obs$trait[mask_idx] <- NA

  for (method in c("conformal", "mc_dropout")) {
    mi <- multi_impute(
      traits        = df_obs,
      tree          = tree,
      species_col   = "species",
      m             = 2L,
      draws_method  = method,
      epochs        = 30L,
      missing_frac  = 0.0,
      verbose       = FALSE,
      seed          = 1L
    )

    for (k in seq_along(mi$datasets)) {
      d <- mi$datasets[[k]]
      observed_idx <- setdiff(seq_len(nrow(df)), mask_idx)
      # Observed rows preserved exactly (sanity check).
      expect_equal(d$trait[observed_idx], df$trait[observed_idx],
                    tolerance = 1e-8,
                    label = sprintf("draws_method=%s, dataset %d, observed", method, k))

      # Imputed rows must align with their own species.  We allow
      # generous tolerance because mc_dropout adds Gaussian noise.
      imp_pred <- d$trait[mask_idx]
      imp_truth <- sp_truth[as.character(df$species[mask_idx])]
      r <- cor(imp_pred, imp_truth)
      expect_gt(r, 0.85,
                 label = sprintf("cor(pred, truth) for draws_method=%s dataset %d (got %.3f)",
                                  method, k, r))
    }
  }
})

test_that("multi_impute(draws_method='conformal') produces non-zero between-draw variance on shuffled input", {
  # Regression test for the row-alignment bug in .sample_conformal_draw
  # found by adversarial review (commit fixing it: pending).
  # The bug: imputed_mask is in user-input row order, but pred$imputed is
  # in internal (tree-tip) order.  .sample_conformal_draw was indexing
  # imp[[nm]][rows] with rows from imputed_mask, so on shuffled-input
  # data it perturbed the WRONG cells, and build_completed re-aligned
  # the unperturbed (deterministic) values back into the user's missing
  # cells — producing identical datasets across all M draws (zero
  # between-imputation variance).
  set.seed(20260428)
  tree <- ape::rtree(15)
  # Single-obs but ROWNAMES are SHUFFLED relative to tree$tip.label so
  # preprocess_traits has to reorder.
  shuffled_sp <- sample(tree$tip.label)
  df <- data.frame(
    row.names = shuffled_sp,
    tr = abs(rnorm(15)) + 0.5
  )
  # Mask 5 of 15 cells.
  miss_idx <- sample(15, 5)
  df_obs <- df
  df_obs$tr[miss_idx] <- NA

  mi <- multi_impute(
    traits        = df_obs,
    tree          = tree,
    m             = 5L,
    draws_method  = "conformal",
    epochs        = 30L,
    missing_frac  = 0.0,
    verbose       = FALSE,
    seed          = 1L
  )

  # Pull values at the originally-missing cells from each of the 5 datasets.
  draws_at_miss <- sapply(mi$datasets, function(d) d$tr[miss_idx])
  # Compute SD across the 5 draws at each masked cell.
  per_cell_sd <- apply(draws_at_miss, 1L, stats::sd)

  # Pre-fix bug behaviour: per_cell_sd ≡ 0 for all cells (zero variance
  # because the wrong cells were perturbed).  Post-fix: the SD should be
  # strictly positive at every masked cell because conformal-width
  # Gaussian sampling was correctly applied.
  expect_true(all(per_cell_sd > 1e-8),
              info = sprintf(
                "All 5 masked cells should have non-zero SD across draws. Got: %s",
                paste(round(per_cell_sd, 6), collapse = ", ")))
})
