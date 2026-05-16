# Unit tests for script/gha/ helpers.
# These exercise the helpers via testthat so they're covered by the
# regular `devtools::test()` run, even though they live under script/.

.gha_path <- function(...) {
  # tests/testthat run from tests/testthat/; the helpers are at
  # ../../script/gha/<name>.R relative to the test file.
  file.path("..", "..", "script", "gha", ...)
}

test_that("[gha-config] PIGAUTO_CI_CONFIG loads with documented defaults", {
  e <- new.env()
  source(.gha_path("_ci_config.R"), local = e)
  cfg <- e$PIGAUTO_CI_CONFIG
  expect_type(cfg, "list")
  expect_equal(cfg$subset_n,          2000L)
  expect_equal(cfg$n_imputations,     10L)
  expect_equal(cfg$missing_frac,      0.30)
  expect_equal(cfg$seed,              2026L)
  expect_equal(cfg$pool_method,       "median")
  expect_false(cfg$clamp_outliers)
  expect_true(cfg$phylo_signal_gate)
  expect_equal(cfg$mc_cores,          4L)
})

test_that("[gha-config] .normalize_eval() projects to the canonical schema", {
  e <- new.env()
  source(.gha_path("_ci_config.R"), local = e)
  norm <- e$.normalize_eval

  raw <- data.frame(
    trait          = c("mass", "diet"),
    type           = c("continuous", "categorical"),
    imputation_idx = c(1L, 1L),
    rmse           = c(0.41, NA_real_),
    mae            = c(0.30, NA_real_),
    pearson_r      = c(0.85, NA_real_),
    accuracy       = c(NA_real_, 0.72),
    brier          = c(NA_real_, 0.18),
    time_sec       = c(120, 120),
    stringsAsFactors = FALSE
  )

  out <- norm(raw, dataset = "avonet", method = "pigauto_ci")
  expected_cols <- c("dataset", "trait", "type", "method",
                     "imputation_idx", "rmse", "mae", "pearson_r",
                     "accuracy", "brier", "time_sec")
  expect_equal(colnames(out), expected_cols)
  expect_equal(out$dataset, c("avonet", "avonet"))
  expect_equal(out$method,  c("pigauto_ci", "pigauto_ci"))
  expect_equal(nrow(out), 2L)
})

test_that("[gha-mask] latent-cell -> user-col mask projection is correct", {
  skip_if_not_installed("ape")
  e <- new.env()
  source(.gha_path("_ci_config.R"), local = e)

  # Tiny dataset: 8 species, mix of continuous + factor (categorical)
  set.seed(2026L)
  n <- 8L
  df <- data.frame(
    mass = stats::rnorm(n),
    diet = factor(sample(c("a","b","c"), n, replace = TRUE)),
    row.names = paste0("sp", seq_len(n))
  )
  tree <- ape::rcoal(n, tip.label = rownames(df))
  pd <- preprocess_traits(df, tree, log_transform = FALSE)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, val_frac = 0.5,
                                 seed = 1L, trait_map = pd$trait_map)
  # Reproduce the wrapper's latent -> user mask logic
  mask_latent <- matrix(FALSE, nrow = nrow(pd$X_scaled), ncol = ncol(pd$X_scaled))
  mask_latent[splits$test_idx] <- TRUE
  user_mask <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                       dimnames = list(rownames(df), colnames(df)))
  for (k in seq_along(pd$trait_map)) {
    tm <- pd$trait_map[[k]]
    user_cols   <- if (!is.null(tm$input_cols)) tm$input_cols else tm$name
    latent_cols <- tm$latent_cols
    hit_rows <- apply(mask_latent[, latent_cols, drop = FALSE], 1L, any)
    for (uc in user_cols) user_mask[, uc] <- user_mask[, uc] | hit_rows
  }
  # categorical 'diet' expands to multiple latent cols; the projection must
  # collapse them into a single user col
  expect_true(all(colnames(user_mask) == c("mass", "diet")))
  expect_equal(dim(user_mask), c(n, 2L))
  # at least some test cells should have been mapped through
  expect_true(sum(user_mask) > 0L)
})

test_that("[gha-eval] .eval_per_imputation() emits canonical long-format rows", {
  e <- new.env()
  source(.gha_path("_ci_config.R"), local = e)
  eval_fn <- e$.eval_per_imputation

  set.seed(7L)
  truth <- data.frame(
    mass = c(1.0, 2.0, 3.0, 4.0),
    diet = factor(c("a", "b", "a", "b")),
    stringsAsFactors = FALSE
  )
  mask <- matrix(c(FALSE, TRUE,  TRUE,  FALSE,
                   TRUE,  FALSE, FALSE, TRUE),
                 nrow = 4L, ncol = 2L,
                 dimnames = list(NULL, c("mass", "diet")))
  imp1 <- data.frame(
    mass = c(1.0, 2.1, 2.8, 4.0),
    diet = factor(c("a", "b", "a", "b")),
    stringsAsFactors = FALSE
  )
  imp2 <- data.frame(
    mass = c(1.0, 2.4, 3.3, 4.0),
    diet = factor(c("a", "b", "a", "a")),
    stringsAsFactors = FALSE
  )
  trait_map <- list(
    list(name = "mass", type = "continuous", input_cols = "mass"),
    list(name = "diet", type = "categorical", input_cols = "diet")
  )

  out <- eval_fn(list(imp1, imp2), truth, mask, trait_map, t_fit_sec = 99)
  expect_equal(sort(unique(out$trait)), c("diet", "mass"))
  expect_equal(sort(unique(out$imputation_idx)), c(1L, 2L))
  expect_equal(nrow(out), 4L)  # 2 traits × 2 imputations
  # mass imp1: rows 2,3 → t=(2,3), p=(2.1,2.8) → rmse=sqrt((0.01+0.04)/2)≈0.158
  mass_imp1 <- out[out$trait == "mass" & out$imputation_idx == 1L, ]
  expect_equal(mass_imp1$rmse, sqrt((0.01 + 0.04) / 2), tolerance = 1e-6)
  expect_true(is.na(mass_imp1$accuracy))
  # diet imp1: rows 1,4 → t=(a,b), p=(a,b) → acc=1
  diet_imp1 <- out[out$trait == "diet" & out$imputation_idx == 1L, ]
  expect_equal(diet_imp1$accuracy, 1.0)
  # diet imp2: rows 1,4 → t=(a,b), p=(a,a) → acc=0.5
  diet_imp2 <- out[out$trait == "diet" & out$imputation_idx == 2L, ]
  expect_equal(diet_imp2$accuracy, 0.5)
  expect_true(all(out$time_sec == 99))
})
