# tests/testthat/test-safety-nets.R

.load_avonet <- function() {
  data(avonet300, tree300, envir = environment())
  df <- avonet300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  list(df = df, tree = tree300)
}

test_that("Phase B3 warning fires on imbalanced AVONET300 Migration trait", {
  skip_if_no_libtorch()
  withr::local_seed(1L)
  dat <- .load_avonet()
  
  expect_warning(
    impute(dat$df, dat$tree, n_imputations = 1L, pool_method = "median", epochs = 1L),
    regexp = "Imbalanced K>=3 ordinal trait 'Migration'"
  )
})

test_that("Phase B3 warning stays silent for balanced K=5 ordinal trait", {
  skip_if_no_libtorch()
  withr::local_seed(2L)
  dat <- .load_avonet()
  
  k5_levels <- c("A", "B", "C", "D", "E")
  clean_df <- data.frame(
    Dummy_Size = factor(rep(k5_levels, length.out = nrow(dat$df)), 
                        levels = k5_levels, ordered = TRUE),
    row.names = rownames(dat$df)
  )
  
  clean_df$Dummy_Size[c(1, 2, 3)] <- NA
  
  # CRITICAL FIX: Instead of expect_no_warning (which fails on ANY unrelated warning),
  # we use a custom handler to ensure our specific warning DOES NOT fire.
  warning_fired <- FALSE
  withCallingHandlers({
    impute(clean_df, dat$tree, n_imputations = 1L, pool_method = "median", epochs = 1L)
  }, warning = function(w) {
    if (grepl("Imbalanced K>=3 ordinal trait", w$message)) {
      warning_fired <<- TRUE
    }
  })
  
  expect_false(warning_fired)
})

test_that("Compute scale estimator fires correctly on large N", {
  skip_if_no_libtorch()
  withr::local_seed(3L)
  
  big_tree <- ape::rcoal(5000)
  big_df <- data.frame(
    Dummy_Trait = rnorm(5000),
    row.names = big_tree$tip.label
  )
  big_df$Dummy_Trait[1] <- NA
  
  expect_error(
    withCallingHandlers({
      impute(big_df, big_tree, n_imputations = 10L)
    }, message = function(m) {
      if (grepl("longer", m$message, ignore.case = TRUE)) {
        stop("RIPCORD_PULLED")
      }
    }),
    regexp = "RIPCORD_PULLED"
  )
})
skip_if_no_libtorch()
