# tests/testthat/test-safety-nets.R

test_that("Phase B3 warning fires on imbalanced K>=3 ordinal", {
  data(avonet300, tree300)
  df <- avonet300
  rownames(df) <- df$Species_Key; df$Species_Key <- NULL
  df$Migration <- factor(df$Migration, ordered = TRUE)
  
  set.seed(1L)
  hide <- sample(which(!is.na(df$Migration)), 30L)
  df$Migration[hide] <- NA
  
  # We expect a warning matching our regex (epochs=1 so the test runs instantly)
  expect_warning(
    impute(df, tree300, n_imputations = 1L, pool_method = "median", epochs = 1L),
    regexp = "Imbalanced K>=3 ordinal trait"
  )
})

test_that("Phase B3 warning stays silent for balanced K=5 ordinal", {
  data(avonet300, tree300)
  df <- avonet300
  rownames(df) <- df$Species_Key; df$Species_Key <- NULL
  
  # Create a perfectly balanced K=5 ordinal factor
  k5_levels <- c("Tiny", "Small", "Medium", "Large", "Giant")
  df$Dummy_Size <- factor(rep(k5_levels, each = 60), levels = k5_levels, ordered = TRUE)
  
  set.seed(2L)
  hide <- sample(1:300, 30L)
  df$Dummy_Size[hide] <- NA
  
  # We expect NO warning
  expect_no_warning(
    impute(df, tree300, n_imputations = 1L, pool_method = "median", epochs = 1L)
  )
})

test_that("Compute scale estimator message fires on massive dataset", {
  data(avonet300, tree300)
  
  # Create a dummy dataframe with 5000 rows
  df_massive <- data.frame(Mass = rnorm(5000))
  rownames(df_massive) <- paste0("FakeBird_", 1:5000)
  
  # We expect the message to fire (even if the function errors out after due to tree mismatch)
  expect_message(
    tryCatch(impute(df_massive, tree300, n_imputations = 10L, epochs = 1L), error = function(e) NULL),
    regexp = "longer than a default run"
  )
})