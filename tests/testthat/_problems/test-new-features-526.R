# Extracted from test-new-features.R:526

# prequel ----------------------------------------------------------------------
make_test_data_new <- function(n = 40, p = 2, seed = 100) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label
  df   <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5
  )
  list(tree = tree, df = df)
}
build_quick_fit <- function(seed = 100) {
  td  <- make_test_data_new(n = 40, seed = seed)
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = seed, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = seed)
  list(fit = fit, pd = pd, spl = spl, tree = td$tree)
}

# test -------------------------------------------------------------------------
set.seed(602)
tree <- ape::rtree(10)
df <- data.frame(
    species = rep(tree$tip.label, each = 4),
    y = abs(rnorm(40)) + 0.5
  )
covs <- data.frame(
    temp = rep(c(10, 20, 30, 40), times = 10)  # systematic within-species variation
  )
result <- impute(df, tree, species_col = "species",
                   covariates = covs,
                   epochs = 30L, verbose = FALSE, seed = 602L,
                   eval_every = 10L, patience = 5L,
                   missing_frac = 0.25)
pred <- result$prediction$imputed
for (sp in tree$tip.label[1:3]) {
    sp_rows <- which(df$species == sp)
    sp_preds <- pred$y[sp_rows]
    # All predictions within species should be identical (within numeric tolerance)
    expect_equal(max(sp_preds) - min(sp_preds), 0,
                 tolerance = 1e-6,
                 label = paste("Within-species prediction range for", sp))
  }
