make_evaluation_fixture <- function() {
  set.seed(42)
  tree <- ape::rtree(60)
  traits <- data.frame(row.names = tree$tip.label,
    mass = exp(rnorm(60, 3, .5)),
    clutch = as.integer(rpois(60, 3) + 1L),
    migr = factor(sample(c('no', 'yes'), 60, TRUE)),
    diet = factor(sample(c('herb', 'carn', 'omni'), 60, TRUE)),
    threat = ordered(sample(c('LC', 'VU', 'EN'), 60, TRUE),
                     levels = c('LC', 'VU', 'EN')))
  pd <- suppressMessages(preprocess_traits(traits, tree))
  spl <- make_missing_splits(pd$X_scaled, missing_frac = .2,
                             seed = 1, trait_map = pd$trait_map)
  fit <- suppressWarnings(fit_pigauto(pd, tree, splits = spl,
                                     gnn = FALSE, verbose = FALSE, seed = 1))
  list(pd = pd, spl = spl, fit = fit)
}

test_that('production predictions cannot silently score their observed cells', {
  z <- make_evaluation_fixture()
  pred <- predict(z$fit)
  expect_error(evaluate_imputation(pred, z$pd$X_scaled, z$spl),
               'evaluation')
})

test_that('public evaluation mode agrees with masked prediction', {
  z <- make_evaluation_fixture()
  held <- c(z$spl$val_idx, z$spl$test_idx)
  pred <- predict(z$fit, evaluation = TRUE)
  masked <- predict(z$fit, .mask_observed_idx = held)
  ev <- suppressWarnings(evaluate_imputation(pred, z$pd$X_scaled, z$spl))
  internal <- suppressWarnings(evaluate(z$fit))
  for (trait in c('mass', 'clutch', 'threat')) {
    actual <- ev$rmse[ev$split == 'test' & ev$trait == trait]
    expected <- internal$value[internal$method == 'pigauto' &
                               internal$trait == trait & internal$metric == 'rmse']
    expect_equal(actual, expected)
  }
  for (trait in c('migr', 'diet')) {
    actual <- ev$accuracy[ev$split == 'test' & ev$trait == trait]
    expected <- internal$value[internal$method == 'pigauto' &
                               internal$trait == trait & internal$metric == 'accuracy']
    expect_equal(actual, expected)
  }
  expect_equal(pred$imputed_latent, masked$imputed_latent)
  expect_gt(max(abs(pred$imputed_latent[held] - z$pd$X_scaled[held])), .1)
  expect_no_error(suppressWarnings(evaluate_imputation(pred, z$pd$X_scaled, z$spl)))
  partial <- predict(z$fit, .mask_observed_idx = z$spl$test_idx)
  expect_error(evaluate_imputation(partial, z$pd$X_scaled, z$spl), 'evaluation')
  old <- pred
  old$evaluation_safe_idx <- NULL
  expect_no_error(suppressWarnings(evaluate_imputation(old, z$pd$X_scaled, z$spl)))
  no_splits <- z$fit
  no_splits$splits <- NULL
  expect_error(predict(no_splits, evaluation = TRUE), 'splits')
  expect_error(predict(z$fit, evaluation = NA), 'evaluation')
  expect_error(predict(z$fit, newdata = z$pd, evaluation = TRUE), 'newdata')
})
