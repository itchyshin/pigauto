# Issue #209: result helpers accept pigauto_result; raw saveRDS names save_pigauto.

tiny_result <- function(gnn = FALSE) {
  set.seed(209L)
  tree <- ape::rtree(12)
  df <- data.frame(
    mass = abs(as.numeric(ape::rTraitCont(tree))) + 0.5,
    row.names = tree$tip.label
  )
  df$mass[1:2] <- NA_real_
  impute(df, tree, epochs = 8L, verbose = FALSE, seed = 209L, gnn = gnn)
}

test_that("plot, plot_uncertainty, and save_pigauto accept a pigauto_result", {
  res <- tiny_result(gnn = FALSE)
  expect_s3_class(res, "pigauto_result")
  pdf(tempfile(fileext = ".pdf"))
  expect_silent(plot(res, type = "gates"))
  dev.off()
  expect_s3_class(plot_uncertainty(res, trait_name = "mass"), "ggplot")
  tmp <- tempfile(fileext = ".pigauto")
  on.exit(unlink(tmp), add = TRUE)
  save_pigauto(res, tmp)
  expect_true(file.exists(tmp))
})

test_that("load_pigauto() then predict() works", {
  res <- tiny_result(gnn = FALSE)
  tmp <- tempfile(fileext = ".pigauto")
  on.exit(unlink(tmp), add = TRUE)
  save_pigauto(res, tmp)
  pred <- predict(load_pigauto(tmp), return_se = FALSE)
  expect_s3_class(pred, "pigauto_pred")
  expect_equal(nrow(pred$imputed), nrow(res$completed))
})

test_that("raw saveRDS/readRDS then predict names save_pigauto", {
  skip_if_no_libtorch()
  res <- tiny_result(gnn = TRUE)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(res$fit, tmp)
  fit2 <- readRDS(tmp)
  # Same-session tensors can still look alive; a raw saveRDS object is not
  # a portable fit. Break the first weight so predict must take the named path.
  if (length(fit2$model_state)) {
    fit2$model_state[[1L]] <- "broken-after-saveRDS"
  }
  expect_error(predict(fit2), "save_pigauto")
})
