test_that("lazy torch factories reuse definitions but create separate modules", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed())
  for (name in c("GraphTransformerBlock", "ResidualPhyloDAE")) {
    getter <- get(paste0(name, "_module"), envir = asNamespace("pigauto"))
    expect_identical(getter(), getter())
    expect_s3_class(attr(getter(), "module"), "R6ClassGenerator")
  }
  block <- get("GraphTransformerBlock", envir = asNamespace("pigauto"))
  expect_false(identical(block(hidden_dim = 4L), block(hidden_dim = 4L)))
})
