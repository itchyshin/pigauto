test_that("finite covariates that overflow scaling get a useful error", {
  tree <- ape::read.tree(text = paste0(
    "(((a:1,b:1):1,(c:1,d:1):1):1,",
    "((e:1,f:1):1,(g:1,h:1):1):1);"
  ))
  traits <- data.frame(x = as.numeric(seq_len(8)), row.names = tree$tip.label)
  extreme_mean <- data.frame(huge = c(1.7e308, rep(-1.7e308, 7)),
                             row.names = tree$tip.label)
  extreme_sd <- data.frame(huge = c(1e200, -1e200, rep(0, 6)),
                           row.names = tree$tip.label)
  expect_true(all(is.finite(extreme_mean$huge)))
  expect_true(all(is.finite(extreme_sd$huge)))
  expect_error(preprocess_traits(traits, tree, covariates = extreme_mean),
               "(?i)huge.*rescal")
  expect_error(preprocess_traits(traits, tree, covariates = extreme_sd),
               "(?i)huge.*rescal")
})

test_that("ordinary and constant numeric covariate encodings are unchanged", {
  tree <- ape::read.tree(text = paste0(
    "(((a:1,b:1):1,(c:1,d:1):1):1,",
    "((e:1,f:1):1,(g:1,h:1):1):1);"
  ))
  traits <- data.frame(x = as.numeric(seq_len(8)), row.names = tree$tip.label)
  cov <- data.frame(z = seq_len(8), constant = rep(5, 8),
                    group = factor(rep(c("a", "b"), 4)),
                    row.names = tree$tip.label)
  result <- preprocess_traits(traits, tree, covariates = cov)
  expect_s3_class(result, "pigauto_data")
  expect_equal(unname(result$covariates[, "z"]), (seq_len(8) - 4.5) / sqrt(6))
  expect_equal(unname(result$covariates[, "constant"]), rep(0, 8))
  expect_equal(unname(result$cov_means[c("z", "constant")]), c(4.5, 5))
  expect_equal(unname(result$cov_sds[c("z", "constant")]), c(sqrt(6), 1))
  expect_equal(unname(result$covariates[, "z"] * result$cov_sds[["z"]] +
                        result$cov_means[["z"]]), as.numeric(seq_len(8)))
  expect_equal(unname(result$covariates[, "group.a"]), rep(c(1, 0), 4))
  expect_equal(unname(result$covariates[, "group.b"]), rep(c(0, 1), 4))
  expect_true(all(is.finite(result$covariates)))
})

test_that("finite tiny covariate spread retains the existing unit SD fallback", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  traits <- data.frame(x = c(1, 2, 3), row.names = tree$tip.label)
  cov <- data.frame(tiny = c(-1e-13, 0, 1e-13), row.names = tree$tip.label)
  result <- preprocess_traits(traits, tree, covariates = cov)
  expect_identical(unname(result$cov_means[["tiny"]]), 0)
  expect_identical(unname(result$cov_sds[["tiny"]]), 1)
  expect_equal(unname(result$covariates[, "tiny"]), cov$tiny, tolerance = 1e-15)
})
