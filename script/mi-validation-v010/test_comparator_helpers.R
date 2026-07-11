#!/usr/bin/env Rscript

source("script/mi-validation-v010/0_prepare.R")
source("script/mi-validation-v010/comparator_helpers.R")

for (dgp_name in c("lm", "glm", "lmer")) {
  dgp <- simulate_validation_dgp(dgp_name, "auxiliary", seed = 20260711L)
  draws_a <- .draw_oracle_conditional_imputations(dgp, m = 5L, seed = 42L)
  draws_b <- .draw_oracle_conditional_imputations(dgp, m = 5L, seed = 42L)
  draws_c <- .draw_oracle_conditional_imputations(dgp, m = 5L, seed = 43L)
  missing <- !is.finite(dgp$observed$x)
  observed <- !missing
  stopifnot(
    identical(draws_a, draws_b),
    !identical(draws_a, draws_c),
    all(vapply(draws_a, function(x) all(is.finite(x$x)), logical(1))),
    all(vapply(draws_a, function(x) {
      identical(x$x[observed], dgp$observed$x[observed])
    }, logical(1))),
    stats::sd(vapply(draws_a, function(x) x$x[which(missing)[1]], numeric(1))) > 0
  )
}

cat("comparator helper tests: PASS\n")
