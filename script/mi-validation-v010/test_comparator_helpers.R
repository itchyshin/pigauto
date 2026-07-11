#!/usr/bin/env Rscript

source("script/mi-validation-v010/0_prepare.R")
source("script/mi-validation-v010/comparator_helpers.R")

config <- mi_validation_config(c(
  "--profile=smoke", "--smcfcs-numit=33", "--jomo-nburn=444",
  "--jomo-nbetween=555"
))
default_config <- mi_validation_config("--profile=smoke")
manifest <- make_manifest(config)
stopifnot(
  default_config$jomo_nbetween == 100L,
  all(manifest$smcfcs_numit == 33L),
  all(manifest$jomo_nburn == 444L),
  all(manifest$jomo_nbetween == 555L)
)

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
  .validate_completed_datasets(draws_a, dgp, 5L, "oracle-test")
  wrong_count <- tryCatch(
    .validate_completed_datasets(draws_a[-1], dgp, 5L, "oracle-test"),
    error = conditionMessage
  )
  changed <- draws_a
  changed[[1]]$x[which(observed)[1]] <- changed[[1]]$x[which(observed)[1]] + 1
  changed_observed <- tryCatch(
    .validate_completed_datasets(changed, dgp, 5L, "oracle-test"),
    error = conditionMessage
  )
  stopifnot(
    grepl("returned 4 datasets", wrong_count, fixed = TRUE),
    grepl("changed observed x", changed_observed, fixed = TRUE)
  )
}

cat("comparator helper tests: PASS\n")
