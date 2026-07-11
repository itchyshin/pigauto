#!/usr/bin/env Rscript

source("script/mi-validation-v010/0_prepare.R")
source("script/mi-validation-v010/comparator_helpers.R")

config <- mi_validation_config(c(
  "--profile=smoke", "--smcfcs-numit=33", "--jomo-nburn=444",
  "--smcfcs-rjlimit=22222", "--jomo-nbetween=555"
))
default_config <- mi_validation_config("--profile=smoke")
manifest <- make_manifest(config)
stopifnot(
  default_config$jomo_nbetween == 100L,
  default_config$smcfcs_rjlimit == 10000L,
  all(manifest$smcfcs_numit == 33L),
  all(manifest$smcfcs_rjlimit == 22222L),
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

lm_dgp <- simulate_validation_dgp("lm", "auxiliary", seed = 20260711L)
lm_a <- .draw_bayes_norm_lm_imputations(lm_dgp, m = 50L, seed = 42L)
lm_b <- .draw_bayes_norm_lm_imputations(lm_dgp, m = 50L, seed = 42L)
lm_c <- .draw_bayes_norm_lm_imputations(lm_dgp, m = 50L, seed = 43L)
lm_missing <- !is.finite(lm_dgp$observed$x)
lm_observed <- !lm_missing
stopifnot(
  identical(lm_a$datasets, lm_b$datasets),
  !identical(lm_a$datasets, lm_c$datasets),
  identical(lm_a$diagnostics$engine, "bayes_norm_lm"),
  identical(lm_a$diagnostics$imputation_formula, "x ~ y + z + I(z^2)"),
  length(lm_a$diagnostics$warnings) == 0L,
  all(vapply(lm_a$datasets, function(x) all(is.finite(x$x)), logical(1))),
  all(vapply(lm_a$datasets, function(x) {
    identical(x$x[lm_observed], lm_dgp$observed$x[lm_observed])
  }, logical(1))),
  stats::sd(vapply(lm_a$datasets, function(x) {
    x$x[which(lm_missing)[1]]
  }, numeric(1))) > 0
)

wrong_dispatch <- tryCatch(
  .draw_smcfcs_imputations(lm_dgp, m = 2L, seed = 42L),
  error = conditionMessage
)
stopifnot(grepl("only for the logit glm", wrong_dispatch, fixed = TRUE))

cat("comparator helper tests: PASS\n")
