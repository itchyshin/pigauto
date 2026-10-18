#!/usr/bin/env Rscript
root <- normalizePath(Sys.getenv("PIGAUTO_ROOT", "."), mustWork = TRUE)
source(file.path(root, "script/validation_20261018/manifest.R"))
source(file.path(root, "script/validation_20261018/evaluate.R"))
must_fail <- function(x) stopifnot(inherits(try(force(x), silent = TRUE), "try-error"))

m <- validation_manifest(root)
stopifnot(isTRUE(validate_manifest(m)))
stopifnot(sum(m$slice == "prediction_core") == 18L,
          sum(m$slice == "prediction_factorial") == 56L,
          sum(m$slice == "prediction_low_signal") == 8L,
          sum(m$slice == "zi_count") == 7L,
          sum(m$slice == "multi_proportion") == 8L,
          sum(m$slice == "posterior") == 8L,
          sum(m$slice == "real") == 4L,
          all(m$reps[m$slice == "posterior"] == 500L),
          all(m$reps[m$slice == "analysis"] == 500L),
          all(m$reps[m$slice == "real"] == 20L))
j <- validation_jobs(m)
stopifnot(nrow(j) == sum(m$reps), nrow(j) == 32240L,
          !anyDuplicated(j$job_id),
          !anyDuplicated(j$seed), all(j$rep >= 1L),
          all(j$rep <= j$reps))
bad <- m; bad$scenario_id[2] <- bad$scenario_id[1]; must_fail(validate_manifest(bad))
bad <- m; bad$reps[1] <- 0L; must_fail(validate_manifest(bad))

small <- j[1:4, ]
small$claim_tier <- "inferential"
make_env <- function(k, status = "ok", truth = 0, estimate = 0,
                     se = 1, lower = -2, upper = 2,
                     complete_reference = 0) {
  x <- small[k, ]
  list(job_id = x$job_id, scenario_id = x$scenario_id,
       rep = x$rep, seed = x$seed, status = status,
       source_identity = "fresh",
       results = data.frame(job_id = x$job_id, method = "pigauto",
                            estimand = "slope", status = status, truth = truth,
                            estimate = estimate, se = se, lower = lower,
                            upper = upper, complete_reference = complete_reference,
                            failure_reason = if (status == "ok") NA_character_ else "fit failed"))
}
e <- lapply(1:4, make_env)
v <- validate_records(e, small, identity = "fresh")
stopifnot(nrow(v$attempts) == 4L, nrow(v$results) == 4L)
must_fail(validate_records(c(e, e[1]), small, identity = "fresh"))
must_fail(validate_records(e[-1], small, identity = "fresh"))
stale <- e; stale[[1]]$source_identity <- "old";
must_fail(validate_records(stale, small, identity = "fresh"))
wrong <- e; wrong[[1]]$seed <- wrong[[1]]$seed + 1L;
must_fail(validate_records(wrong, small, identity = "fresh"))
wrong_row <- e; wrong_row[[1]]$results$scenario_id <- "stale_scenario";
must_fail(validate_records(wrong_row, small, identity = "fresh"))
malformed <- e; malformed[[1]]$results$lower <- 3;
must_fail(validate_records(malformed, small, identity = "fresh"))
malformed <- e; malformed[[1]]$results$complete_reference <- NA_real_;
must_fail(validate_records(malformed, small, identity = "fresh"))
failed <- e; failed[[1]] <- make_env(1, status = "failed", estimate = NA_real_,
                                    se = NA_real_, lower = NA_real_, upper = NA_real_)
vf <- validate_records(failed, small, identity = "fresh")
sf <- summarise_records(vf)
stopifnot(sf$n_attempts[1] == 4L, sf$n_failed[1] == 1L,
          sf$n_unavailable[1] == 0L, sf$n_missing[1] == 0L,
          sf$classification[1] == "indeterminate")
partial <- e
partial[[1]]$status <- "failed"
other <- partial[[1]]$results
other$method <- "other"
other$status <- "failed"
other$estimate <- other$se <- other$lower <- other$upper <- NA_real_
other$failure_reason <- "second arm failed"
partial[[1]]$results <- rbind(partial[[1]]$results, other)
vp <- validate_records(partial, small, identity = "fresh")
sp <- summarise_records(vp)
stopifnot(sp$n_ok[sp$method == "pigauto"] == 4L,
          sp$n_failed[sp$method == "other"] == 1L,
          sp$n_missing[sp$method == "other"] == 3L)
unsupported <- e
unsupported[[1]] <- make_env(1, status = "unavailable", estimate = NA_real_,
                             se = NA_real_, lower = NA_real_, upper = NA_real_)
su <- summarise_records(validate_records(unsupported, small, identity = "fresh"))
stopifnot(su$n_attempts[1] == 4L, su$n_ok[1] == 3L,
          su$n_failed[1] == 0L, su$n_unavailable[1] == 1L,
          su$n_missing[1] == 0L, su$classification[1] == "indeterminate")
supported_scope <- e
unavailable_row <- supported_scope[[1]]$results
unavailable_row$estimand <- "unsupported_interval"
unavailable_row$status <- "unavailable"
unavailable_row$estimate <- unavailable_row$se <- unavailable_row$lower <-
  unavailable_row$upper <- NA_real_
unavailable_row$failure_reason <- "point-only method has no interval"
supported_scope[[1]]$results <- rbind(supported_scope[[1]]$results, unavailable_row)
ss <- summarise_records(validate_records(supported_scope, small, identity = "fresh"))
stopifnot(ss$n_unavailable[ss$estimand == "unsupported_interval"] == 1L,
          ss$n_missing[ss$estimand == "unsupported_interval"] == 3L,
          ss$n_failed[ss$estimand == "unsupported_interval"] == 0L)
invalid_scope <- supported_scope
invalid_scope[[1]]$results$status[2L] <- "failed"
must_fail(validate_records(invalid_scope, small, identity = "fresh"))
truthless <- e; for (i in seq_along(truthless)) truthless[[i]]$results$truth <- NA_real_
descriptive_jobs <- small; descriptive_jobs$claim_tier <- "descriptive"
st <- summarise_records(validate_records(truthless, descriptive_jobs, identity = "fresh"))
stopifnot(st$classification[1] == "descriptive")
all_failed <- lapply(seq_len(4), function(i) {
  x <- make_env(i, status = "failed", estimate = NA_real_, se = NA_real_,
                lower = NA_real_, upper = NA_real_)
  x$results <- x$results[FALSE, ]
  x
})
af <- summarise_records(validate_records(all_failed, small, identity = "fresh"))
stopifnot(af$n_attempts[1] == 4L, af$n_failed[1] == 4L,
          af$n_unavailable[1] == 0L, af$n_missing[1] == 0L,
          af$classification[1] == "indeterminate")

known <- e
for (i in seq_along(known)) {
  known[[i]]$results$estimate <- i + 1
  known[[i]]$results$complete_reference <- 1
  known[[i]]$results$lower <- i + 1 - i
  known[[i]]$results$upper <- i + 1 + i
}
sk <- summarise_records(validate_records(known, small, identity = "fresh"))
widths <- 2 * seq_len(4)
diffs <- seq_len(4)
stopifnot(isTRUE(all.equal(sk$mean_interval_width[1], mean(widths))),
          isTRUE(all.equal(sk$interval_width_mcse[1], stats::sd(widths) / 2)),
          isTRUE(all.equal(sk$interval_width_ci_low[1],
                           mean(widths) - 1.96 * stats::sd(widths) / 2)),
          isTRUE(all.equal(sk$interval_width_ci_high[1],
                           mean(widths) + 1.96 * stats::sd(widths) / 2)),
          isTRUE(all.equal(sk$paired_complete_diff[1], mean(diffs))),
          isTRUE(all.equal(sk$paired_complete_diff_mcse[1], stats::sd(diffs) / 2)),
          isTRUE(all.equal(sk$paired_complete_diff_ci_low[1],
                           mean(diffs) - 1.96 * stats::sd(diffs) / 2)),
          isTRUE(all.equal(sk$paired_complete_diff_ci_high[1],
                           mean(diffs) + 1.96 * stats::sd(diffs) / 2)))

paired <- e
for (i in seq_along(paired)) {
  paired[[i]]$results$estimate <- i
  paired[[i]]$results$lower <- i - 2
  paired[[i]]$results$upper <- i + 2
}
for (i in 1:3) {
  comparison <- paired[[i]]$results
  comparison$method <- "comparison"
  comparison$estimate <- i + c(10, 9, 8)[i]
  comparison$lower <- comparison$estimate - 2
  comparison$upper <- comparison$estimate + 2
  if (i == 3L) {
    comparison$status <- "failed"
    comparison$estimate <- comparison$se <- comparison$lower <- comparison$upper <- NA_real_
    comparison$failure_reason <- "comparison fit failed"
    paired[[i]]$status <- "failed"
  }
  paired[[i]]$results <- rbind(paired[[i]]$results, comparison)
}
pc <- paired_method_contrasts(validate_records(paired, small, identity = "fresh"))
stopifnot(nrow(pc) == 1L, pc$method_a == "comparison", pc$method_b == "pigauto",
          pc$n_attempts == 4L, pc$n_a_ok == 2L, pc$n_b_ok == 4L,
          pc$n_pairs == 2L, pc$n_missing_pairs == 2L,
          isTRUE(all.equal(pc$mean_difference, 9.5)),
          isTRUE(all.equal(pc$difference_mcse, 0.5)),
          isTRUE(all.equal(pc$difference_ci_low, 9.5 - 1.96 * 0.5)),
          isTRUE(all.equal(pc$difference_ci_high, 9.5 + 1.96 * 0.5)))

big <- j[j$slice == "posterior" & j$scenario_id == "posterior_001", ]
stopifnot(nrow(big) == 500L)
z <- stats::qnorm((seq_len(500) - 0.5) / 500)
sd_z <- stats::sd(z)
make_big <- function(i, width = 1.96) {
  b <- big[i, ]
  list(job_id = b$job_id, scenario_id = b$scenario_id, rep = b$rep,
       seed = b$seed, status = "ok", source_identity = "fresh",
       results = data.frame(job_id = b$job_id, method = "pigauto",
                            estimand = "slope", status = "ok", truth = 0,
                            estimate = z[i], se = sd_z,
                            lower = z[i] - width, upper = z[i] + width,
                            complete_reference = z[i], failure_reason = NA_character_))
}
set.seed(42)
before_rng <- .Random.seed
good <- summarise_records(validate_records(lapply(seq_len(500), make_big),
                                            big, identity = "fresh"))
stopifnot(identical(before_rng, .Random.seed))
stopifnot(good$classification[1] == "supported", good$n_attempts[1] == 500L,
          good$coverage[1] >= 0.925, good$coverage[1] <= 0.975)
wide <- summarise_records(validate_records(lapply(seq_len(500), make_big, width = 100),
                                            big, identity = "fresh"))
stopifnot(wide$classification[1] == "failed", wide$coverage[1] == 1)
cat("HARNESS_PASS\n")
