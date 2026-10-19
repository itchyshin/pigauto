# Evidence accounting is independent of fit code. Each planned job must have
# exactly one fresh envelope, including failures and unavailable estimates.
validate_records <- function(envelopes, jobs, identity = NULL) {
  if (!is.list(envelopes) || !is.data.frame(jobs) ||
      !all(c("job_id", "scenario_id", "rep", "seed", "claim_tier") %in% names(jobs)) ||
      anyDuplicated(jobs$job_id)) stop("invalid planned jobs")
  if (length(envelopes) != nrow(jobs)) stop("missing or surplus job envelopes")
  required <- c("job_id", "scenario_id", "rep", "seed", "status", "results")
  ids <- vapply(envelopes, function(x) {
    if (!is.list(x) || !all(required %in% names(x)) ||
        length(x$job_id) != 1L || is.na(x$job_id)) stop("malformed envelope")
    as.character(x$job_id)
  }, character(1))
  if (anyDuplicated(ids) || !setequal(ids, jobs$job_id))
    stop("duplicate, unknown, or missing job identity")
  envelopes <- envelopes[match(jobs$job_id, ids)]
  valid_status <- c("ok", "failed", "unavailable")
  row_required <- c("job_id", "method", "estimand", "status", "truth",
                    "estimate", "se", "lower", "upper", "complete_reference",
                    "failure_reason")
  result_rows <- vector("list", nrow(jobs))
  attempt_rows <- vector("list", nrow(jobs))
  for (i in seq_len(nrow(jobs))) {
    e <- envelopes[[i]]
    j <- jobs[i, , drop = FALSE]
    if (!identical(as.character(e$scenario_id), as.character(j$scenario_id)) ||
        !identical(as.integer(e$rep), as.integer(j$rep)) ||
        !identical(as.integer(e$seed), as.integer(j$seed)) ||
        length(e$status) != 1L || !e$status %in% valid_status)
      stop("envelope does not match planned job ", j$job_id)
    if (!is.null(identity) &&
        (is.null(e$source_identity) || !identical(e$source_identity, identity)))
      stop("stale or absent source identity for ", j$job_id)
    attempt_rows[[i]] <- data.frame(job_id = j$job_id, scenario_id = j$scenario_id,
                                    rep = j$rep, seed = j$seed, status = e$status,
                                    stringsAsFactors = FALSE)
    r <- e$results
    if (!is.data.frame(r) || !all(row_required %in% names(r)))
      stop("malformed result rows for ", j$job_id)
    if (!nrow(r)) {
      if (e$status == "ok") stop("successful envelope without results")
      r <- data.frame(job_id = j$job_id, method = "unknown", estimand = "unknown",
                      status = e$status, truth = NA_real_, estimate = NA_real_,
                      se = NA_real_, lower = NA_real_, upper = NA_real_,
                      complete_reference = NA_real_, failure_reason =
                        if (is.null(e$failure_reason)) "job failed before estimand" else
                          as.character(e$failure_reason), stringsAsFactors = FALSE)
    }
    if (anyNA(r$job_id) || any(as.character(r$job_id) != j$job_id) ||
        anyNA(r$method) || any(!nzchar(as.character(r$method))) ||
        anyNA(r$estimand) || any(!nzchar(as.character(r$estimand))) ||
        anyNA(r$status) || any(!r$status %in% valid_status) ||
        anyDuplicated(paste(r$method, r$estimand, sep = "\034")))
      stop("invalid or duplicate result identity for ", j$job_id)
    for (field in intersect(c("scenario_id", "rep", "seed", "slice", "adapter"),
                            names(r))) {
      if (field %in% names(j) &&
          (anyNA(r[[field]]) || any(as.character(r[[field]]) !=
                                       as.character(j[[field]]))))
        stop("result provenance mismatch for ", j$job_id, ": ", field)
    }
    if (e$status == "ok" && any(r$status == "failed"))
      stop("successful envelope contains failed result")
    complete <- r$status == "ok"
    numeric_fields <- c("truth", "estimate", "se", "lower", "upper",
                        "complete_reference")
    if (!all(vapply(r[numeric_fields], is.numeric, logical(1))))
      stop("result numeric field type invalid")
    if (any(complete & !is.finite(r$estimate)))
      stop("non-finite completed estimate")
    if (any(complete & j$claim_tier == "inferential" &
            (!is.finite(r$truth) | !is.finite(r$complete_reference) |
             !is.finite(r$se) | r$se <= 0 | !is.finite(r$lower) |
             !is.finite(r$upper) | r$lower > r$estimate |
             r$estimate > r$upper)))
      stop("malformed inferential completion")
    has_interval <- complete & is.finite(r$lower) & is.finite(r$upper)
    if (any(has_interval & (r$lower > r$upper | r$lower > r$estimate |
                            r$estimate > r$upper)))
      stop("malformed completed interval")
    r$scenario_id <- j$scenario_id
    r$rep <- j$rep
    r$seed <- j$seed
    r$claim_tier <- j$claim_tier
    result_rows[[i]] <- r[, c("job_id", "scenario_id", "rep", "seed", "claim_tier",
                              row_required[-1]), drop = FALSE]
  }
  out <- list(attempts = do.call(rbind, attempt_rows),
              results = do.call(rbind, result_rows), jobs = jobs,
              source_identity = identity)
  rownames(out$attempts) <- NULL
  rownames(out$results) <- NULL
  out
}

.validation_mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) return(c(mean = NA_real_, mcse = NA_real_, lo = NA_real_, hi = NA_real_))
  avg <- mean(x)
  se <- if (n >= 2L) stats::sd(x) / sqrt(n) else NA_real_
  c(mean = avg, mcse = se, lo = avg - 1.96 * se, hi = avg + 1.96 * se)
}

# Compare arms on the same simulated replicate, not as independent samples.
# Differences are method_a minus method_b, with method names in lexical order.
paired_method_contrasts <- function(validated) {
  if (!is.list(validated) || !all(c("results", "jobs") %in% names(validated)) ||
      !is.data.frame(validated$results) || !is.data.frame(validated$jobs))
    stop("pass validated records")
  empty <- data.frame(scenario_id = character(), estimand = character(),
    method_a = character(), method_b = character(), n_attempts = integer(),
    n_a_ok = integer(), n_b_ok = integer(), n_pairs = integer(),
    n_missing_pairs = integer(), mean_difference = numeric(),
    difference_mcse = numeric(), difference_ci_low = numeric(),
    difference_ci_high = numeric(), stringsAsFactors = FALSE)
  r <- validated$results
  keys <- unique(r[, c("scenario_id", "estimand"), drop = FALSE])
  keys <- keys[order(keys$scenario_id, keys$estimand), , drop = FALSE]
  out <- list()
  for (i in seq_len(nrow(keys))) {
    scenario <- keys$scenario_id[[i]]; estimand <- keys$estimand[[i]]
    sub <- r[r$scenario_id == scenario & r$estimand == estimand, , drop = FALSE]
    methods <- sort(unique(as.character(sub$method)))
    if (length(methods) < 2L) next
    planned <- validated$jobs[validated$jobs$scenario_id == scenario, , drop = FALSE]
    for (pair in utils::combn(methods, 2L, simplify = FALSE)) {
      a <- sub[sub$method == pair[[1L]] & sub$status == "ok", , drop = FALSE]
      b <- sub[sub$method == pair[[2L]] & sub$status == "ok", , drop = FALSE]
      matched <- planned$job_id[planned$job_id %in% a$job_id &
                                  planned$job_id %in% b$job_id]
      difference <- a$estimate[match(matched, a$job_id)] -
        b$estimate[match(matched, b$job_id)]
      difference <- difference[is.finite(difference)]
      point <- .validation_mean_ci(difference)
      out[[length(out) + 1L]] <- data.frame(
        scenario_id = scenario, estimand = estimand,
        method_a = pair[[1L]], method_b = pair[[2L]],
        n_attempts = nrow(planned), n_a_ok = nrow(a), n_b_ok = nrow(b),
        n_pairs = length(difference), n_missing_pairs = nrow(planned) - length(difference),
        mean_difference = point["mean"], difference_mcse = point["mcse"],
        difference_ci_low = point["lo"], difference_ci_high = point["hi"],
        row.names = NULL)
    }
  }
  if (!length(out)) return(empty)
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

.validation_binom_ci <- function(k, n) {
  if (n < 1L) return(c(lo = NA_real_, hi = NA_real_))
  z <- 1.96
  p <- k / n
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(lo = center - half, hi = center + half)
}

.validation_joint_boot <- function(estimate, se, truth, key, B = 2000L) {
  n <- length(estimate)
  empty <- c(norm_lo = NA_real_, norm_hi = NA_real_, norm_mcse = NA_real_,
             ratio_lo = NA_real_, ratio_hi = NA_real_, ratio_mcse = NA_real_)
  if (n < 3L || any(!is.finite(estimate), !is.finite(se), !is.finite(truth)))
    return(empty)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  chars <- utf8ToInt(key)
  local_seed <- as.integer((sum(as.double(chars) * seq_along(chars)) + 20261018) %%
                             .Machine$integer.max)
  set.seed(local_seed)
  draws <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
  est <- matrix(estimate[draws], nrow = n)
  ses <- matrix(se[draws], nrow = n)
  err <- matrix((estimate - truth)[draws], nrow = n)
  em <- colMeans(est)
  esd <- sqrt(pmax(0, (colSums(est^2) - n * em^2) / (n - 1L)))
  ratio <- colMeans(ses) / esd
  norm <- colMeans(err) / esd
  ratio <- ratio[is.finite(ratio)]
  norm <- norm[is.finite(norm)]
  if (length(ratio) < B * 0.95 || length(norm) < B * 0.95) return(empty)
  c(norm_lo = unname(stats::quantile(norm, 0.025)),
    norm_hi = unname(stats::quantile(norm, 0.975)),
    norm_mcse = stats::sd(norm),
    ratio_lo = unname(stats::quantile(ratio, 0.025)),
    ratio_hi = unname(stats::quantile(ratio, 0.975)),
    ratio_mcse = stats::sd(ratio))
}

summarise_records <- function(validated) {
  if (!is.list(validated) || !all(c("attempts", "results", "jobs") %in% names(validated)))
    stop("pass validated records")
  r <- validated$results
  keys <- unique(r[, c("scenario_id", "method", "estimand", "claim_tier"), drop = FALSE])
  if (!nrow(keys)) stop("no result rows")
  out <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    k <- keys[i, ]
    sub <- r[r$scenario_id == k$scenario_id & r$method == k$method &
               r$estimand == k$estimand, , drop = FALSE]
    planned <- validated$jobs[validated$jobs$scenario_id == k$scenario_id, , drop = FALSE]
    n_attempts <- nrow(planned)
    ok <- sub[sub$status == "ok", , drop = FALSE]
    n_ok <- nrow(ok)
    n_failed <- sum(sub$status == "failed")
    n_unavailable <- sum(sub$status == "unavailable")
    n_missing <- n_attempts - nrow(sub)
    if (n_missing < 0L || n_ok + n_failed + n_unavailable + n_missing != n_attempts)
      stop("result key accounting does not match planned attempts")
    point <- .validation_mean_ci(ok$estimate)
    complete <- .validation_mean_ci(ok$complete_reference)
    interval_width <- .validation_mean_ci(ok$upper - ok$lower)
    paired_complete_diff <- .validation_mean_ci(ok$estimate - ok$complete_reference)
    truth_available <- n_ok > 0L && all(is.finite(ok$truth))
    empirical_sd <- if (n_ok >= 2L) stats::sd(ok$estimate) else NA_real_
    bias <- if (truth_available) .validation_mean_ci(ok$estimate - ok$truth) else
      c(mean = NA_real_, mcse = NA_real_, lo = NA_real_, hi = NA_real_)
    norm_bias <- bias["mean"] / empirical_sd
    se_ratio <- if (is.finite(empirical_sd) && empirical_sd > 0)
      mean(ok$se) / empirical_sd else NA_real_
    boot <- if (truth_available)
      .validation_joint_boot(ok$estimate, ok$se, ok$truth,
                             paste(k$scenario_id, k$method, k$estimand, sep = "/"))
    else c(norm_lo = NA_real_, norm_hi = NA_real_, norm_mcse = NA_real_,
           ratio_lo = NA_real_, ratio_hi = NA_real_, ratio_mcse = NA_real_)
    norm_bias_ci <- c(lo = unname(boot["norm_lo"]), hi = unname(boot["norm_hi"]))
    se_ratio_ci <- c(lo = unname(boot["ratio_lo"]), hi = unname(boot["ratio_hi"]))
    se_ratio_mcse <- boot["ratio_mcse"]
    covered <- if (truth_available && all(is.finite(ok$lower), is.finite(ok$upper)))
      ok$lower <= ok$truth & ok$truth <= ok$upper else logical(0)
    coverage <- if (length(covered)) mean(covered) else NA_real_
    coverage_mcse <- if (length(covered)) sqrt(coverage * (1 - coverage) / length(covered)) else NA_real_
    coverage_ci <- .validation_binom_ci(sum(covered), length(covered))
    classification <- "indeterminate"
    if (k$claim_tier != "inferential") {
      classification <- "descriptive"
    } else if (n_ok == 0L) {
      classification <- "indeterminate"
    } else if (!truth_available) {
      classification <- "descriptive"
    } else if (n_ok == n_attempts && n_ok >= 100L &&
               all(is.finite(c(norm_bias_ci, se_ratio_ci, coverage_ci)))) {
      if (coverage_ci["hi"] < 0.925 || coverage_ci["lo"] > 0.975 ||
          se_ratio_ci["hi"] < 0.9 || se_ratio_ci["lo"] > 1.1 ||
          norm_bias_ci["lo"] > 0.1 || norm_bias_ci["hi"] < -0.1) {
        classification <- "failed"
      } else if (coverage_ci["lo"] >= 0.925 && coverage_ci["hi"] <= 0.975 &&
                 se_ratio_ci["lo"] >= 0.9 && se_ratio_ci["hi"] <= 1.1 &&
                 norm_bias_ci["lo"] >= -0.1 && norm_bias_ci["hi"] <= 0.1) {
        classification <- "supported"
      }
    }
    out[[i]] <- data.frame(scenario_id = k$scenario_id, method = k$method,
       estimand = k$estimand, claim_tier = k$claim_tier,
       n_attempts = n_attempts, n_ok = n_ok, n_failed = n_failed,
       n_unavailable = n_unavailable, n_missing = n_missing,
       mean_estimate = point["mean"], mean_mcse = point["mcse"],
       mean_ci_low = point["lo"], mean_ci_high = point["hi"],
       complete_mean = complete["mean"], complete_mcse = complete["mcse"],
       mean_interval_width = interval_width["mean"],
       interval_width_mcse = interval_width["mcse"],
       interval_width_ci_low = interval_width["lo"],
       interval_width_ci_high = interval_width["hi"],
       paired_complete_diff = paired_complete_diff["mean"],
       paired_complete_diff_mcse = paired_complete_diff["mcse"],
       paired_complete_diff_ci_low = paired_complete_diff["lo"],
       paired_complete_diff_ci_high = paired_complete_diff["hi"],
       empirical_sd = empirical_sd, mean_se = if (n_ok) mean(ok$se) else NA_real_,
       bias = bias["mean"], bias_mcse = bias["mcse"],
       normalized_bias = norm_bias, normalized_bias_mcse = boot["norm_mcse"],
       normalized_bias_ci_low = norm_bias_ci["lo"],
       normalized_bias_ci_high = norm_bias_ci["hi"],
       se_ratio = se_ratio, se_ratio_mcse = se_ratio_mcse,
       se_ratio_ci_low = se_ratio_ci["lo"], se_ratio_ci_high = se_ratio_ci["hi"],
       coverage = coverage, coverage_mcse = coverage_mcse,
       coverage_ci_low = coverage_ci["lo"], coverage_ci_high = coverage_ci["hi"],
       classification = classification, row.names = NULL)
  }
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}
