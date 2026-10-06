#!/usr/bin/env Rscript
# Usage: Rscript 02_summarise_masked_confirmation.R run_dir summary.rds
# Fail-closed: any retained method receipt with status != "ok" stops the run.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("expected: run_dir summary.rds", call. = FALSE)

alpha <- 0.05
winkler <- function(truth, lo, hi, alpha) {
  w <- hi - lo
  below <- truth < lo
  above <- truth > hi
  w + (2 / alpha) * (lo - truth) * below + (2 / alpha) * (truth - hi) * above
}

x <- lapply(c("split", "mondrian"), function(m) readRDS(file.path(args[[1L]], paste0(m, ".rds"))))
names(x) <- c("split", "mondrian")
if (!all(vapply(x, function(z) identical(z$status, "ok"), logical(1)))) {
  stop("retained method error: ", paste(vapply(x, function(z) z$status, character(1)), collapse = ", "), call. = FALSE)
}

tab <- do.call(rbind, Map(function(z, method) transform(z$metrics, method = method), x, names(x)))

mondrian_info <- x$mondrian$mondrian
n_val  <- if (!is.null(mondrian_info)) vapply(mondrian_info, function(m) m$n_val %||% NA_integer_, numeric(1)) else NULL
n_near <- if (!is.null(mondrian_info)) vapply(mondrian_info, function(m) m$n_near %||% NA_integer_, numeric(1)) else NULL
n_far  <- if (!is.null(mondrian_info)) vapply(mondrian_info, function(m) m$n_far %||% NA_integer_, numeric(1)) else NULL
fb     <- if (!is.null(mondrian_info)) vapply(mondrian_info, function(m) isTRUE(m$fallback), logical(1)) else NULL

`%||%` <- function(a, b) if (is.null(a)) b else a

per_stratum <- function(cells, method) {
  if (is.null(cells) || !nrow(cells)) return(NULL)
  cells$stratum[is.na(cells$stratum)] <- "all"
  agg <- do.call(rbind, lapply(split(cells, list(cells$trait, cells$stratum), drop = TRUE), function(g) {
    trait <- g$trait[1L]; stratum <- g$stratum[1L]
    n_test <- nrow(g)
    cov <- mean(g$truth >= g$lo & g$truth <= g$hi)
    hw <- (g$hi - g$lo) / 2
    n_s <- if (!is.null(n_val)) {
      if (identical(stratum, "near") && !is.null(n_near)) n_near[[trait]] %||% NA_real_
      else if (identical(stratum, "far") && !is.null(n_far)) n_far[[trait]] %||% NA_real_
      else n_val[[trait]] %||% NA_real_
    } else NA_real_
    mcse <- sqrt(cov * (1 - cov) / n_test + alpha * (1 - alpha) / (n_s + 2))
    data.frame(
      method = method, trait = trait, stratum = stratum, n_test = n_test,
      coverage = cov, median_half_width = stats::median(hw),
      winkler = mean(winkler(g$truth, g$lo, g$hi, alpha)),
      fallback = if (!is.null(fb)) isTRUE(fb[[trait]]) else NA,
      n_val = n_s, n_near = if (!is.null(n_near)) n_near[[trait]] %||% NA_real_ else NA_real_,
      n_far = if (!is.null(n_far)) n_far[[trait]] %||% NA_real_ else NA_real_,
      mcse = mcse, stringsAsFactors = FALSE
    )
  }))
  rownames(agg) <- NULL
  agg
}

strat_split <- per_stratum(x$split$cells, "split")
strat_mond  <- per_stratum(x$mondrian$cells, "mondrian")
strat <- rbind(strat_split, strat_mond)

paired <- NULL
if (!is.null(strat_split) && !is.null(strat_mond)) {
  key <- c("trait", "stratum")
  m <- merge(strat_mond, strat_split, by = key, suffixes = c("_mondrian", "_split"))
  m$coverage_gain <- m$coverage_mondrian - m$coverage_split
  m$width_ratio <- m$median_half_width_mondrian / m$median_half_width_split
  paired <- m[, c("trait", "stratum", "coverage_gain", "width_ratio",
                   "n_test_mondrian", "n_test_split")]
}

out <- list(methods = x, metrics = tab, stratum = strat, paired = paired)
saveRDS(out, args[[2L]])
print(strat)
if (!is.null(paired)) print(paired)
