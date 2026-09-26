# script/rubin_study/build_data.R
#
# Collects every per-fit result of the freq-vs-BACE Rubin campaign into flat tables under
# script/rubin_study/data/, next to the campaign aggregates. Run from the worktree root after
# script/rubin_pool.sh and script/rubin_campaign_aggregate.R:
#   Rscript script/rubin_study/build_data.R [POOL_DIR]
# Outputs (one row per ...):
#   fits.csv                fit (dataset x arm set): cell, seed, cluster, node, time, wall times, versions, errors
#   fit_cells.csv.gz        fit x arm x trait: per-value coverage, width, z-RMSE, interval score
#   fit_estimands.csv.gz    fit x arm x estimand: pooled estimate, SE, interval, df, FMI, truth, covered
#   bace_diagnostics.csv    BACE fit: convergence verdict, drift, ESS, input cleaning
#   agg/                    copies of the aggregator's tables (cells, down, paired, failures, report.json)
a <- commandArgs(trailingOnly = TRUE)
pool <- if (length(a)) a[1] else file.path(Sys.getenv("HOME"), "pigauto_rubin_pool")
out <- file.path("script", "rubin_study", "data")
dir.create(file.path(out, "agg"), recursive = TRUE, showWarnings = FALSE)

files <- list.files(file.path(pool, c("bace", "freq")), "\\.rds$", recursive = TRUE, full.names = TRUE)
set_of <- function(f) basename(dirname(dirname(f)))            # bace | freq
cluster_of <- function(f) basename(dirname(f))                 # nibi | fir | rorqual | totoro
# the pool can hold the same fit from two hosts (duplicate submissions); keep one per (set, tag)
key <- paste(set_of(files), sub("_dup[0-9]+\\.rds$", ".rds", basename(files)))
files <- files[!duplicated(key)]

pkgv <- function(si, p) {
  v <- si$otherPkgs[[p]]$Version %||% si$loadedOnly[[p]]$Version
  if (is.null(v)) NA_character_ else v
}
`%||%` <- function(x, y) if (is.null(x)) y else x

fits <- vector("list", length(files)); cells <- fits; est <- fits; diag <- fits
for (i in seq_along(files)) {
  x <- readRDS(files[i]); si <- x$sessionInfo
  id <- data.frame(set = set_of(files[i]), tag = x$tag, n = x$n, lambda = x$lambda, rho = x$rho, seed = x$seed,
                   stringsAsFactors = FALSE)
  w <- x$walls %||% numeric(0)
  fits[[i]] <- cbind(id, data.frame(
    cluster = cluster_of(files[i]), node = x$host %||% NA, time = format(x$time %||% NA),
    M = x$M, arms = paste(x$arms, collapse = ";"), realised_frac = x$realised_frac,
    wall_s_total = sum(w), wall_s = paste(names(w), round(w, 1), sep = "=", collapse = ";"),
    n_errors = length(x$errors), errors = paste(names(x$errors), collapse = ";"),
    R = if (!is.null(si)) paste(si$R.version$major, si$R.version$minor, sep = ".") else NA,
    BACE = pkgv(si, "BACE"), MCMCglmm = pkgv(si, "MCMCglmm"), Rphylopars = pkgv(si, "Rphylopars"),
    git_hash = x$git_hash %||% NA, stringsAsFactors = FALSE))
  if (!is.null(x$cells) && nrow(x$cells)) cells[[i]] <- cbind(id, x$cells)
  if (!is.null(x$estimands) && nrow(x$estimands)) est[[i]] <- cbind(id, x$estimands)
  d <- x$diag$bace
  if (!is.null(d)) diag[[i]] <- cbind(id, data.frame(
    converged = d$converged %||% NA, n_attempts = d$n_attempts %||% NA, drift = d$drift %||% NA,
    ess_min = d$ess_min %||% NA, ess_med = d$ess_med %||% NA, ess_frac_low = d$ess_frac_low %||% NA,
    input_fix = paste(d$input_fix %||% character(0), collapse = ";"), stringsAsFactors = FALSE))
}
bind <- function(l) do.call(rbind, l[!vapply(l, is.null, logical(1))])
F <- bind(fits); C <- bind(cells); E <- bind(est); D <- bind(diag)
utils::write.csv(F, file.path(out, "fits.csv"), row.names = FALSE)
utils::write.csv(C, gzfile(file.path(out, "fit_cells.csv.gz")), row.names = FALSE)
utils::write.csv(E, gzfile(file.path(out, "fit_estimands.csv.gz")), row.names = FALSE)
utils::write.csv(D, file.path(out, "bace_diagnostics.csv"), row.names = FALSE)

# as-shipped BACE failures (first pass, before input cleaning): one row per record
fl <- list.files(file.path(pool, c("bace_asshipped_failed", "bace_failed_later")), "\\.rds$", recursive = TRUE,
                 full.names = TRUE)
FR <- do.call(rbind, lapply(fl, function(f) { x <- readRDS(f)
  data.frame(folder = basename(dirname(dirname(f))), cluster = basename(dirname(f)), tag = x$tag %||% basename(f),
             n = x$n %||% NA, lambda = x$lambda %||% NA, rho = x$rho %||% NA, seed = x$seed %||% NA,
             errors = paste(vapply(x$errors %||% list(), function(e) substr(as.character(e), 1, 120), ""),
                            collapse = " | "), stringsAsFactors = FALSE) }))
if (!is.null(FR)) utils::write.csv(FR, file.path(out, "failure_records.csv"), row.names = FALSE)

agg <- file.path(pool, "agg")
for (f in c("cells.csv", "cells_n.csv", "cells_l.csv", "down.csv", "down_n.csv", "down_l.csv", "paired.csv",
            "failures.csv", "report.json", "sanity.txt"))
  file.copy(file.path(agg, f), file.path(out, "agg", f), overwrite = TRUE)

cat(sprintf("fits %d (bace %d, freq %d); cell rows %d; estimand rows %d; BACE diagnostics %d; failure records %d\n",
            nrow(F), sum(F$set == "bace"), sum(F$set == "freq"), nrow(C), nrow(E), nrow(D), NROW(FR)))
