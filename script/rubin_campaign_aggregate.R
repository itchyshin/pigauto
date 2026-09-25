# script/rubin_campaign_aggregate.R
#
# Aggregate the pooled rubin-freq-bace campaign (script/rubin_pool.sh) into the tables behind the final report.
#   Rscript script/rubin_campaign_aggregate.R [POOL_DIR] [OUT_DIR]
#
# Unit of replication: one simulated dataset (n, lambda, rho, seed). Frequentist and BACE arms were run on the
# same datasets (same seeds, same make_cell), so paired differences are taken per dataset. Every Monte Carlo SE
# is over datasets. Per-cell coverage is pooled over c1 and c2 (cell-weighted within a dataset); prp is
# reported separately on the logit scale.
# Outputs (CSV + one JSON for the page): cells.csv, down.csv, paired.csv, failures.csv, sanity.txt, report.json.

suppressMessages(library(jsonlite))
a <- commandArgs(trailingOnly = TRUE)
pool <- if (length(a) >= 1) a[1] else file.path(Sys.getenv("HOME"), "pigauto_rubin_pool")
out <- if (length(a) >= 2) a[2] else file.path(pool, "agg")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

read_dir <- function(d) {
  fs <- list.files(d, "^rubin_.*\\.rds$", recursive = TRUE, full.names = TRUE)
  fs <- fs[!duplicated(basename(fs))]                         # one copy per dataset x arm set
  lapply(fs, function(f) { x <- tryCatch(readRDS(f), error = function(e) NULL); if (!is.null(x)) x$file <- basename(f); x })
}
B <- Filter(Negate(is.null), read_dir(file.path(pool, "bace")))
F <- Filter(Negate(is.null), read_dir(file.path(pool, "freq")))
A <- Filter(Negate(is.null), read_dir(file.path(pool, "bace_asshipped_failed")))
key <- function(x) sprintf("n%d_l%s_r%s_s%d", x$n, format(x$lambda), format(x$rho), x$seed)

# ---- long tables ----------------------------------------------------------------------------------
cells_long <- function(L) do.call(rbind, lapply(L, function(x) if (!is.null(x$cells)) cbind(ds = key(x), n = x$n, lambda = x$lambda,
                                                   rho = x$rho, seed = x$seed, x$cells)))
est_long <- function(L) do.call(rbind, lapply(L, function(x) if (!is.null(x$estimands)) cbind(ds = key(x), n = x$n, lambda = x$lambda,
                                                   rho = x$rho, seed = x$seed, x$estimands)))
cl <- rbind(cells_long(B), cells_long(F)); el <- rbind(est_long(B), est_long(F))
cl <- cl[!(cl$arm == "complete"), ]
# the complete-data reference must agree between the two arm sets on the same dataset (same simulated data)
cb <- el[el$arm == "complete" & el$ds %in% vapply(B, key, ""), ]; cf <- el[el$arm == "complete" & el$ds %in% vapply(F, key, ""), ]
mm <- merge(cb[!duplicated(paste(cb$ds, cb$estimand)), c("ds", "estimand", "estimate")],
            cf[!duplicated(paste(cf$ds, cf$estimand)), c("ds", "estimand", "estimate")], by = c("ds", "estimand"))
sanity <- sprintf("complete-data estimate, BACE files vs freq files, same dataset: %d pairs, max |diff| %.2e",
                  nrow(mm), if (nrow(mm)) max(abs(mm$estimate.x - mm$estimate.y)) else NA)
el <- el[!duplicated(paste(el$ds, el$arm, el$estimand)), ]

# ---- per-cell (c1, c2 pooled within a dataset) ----------------------------------------------------------
cc <- cl[cl$trait %in% c("c1", "c2"), ]
per_ds <- aggregate(cbind(cv = coverage * n_cells, wd = width * n_cells, zr = zRMSE * n_cells, is = interval_score * n_cells, w = n_cells) ~
                      ds + n + lambda + rho + arm, data = cc, FUN = sum)
per_ds <- transform(per_ds, coverage = cv / w, width = wd / w, zrmse = zr / w, iscore = is / w)
se <- function(x) sd(x) / sqrt(length(x))
summ <- function(d, by, vars) do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(g) {
  r <- g[1, by, drop = FALSE]
  for (v in vars) { r[[v]] <- mean(g[[v]]); r[[paste0(v, "_se")]] <- se(g[[v]]) }
  r$width_median <- if ("width" %in% names(g)) median(g$width) else NA
  r$datasets <- nrow(g); r }))
cells <- summ(per_ds, c("arm", "n", "lambda", "rho"), c("coverage", "width", "zrmse", "iscore"))
cells_n <- summ(per_ds, c("arm", "n"), c("coverage", "width", "zrmse", "iscore"))
cells_l <- summ(per_ds, c("arm", "n", "lambda"), c("coverage", "width", "zrmse", "iscore"))
prp <- cl[cl$trait == "prp", ]
prp_s <- if (nrow(prp)) summ(transform(prp, iscore = interval_score, zrmse = zRMSE), c("arm", "n"), c("coverage", "width", "zrmse")) else NULL

# ---- downstream ---------------------------------------------------------------------------------
dd <- transform(el, err_rho = estimate - truth, err_cd = estimate - complete_data, ciw = upper - lower, covered = as.numeric(covered))
down <- do.call(rbind, lapply(split(dd, dd[c("arm", "n", "lambda", "rho", "estimand")], drop = TRUE), function(g) data.frame(
  g[1, c("arm", "n", "lambda", "rho", "estimand")], coverage = mean(g$covered), coverage_se = se(g$covered),
  bias = mean(g$err_rho), bias_se = se(g$err_rho), rmse = sqrt(mean(g$err_rho^2)), bias_vs_complete = mean(g$err_cd),
  ci_width = mean(g$ciw), fmi_median = if (all(is.na(g$fmi))) NA else median(g$fmi, na.rm = TRUE), datasets = nrow(g), row.names = NULL)))
down_l <- do.call(rbind, lapply(split(dd, dd[c("arm", "n", "lambda", "estimand")], drop = TRUE), function(g) data.frame(
  g[1, c("arm", "n", "lambda", "estimand")], coverage = mean(g$covered), coverage_se = se(g$covered), bias = mean(g$err_rho),
  rmse = sqrt(mean(g$err_rho^2)), ci_width = mean(g$ciw), fmi_median = if (all(is.na(g$fmi))) NA else median(g$fmi, na.rm = TRUE),
  datasets = nrow(g), row.names = NULL)))
down_n <- do.call(rbind, lapply(split(dd, dd[c("arm", "n", "estimand")], drop = TRUE), function(g) data.frame(
  g[1, c("arm", "n", "estimand")], coverage = mean(g$covered), coverage_se = se(g$covered), bias = mean(g$err_rho), bias_se = se(g$err_rho),
  rmse = sqrt(mean(g$err_rho^2)), ci_width = mean(g$ciw), fmi_median = if (all(is.na(g$fmi))) NA else median(g$fmi, na.rm = TRUE),
  datasets = nrow(g), row.names = NULL)))

# ---- paired differences on the same datasets ------------------------------------------------------------
pair <- function(a1, a2, what) {
  if (what == "cell") { x <- per_ds[per_ds$arm == a1, c("ds", "n", "lambda", "coverage")]; y <- per_ds[per_ds$arm == a2, c("ds", "coverage")]; v <- "coverage" }
  else { x <- dd[dd$arm == a1 & dd$estimand == what, c("ds", "n", "lambda", "covered")]; y <- dd[dd$arm == a2 & dd$estimand == what, c("ds", "covered")]; v <- "covered" }
  m <- merge(x, y, by = "ds"); if (!nrow(m)) return(NULL)
  m$d <- m[[paste0(v, ".x")]] - m[[paste0(v, ".y")]]
  do.call(rbind, lapply(split(m, m[c("n", "lambda")], drop = TRUE), function(g) data.frame(
    contrast = paste(a1, "-", a2), measure = what, n = g$n[1], lambda = g$lambda[1], diff = mean(g$d), se = se(g$d), datasets = nrow(g))))
}
paired <- do.call(rbind, c(
  lapply(list(c("bace_chain", "bace"), c("freqA", "bace_chain"), c("freqA", "bace"), c("freqA", "freqB")), function(p) pair(p[1], p[2], "cell")),
  lapply(list(c("bace_chain", "bace"), c("freqA", "bace_chain"), c("freqA", "bace")), function(p) pair(p[1], p[2], "slope")),
  lapply(list(c("bace_chain", "bace"), c("freqA", "bace_chain"), c("freqA", "bace")), function(p) pair(p[1], p[2], "cor"))))

# ---- failures ledger ------------------------------------------------------------------------------------
fl <- function(L, what) if (length(L)) do.call(rbind, lapply(L, function(x) data.frame(set = what, n = x$n, lambda = x$lambda,
  bace_error = !is.null(x$errors$bace_fit), freq_error = any(grepl("^freq", names(x$errors))),
  cleaned = length(x$diag$bace$input_fix) > 0, fA_fail = x$diag$freqA$n_fail %||% NA, fA_degen = x$diag$freqA$n_degenerate %||% NA)))
fr <- rbind(fl(B, "bace"), fl(F, "freq"), fl(A, "bace_asshipped_failed"))
failures <- do.call(rbind, lapply(split(fr, fr[c("set", "n", "lambda")], drop = TRUE), function(g) data.frame(
  set = g$set[1], n = g$n[1], lambda = g$lambda[1], files = nrow(g), bace_errors = sum(g$bace_error), freq_errors = sum(g$freq_error),
  bace_cleaned = sum(g$cleaned), freqA_refit_fail = sum(g$fA_fail, na.rm = TRUE), freqA_degenerate = sum(g$fA_degen, na.rm = TRUE), row.names = NULL)))

for (nm in c("cells", "cells_n", "cells_l", "down", "down_n", "down_l", "paired", "failures")) write.csv(get(nm), file.path(out, paste0(nm, ".csv")), row.names = FALSE)
writeLines(c(sanity, sprintf("BACE files %d, freq files %d, as-shipped failure records %d", length(B), length(F), length(A))), file.path(out, "sanity.txt"))
writeLines(toJSON(list(cells = cells, cells_n = cells_n, cells_l = cells_l, prp = prp_s, down = down, down_n = down_n, down_l = down_l, paired = paired, failures = failures,
                       sanity = sanity, counts = list(bace = length(B), freq = length(F), asshipped_failed = length(A))),
                  digits = 5, na = "null", auto_unbox = TRUE, dataframe = "rows"), file.path(out, "report.json"))
cat(readLines(file.path(out, "sanity.txt")), sep = "\n"); print(cells_n, digits = 3, row.names = FALSE); print(down_n, digits = 3, row.names = FALSE)
