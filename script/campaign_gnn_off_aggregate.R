# script/campaign_gnn_off_aggregate.R
# Aggregate campaign cells (rds from campaign_gnn_off_cell.R OR campaign_sim_cell.R) into per
# (dgp, n, arm, trait) tables with mean and Monte Carlo SE over seeds, wall-time summaries, coverage,
# dispatch paths, and error counts. Section I additions (only fire when the new-shape fields are
# present -- i.e. rds written by campaign_sim_cell.R): per (cell, arm, trait, metric) MCSE, a paired
# difference vs `--reference` with MCSE of the per-replicate difference (common seeds only), pooled
# ECE (10 equal-mass bins over all masked cells of a cell) with a bootstrap-over-replicates MCSE,
# failure rate, and summary.csv / paired.csv / failures.csv / raw_index.csv + sessionInfo().
# Usage: Rscript script/campaign_gnn_off_aggregate.R <results_dir> <out_prefix> [--reference ARM]
args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
reference_arm <- get_flag("--reference", "floor")
pos <- args[!grepl("^--", args)]; pos <- pos[!(pos %in% c(reference_arm))]
res_dir <- pos[1]; out <- pos[2]
fs <- list.files(res_dir, pattern = "\\.rds$", full.names = TRUE)
fs <- fs[!grepl("_smoke", fs)]
cells <- lapply(fs, readRDS)
cat(length(cells), "cells read\n")
tab <- do.call(rbind, lapply(cells, function(c) { r <- c$results; if (is.null(r) || !nrow(r)) return(NULL); r$dgp <- c$dgp; r$n <- c$n; r$seed <- c$seed; r }))
walls <- do.call(rbind, lapply(cells, function(c) data.frame(dgp = c$dgp, n = c$n, seed = c$seed, arm = names(c$walls), wall = as.numeric(unlist(c$walls)))))
errs <- do.call(rbind, lapply(cells, function(c) if (length(c$errors)) data.frame(dgp = c$dgp, n = c$n, seed = c$seed, arm = names(c$errors), error = unlist(c$errors)) else NULL))
paths <- do.call(rbind, lapply(cells, function(c) { p <- c$paths$gnn_off; if (is.null(p)) return(NULL); data.frame(dgp = c$dgp, n = c$n, seed = c$seed, trait = names(p), path = as.character(p)) }))
mcse <- function(x) stats::sd(x) / sqrt(length(x))
# per trait
per_trait <- do.call(rbind, lapply(split(tab, list(tab$dgp, tab$n, tab$arm, tab$trait, tab$metric), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], trait = d$trait[1], metric = d$metric[1], n_seeds = nrow(d),
             mean = mean(d$value, na.rm = TRUE), mcse = mcse(d$value[!is.na(d$value)]),
             coverage = if (all(is.na(d$coverage))) NA_real_ else mean(d$coverage, na.rm = TRUE),
             coverage_mcse = if (all(is.na(d$coverage))) NA_real_ else mcse(d$coverage[!is.na(d$coverage)]))))
# per cell-arm summary: mean over continuous traits (zRMSE) and discrete traits (accuracy), then over seeds
cell_arm <- do.call(rbind, lapply(split(tab, list(tab$dgp, tab$n, tab$seed, tab$arm), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], seed = d$seed[1], arm = d$arm[1],
             zrmse = mean(d$value[d$metric == "zRMSE"], na.rm = TRUE),
             acc = if (any(d$metric == "accuracy" & !is.na(d$value))) mean(d$value[d$metric == "accuracy"], na.rm = TRUE) else NA_real_,
             coverage = if (all(is.na(d$coverage))) NA_real_ else mean(d$coverage, na.rm = TRUE))))
summary_arm <- do.call(rbind, lapply(split(cell_arm, list(cell_arm$dgp, cell_arm$n, cell_arm$arm), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], n_seeds = nrow(d),
             zrmse = mean(d$zrmse), zrmse_mcse = mcse(d$zrmse),
             acc = mean(d$acc, na.rm = TRUE), acc_mcse = if (all(is.na(d$acc))) NA_real_ else mcse(d$acc[!is.na(d$acc)]),
             coverage = mean(d$coverage, na.rm = TRUE))))
summary_arm <- summary_arm[order(summary_arm$dgp, summary_arm$n, summary_arm$zrmse), ]
wall_sum <- aggregate(wall ~ dgp + n + arm, walls, function(x) c(mean = mean(x), max = max(x)))
path_tab <- if (!is.null(paths)) as.data.frame(table(paths$dgp, paths$n, paths$trait, paths$path)) else NULL
saveRDS(list(tab = tab, per_trait = per_trait, cell_arm = cell_arm, summary_arm = summary_arm, walls = walls,
             wall_sum = wall_sum, errors = errs, paths = paths, n_cells = length(cells)), paste0(out, ".rds"))
write.csv(summary_arm, paste0(out, "_summary_arm.csv"), row.names = FALSE)
write.csv(per_trait, paste0(out, "_per_trait.csv"), row.names = FALSE)
cat("errors:", if (is.null(errs)) 0 else nrow(errs), "\n")
print(summary_arm, digits = 3, row.names = FALSE)

# ---- Section I: campaign_sim_cell.R rds shape (lambda/rho/evo/miss/calib/failed present) ------
is_sim_shape <- length(cells) > 0 && !is.null(cells[[1]]$lambda)
if (is_sim_shape) {
  regime_of <- function(c) data.frame(lambda = c$lambda %||% NA, rho = c$rho %||% NA, evo = c$evo %||% NA,
                                       miss = c$miss %||% NA, thresholds = c$thresholds %||% NA)
  `%||%` <- function(a, b) if (is.null(a)) b else a
  tab2 <- do.call(rbind, lapply(cells, function(c) {
    r <- c$results; if (is.null(r) || !nrow(r)) return(NULL)
    cbind(r, regime_of(c))
  }))
  mcse_boot <- function(x, R = 500L) {
    x <- x[is.finite(x)]; if (length(x) < 2) return(NA_real_)
    b <- replicate(R, mean(sample(x, length(x), replace = TRUE)))
    stats::sd(b)
  }
  cell_key <- c("dgp", "n", "lambda", "rho", "evo", "miss")
  summary_i <- do.call(rbind, lapply(split(tab2, list(tab2$dgp, tab2$n, tab2$lambda, tab2$rho, tab2$evo, tab2$miss, tab2$arm, tab2$trait, tab2$metric), drop = TRUE), function(d)
    data.frame(d[1, cell_key], arm = d$arm[1], trait = d$trait[1], metric = d$metric[1], reps = nrow(d),
               mean = mean(d$value, na.rm = TRUE), mcse = mcse(d$value[!is.na(d$value)]))))
  write.csv(summary_i, paste0(out, "_summary.csv"), row.names = FALSE)

  # paired difference vs reference_arm, on common seeds only
  paired_rows <- list()
  for (grp in split(tab2, list(tab2$dgp, tab2$n, tab2$lambda, tab2$rho, tab2$evo, tab2$miss, tab2$trait, tab2$metric), drop = TRUE)) {
    ref <- grp[grp$arm == reference_arm, ]
    if (!nrow(ref)) next
    for (a in setdiff(unique(grp$arm), reference_arm)) {
      cur <- grp[grp$arm == a, ]
      m <- merge(ref[c("seed", "value")], cur[c("seed", "value")], by = "seed", suffixes = c(".ref", ".arm"))
      if (!nrow(m)) next
      d <- m$value.arm - m$value.ref
      paired_rows[[length(paired_rows) + 1L]] <- data.frame(grp[1, cell_key], trait = grp$trait[1], metric = grp$metric[1],
                                                              arm = a, reference = reference_arm, n_common_seeds = nrow(m),
                                                              mean_diff = mean(d), mcse_diff = mcse(d))
    }
  }
  paired <- do.call(rbind, paired_rows)
  if (!is.null(paired)) write.csv(paired, paste0(out, "_paired.csv"), row.names = FALSE)

  # pooled ECE per (cell, arm), 10 equal-mass bins over all masked cells, bootstrap MCSE over replicates
  calib_all <- do.call(rbind, lapply(cells, function(c) { cb <- c$calib; if (is.null(cb) || !nrow(cb)) return(NULL); cbind(cb, regime_of(c), dgp = c$dgp, n = c$n, seed = c$seed) }))
  ece <- NULL
  if (!is.null(calib_all)) {
    ece_one <- function(conf, correct, n_bins = 10L) {
      ord <- order(conf); conf <- conf[ord]; correct <- correct[ord]
      bins <- cut(seq_along(conf), breaks = n_bins, labels = FALSE)
      sum(vapply(seq_len(n_bins), function(b) {
        idx <- bins == b; if (!any(idx)) return(0)
        (sum(idx) / length(conf)) * abs(mean(correct[idx]) - mean(conf[idx]))
      }, numeric(1)))
    }
    ece <- do.call(rbind, lapply(split(calib_all, list(calib_all$dgp, calib_all$n, calib_all$lambda, calib_all$rho,
                                                        calib_all$evo, calib_all$miss, calib_all$arm), drop = TRUE), function(d) {
      by_seed <- split(d, d$seed)
      per_seed_ece <- vapply(by_seed, function(s) ece_one(s$confidence, s$correct), numeric(1))
      data.frame(d[1, cell_key], arm = d$arm[1], n_reps = length(per_seed_ece),
                 ece = mean(per_seed_ece), ece_mcse = mcse_boot(per_seed_ece))
    }))
    write.csv(ece, paste0(out, "_ece.csv"), row.names = FALSE)
  }

  # failure rate per (cell, arm)
  fail_rows <- do.call(rbind, lapply(cells, function(c) {
    if (!length(c$failed)) return(NULL)
    data.frame(regime_of(c), dgp = c$dgp, n = c$n, seed = c$seed, arm = names(c$failed))
  }))
  failure_rate <- if (!is.null(fail_rows)) do.call(rbind, lapply(split(fail_rows, list(fail_rows$dgp, fail_rows$n, fail_rows$arm), drop = TRUE), function(d)
    data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], n_failed = nrow(d)))) else NULL
  write.csv(if (is.null(failure_rate)) data.frame() else failure_rate, paste0(out, "_failures.csv"), row.names = FALSE)

  raw_index <- do.call(rbind, lapply(seq_along(cells), function(i) {
    c <- cells[[i]]
    data.frame(file = basename(fs[i]), dgp = c$dgp, n = c$n, seed = c$seed, lambda = c$lambda %||% NA,
               rho = c$rho %||% NA, evo = c$evo %||% NA, miss = c$miss %||% NA, host = c$host %||% NA,
               arms = paste(c$arms, collapse = ";"), n_failed = length(c$failed))
  }))
  write.csv(raw_index, paste0(out, "_raw_index.csv"), row.names = FALSE)
  saveRDS(utils::sessionInfo(), paste0(out, "_sessionInfo.rds"))
  cat("Section I outputs written:", paste0(out, c("_summary.csv", "_paired.csv", "_ece.csv", "_failures.csv", "_raw_index.csv")), "\n")
}
