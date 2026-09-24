# script/rubin_prerun_summary.R
#
# Summarise the BACE settings pre-run (docs/dev-log/arc/2026-09-24-rubin-prerun-plan.md) and apply its
# selection rule. Reads every rds under <dir>/runs<r>_nitt<nitt>/ written by script/rubin_cell.R.
#
#   Rscript script/rubin_prerun_summary.R --dir <prerun dir> [--expected 20]
#
# Per setting (runs, nitt) and n: fits found vs expected (a missing rds = timeout or crash, which the runner
# cannot record), BACE fit errors, BACE's convergence pass share, median ESS, median wall times of the fit and
# the chain, and, for bace and bace_chain, mean per-cell coverage of c1 and c2 and slope / correlation
# coverage of rho, each with a Monte Carlo SE computed across replicates (one value per fit).
# Selection rule: the cheapest setting with >= 80% convergence passes, median ESS >= 100, and bace_chain
# per-cell coverage within 2 MCSE of the most expensive setting's. No qualifying setting is reported as such.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1L] }
dir <- get_arg("--dir"); expected <- as.integer(get_arg("--expected", 20L))
stopifnot("--dir is required" = !is.null(dir))

sets <- list.dirs(dir, recursive = FALSE)
sets <- sets[grepl("^runs[0-9]+_nitt[0-9]+$", basename(sets))]
stopifnot("no runs<r>_nitt<nitt> directories" = length(sets) > 0L)

mcse <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x)) }
rows <- list()
for (sd in sets) {
  runs <- as.integer(sub("^runs([0-9]+)_nitt.*$", "\\1", basename(sd)))
  nitt <- as.integer(sub("^.*_nitt([0-9]+)$", "\\1", basename(sd)))
  fs <- list.files(sd, "^rubin_.*\\.rds$", full.names = TRUE)
  if (!length(fs)) next
  xs <- lapply(fs, readRDS)
  for (n in sort(unique(vapply(xs, `[[`, numeric(1), "n")))) {
    xn <- Filter(function(x) x$n == n, xs)
    fit_err <- vapply(xn, function(x) !is.null(x$errors$bace_fit), logical(1))
    ok <- xn[!fit_err]
    conv <- vapply(ok, function(x) isTRUE(x$diag$bace$converged), logical(1))
    ess <- vapply(ok, function(x) as.numeric(x$diag$bace$ess_med %||% NA), numeric(1))
    wall <- function(k) stats::median(vapply(ok, function(x) as.numeric(x$walls[k] %||% NA), numeric(1)), na.rm = TRUE)
    arm_cov <- function(arm) {
      cell <- vapply(ok, function(x) { cc <- x$cells[x$cells$arm == arm & x$cells$trait %in% c("c1", "c2"), ]
                                       if (!nrow(cc)) NA_real_ else stats::weighted.mean(cc$coverage, cc$n_cells) }, numeric(1))
      est <- function(e) vapply(ok, function(x) { r <- x$estimands[x$estimands$arm == arm & x$estimands$estimand == e, ]
                                                   if (!nrow(r)) NA_real_ else as.numeric(r$covered) }, numeric(1))
      sl <- est("slope"); co <- est("cor")
      c(cell = mean(cell, na.rm = TRUE), cell_mcse = mcse(cell), slope = mean(sl, na.rm = TRUE),
        slope_mcse = mcse(sl), cor = mean(co, na.rm = TRUE), cor_mcse = mcse(co))
    }
    b <- arm_cov("bace"); ch <- arm_cov("bace_chain")
    rows[[length(rows) + 1L]] <- data.frame(
      runs = runs, nitt = nitt, n = n, cost_units = nitt * (runs + 40), found = length(xn), expected = expected,
      fit_errors = sum(fit_err), conv_pass = if (length(ok)) mean(conv) else NA_real_,
      ess_med = stats::median(ess, na.rm = TRUE), wall_fit_s = wall("bace_fit"), wall_chain_s = wall("bace_chain"),
      bace_cell = b[["cell"]], bace_cell_mcse = b[["cell_mcse"]], chain_cell = ch[["cell"]],
      chain_cell_mcse = ch[["cell_mcse"]], bace_slope = b[["slope"]], chain_slope = ch[["slope"]],
      bace_cor = b[["cor"]], chain_cor = ch[["cor"]])
  }
}
tab <- do.call(rbind, rows)
tab <- tab[order(tab$n, tab$cost_units), ]
print(tab, digits = 3, row.names = FALSE)

cat("\nSelection rule (per n):\n")
for (n in sort(unique(tab$n))) {
  t <- tab[tab$n == n, ]
  ref <- t[which.max(t$cost_units), ]
  qual <- t[is.finite(t$conv_pass) & t$conv_pass >= 0.8 & is.finite(t$ess_med) & t$ess_med >= 100 &
            abs(t$chain_cell - ref$chain_cell) <= 2 * pmax(t$chain_cell_mcse, ref$chain_cell_mcse, na.rm = TRUE), ]
  if (!nrow(qual)) cat(sprintf("  n = %d: no setting meets the rule\n", n)) else {
    q <- qual[which.min(qual$cost_units), ]
    cat(sprintf("  n = %d: runs %d, nitt %d (conv %.2f, ESS %.0f, chain cell coverage %.3f)\n",
                n, q$runs, q$nitt, q$conv_pass, q$ess_med, q$chain_cell))
  }
}
incomplete <- tab[tab$found < tab$expected, c("runs", "nitt", "n", "found", "expected")]
if (nrow(incomplete)) { cat("\nIncomplete settings (missing rds = timeout or crash):\n"); print(incomplete, row.names = FALSE) }
