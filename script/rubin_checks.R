# script/rubin_checks.R
#
# Gate checks for the rubin-freq-bace lane (.unlazy/rubin-freq-bace/gates/leaf-cell.md). Each gate prints
# "<gate> PASS" or stops loudly.
#
#   Rscript script/rubin_checks.R --gate G-S4a --dir <smoke results dir>
#   Rscript script/rubin_checks.R --gate G-S4b

suppressPackageStartupMessages({ library(ape) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1L] }
gate <- get_arg("--gate"); dir <- get_arg("--dir")
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))

if (identical(gate, "G-S4a")) {
  f <- list.files(dir, "^rubin_.*\\.rds$", full.names = TRUE)
  stopifnot("no rubin rds in --dir" = length(f) >= 1L)
  x <- readRDS(f[1])
  want <- c("freqA", "freqB", "bace", "bace_resid")
  stopifnot("errors recorded" = length(x$errors) == 0L,
            "missing arm in cells" = all(want %in% x$cells$arm),
            "missing arm in estimands" = all(want %in% x$estimands$arm))
  cc <- x$cells
  stopifnot("non-finite per-cell score" = all(is.finite(cc$zRMSE) & is.finite(cc$coverage) & is.finite(cc$width)),
            "zero-width intervals outside bace (as shipped)" = all(cc$frac_B0[cc$arm != "bace"] == 0))
  stopifnot("complete-data reference row missing" = sum(x$estimands$arm == "complete") == 2L)
  e <- x$estimands[x$estimands$arm %in% want, ]
  stopifnot("both estimands per arm" = all(table(e$arm)[want] == 2L),
            "non-finite pooled estimand" = all(is.finite(e$estimate) & is.finite(e$lower) & is.finite(e$upper)),
            "FMI missing or outside [0, 1]" = all(is.finite(e$fmi) & e$fmi >= 0 & e$fmi <= 1),
            "every imputation analysed" = all(e$m_ok == x$M))
  stopifnot("BACE convergence/ESS not recorded" = !is.null(x$diag$bace) && length(x$diag$bace) > 0L)
  stopifnot("partial NA imputations" = all(cc$n_na == 0L))
  cat(sprintf("smoke %s: %d cell rows, %d estimand rows; walls %s\n", basename(f[1]), nrow(cc), nrow(e),
              paste(names(x$walls), round(x$walls), sep = "=", collapse = " ")))
  cat("G-S4a PASS\n")
}

if (identical(gate, "G-S4b")) {
  source(file.path(here, "campaign_gnn_off_lib.R")); source(file.path(here, "rubin_lib.R"))
  cov_run <- function(n, reps = 200L, lambda = 0.7, rho = 0.5) {
    hit <- t(vapply(seq_len(reps), function(i) {
      cell <- make_cell("types_mixed", n, 7000L + i, lambda = lambda, rho = rho)
      s <- est_pgls_slope(cell$truth, cell$tree); co <- est_phylo_cor(cell$truth, cell$tree)
      q <- stats::qt(0.975, s$df_com)
      c(slope = abs(s$estimate - rho) <= q * sqrt(s$variance),
        cor = abs(co$z - atanh(rho)) <= stats::qnorm(0.975) * sqrt(co$variance))
    }, numeric(2)))
    colMeans(hit, na.rm = TRUE)
  }
  r60 <- cov_run(60L); r100 <- cov_run(100L)
  cat(sprintf("complete-data coverage of rho: n=60 slope %.3f cor %.3f | n=100 slope %.3f cor %.3f\n",
              r60["slope"], r60["cor"], r100["slope"], r100["cor"]))
  all_cov <- c(r60, r100)
  if (any(all_cov < 0.92 | all_cov > 0.98)) stop("complete-data coverage outside [0.92, 0.98]")
  cat("G-S4b PASS\n")
}
