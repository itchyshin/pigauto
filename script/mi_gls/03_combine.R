#!/usr/bin/env Rscript
# Combine the main MI-GLS sweep (GNN-based methods) with the fast rerun
# (complete, conditional-BM draws with EM and plug-in Sigma, oracle) into one
# per-cell table, then summarise per regime x method x downstream.
# The fast rerun uses the identical per-cell datasets (same seeds; data are
# simulated before any method runs), so rows combine cell by cell.
#
# Usage: Rscript script/mi_gls/03_combine.R <main_dir> <fast_dir> <out.csv> <out.md>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("expected: main_dir fast_dir out.csv out.md", call. = FALSE)
main_dir <- args[[1L]]; fast_dir <- args[[2L]]; out_csv <- args[[3L]]; out_md <- args[[4L]]

read_dir <- function(d, keep = NULL, rename = NULL) {
  f <- list.files(d, "^regime_[0-9]+_rep_[0-9]+[.]rds$", full.names = TRUE)
  do.call(rbind, lapply(f, function(p) {
    x <- readRDS(p); r <- x$results
    if (!is.null(keep)) r <- r[r$method %in% keep, , drop = FALSE]
    if (!is.null(rename)) r$method <- ifelse(r$method %in% names(rename), rename[r$method], r$method)
    if (!nrow(r)) return(NULL)
    cbind(regime_id = x$regime_id, rep = x$rep, x$regime[rep(1, nrow(r)), c("lambda", "n", "mechanism", "missing")],
          r[, c("method", "downstream", "estimate", "se", "covered")], row.names = NULL)
  }))
}
# Main sweep: GNN-based methods only; its draw_cond used the plug-in Sigma and
# is superseded by the fast rerun's draw_cond_inhouse on the same data.
main <- read_dir(main_dir, keep = c("single", "mi_conf_per_column", "mi_conf_exact", "mi_dropout"))
fast <- read_dir(fast_dir, keep = c("complete", "draw_cond", "draw_cond_inhouse", "oracle"),
                 rename = c(draw_cond = "draw_cond_em"))
d <- rbind(main, fast)
utils::write.csv(d, out_csv, row.names = FALSE)

true_beta <- 0.7
g <- split(d, list(d$regime_id, d$method, d$downstream), drop = TRUE)
s <- do.call(rbind, lapply(g, function(z) {
  est <- z$estimate[is.finite(z$estimate)]; R <- length(est)
  cov <- mean(z$covered, na.rm = TRUE)
  data.frame(regime_id = z$regime_id[1], lambda = z$lambda[1], n = z$n[1],
             mechanism = z$mechanism[1], missing = z$missing[1],
             method = z$method[1], downstream = z$downstream[1], R = R,
             bias = mean(est) - true_beta, bias_mcse = stats::sd(est) / sqrt(R),
             emp_sd = stats::sd(est), mean_se = mean(z$se, na.rm = TRUE),
             se_ratio = mean(z$se, na.rm = TRUE) / stats::sd(est),
             coverage = cov, coverage_mcse = sqrt(cov * (1 - cov) / R))
}))
s <- s[order(s$regime_id, s$downstream, s$method), ]
fmt <- function(v) ifelse(is.numeric(v), formatC(v, digits = 3, format = "f"), v)
lines <- c("# MI under phylogenetic GLS: per-regime summary", "",
           sprintf("Cells: %d rows from %s (GNN methods) and %s (draw methods, complete, oracle).",
                   nrow(d), main_dir, fast_dir), "",
           paste0("| ", paste(names(s), collapse = " | "), " |"),
           paste0("|", paste(rep("---", ncol(s)), collapse = "|"), "|"))
for (i in seq_len(nrow(s))) lines <- c(lines, paste0("| ", paste(vapply(s[i, ], function(v) {
  if (is.numeric(v) && !is.integer(v)) formatC(v, digits = 3, format = "f") else as.character(v) }, ""), collapse = " | "), " |"))
writeLines(lines, out_md)
cat("wrote", out_csv, "and", out_md, "(", nrow(s), "summary rows )\n")
