#!/usr/bin/env Rscript
# HTML generator for Phase 8.1 ρ sweep.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds <- file.path(here, "script", "bench_correlation_sweep.rds")
out_html_s <- file.path(here, "script", "bench_correlation_sweep.html")
out_html_p <- file.path(here, "pkgdown", "assets", "dev",
                          "bench_correlation_sweep.html")
dir.create(dirname(out_html_p), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_rds)) stop("run script/bench_correlation_sweep.R first")

b    <- readRDS(in_rds)
res  <- b$results
cfg  <- b$config

plot_df <- aggregate(value ~ method + rho + metric, data = res, FUN = mean)
p <- ggplot(plot_df, aes(x = rho, y = value,
                           colour = method, linetype = method)) +
  geom_line() + geom_point() +
  facet_wrap(~ metric, scales = "free_y") +
  scale_x_continuous(breaks = cfg$rhos) +
  labs(x = "cross-trait ρ", y = "metric value",
       title = "Phase 8.1: cross-trait-correlation sweep",
       subtitle = sprintf("n=%d, K=%d cont traits, miss=%.2f, reps=%d",
                          cfg$n_species, cfg$n_traits, cfg$miss_frac,
                          cfg$n_reps))
png_path <- file.path(here, "script", "bench_correlation_sweep.png")
ggsave(png_path, p, width = 9, height = 5, dpi = 120)

tbl_html <- paste0(
  "<table border='1' cellpadding='4' style='border-collapse:collapse'>\n",
  "<tr><th>method</th><th>ρ</th><th>metric</th><th>mean value</th></tr>\n",
  paste0(apply(plot_df, 1, function(r)
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%.3f</td></tr>",
            r[["method"]], r[["rho"]], r[["metric"]],
            as.numeric(r[["value"]]))),
    collapse = "\n"),
  "\n</table>")

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>Phase 8.1: ρ sweep</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1100px;",
  "margin:2em auto;padding:0 1em;color:#222}h1{border-bottom:2px solid #ccc}",
  "img{max-width:100%}</style></head><body>",
  "<h1>Phase 8.1: cross-trait-correlation (ρ) sweep</h1>",
  sprintf("<p>Config: n = %d, K = %d continuous traits, %d reps × %d ρ levels × %d methods.</p>",
          cfg$n_species, cfg$n_traits, cfg$n_reps, length(cfg$rhos),
          length(cfg$methods)),
  sprintf("<p>Total wall: %.1f min.</p>", b$script_wall / 60),
  sprintf("<img src='%s' alt='ρ curves'/>", basename(png_path)),
  "<h2>Summary table</h2>", tbl_html, "</body></html>")

writeLines(html, out_html_s)
writeLines(html, out_html_p)
cat("Wrote ", out_html_s, " + ", out_html_p, "\n")
