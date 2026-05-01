#!/usr/bin/env Rscript
#
# script/make_bench_signal_sweep_html.R
#
# Render the signal-sweep benchmark RDS into an HTML page linked from the
# pkgdown validation suite.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds  <- file.path(here, "script", "bench_signal_sweep.rds")
out_html_script   <- file.path(here, "script", "bench_signal_sweep.html")
out_html_pkgdown  <- file.path(here, "pkgdown", "assets", "dev",
                                 "bench_signal_sweep.html")
dir.create(dirname(out_html_pkgdown), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_rds)) {
  stop("Run script/bench_signal_sweep.R first to produce ", in_rds)
}

bench   <- readRDS(in_rds)
results <- bench$results
config  <- bench$config

# ---- Summary plot: one facet per (trait, metric) -----------------------
plot_data <- aggregate(value ~ method + lambda + trait + metric,
                        data = results, FUN = mean)

p <- ggplot(plot_data,
            aes(x = lambda, y = value, colour = method, linetype = method)) +
  geom_line() + geom_point() +
  facet_grid(metric ~ trait, scales = "free_y") +
  scale_x_continuous(breaks = config$lambdas) +
  labs(x = "Pagel's lambda", y = "metric value",
       title = "Phase 8 signal-strength sweep",
       subtitle = sprintf("n = %d, miss_frac = %.2f, n_reps = %d",
                          config$n_species, config$miss_frac, config$n_reps))

png_path <- file.path(here, "script", "bench_signal_sweep.png")
ggsave(png_path, p, width = 10, height = 7, dpi = 120)

# ---- HTML shell --------------------------------------------------------
render_table <- function(tbl) {
  paste0(
    "<table border='1' cellpadding='4' style='border-collapse:collapse'>\n",
    "<tr>", paste0("<th>", names(tbl), "</th>", collapse = ""), "</tr>\n",
    paste0(
      apply(tbl, 1L, function(r) {
        paste0("<tr>",
               paste0("<td>", r, "</td>", collapse = ""),
               "</tr>")
      }),
      collapse = "\n"
    ),
    "\n</table>"
  )
}

tbl_html <- render_table(
  cbind(plot_data[order(plot_data$trait, plot_data$metric,
                         plot_data$lambda, plot_data$method), ],
        rounded = round(plot_data[order(plot_data$trait, plot_data$metric,
                                         plot_data$lambda, plot_data$method),
                                   "value"], 3))
)

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>Phase 8 signal-sweep benchmark</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1100px;",
  "margin:2em auto;padding:0 1em;color:#222} h1{border-bottom:2px solid #ccc}",
  "img{max-width:100%}</style></head><body>",
  "<h1>Phase 8: Pagel's λ signal-strength sweep</h1>",
  sprintf("<p><strong>Config</strong>: n = %d species, %d reps × %d lambda levels × %d methods.</p>",
          config$n_species, config$n_reps, length(config$lambdas),
          length(config$methods)),
  sprintf("<p>Total wall time: %.1f min.</p>", bench$script_wall / 60),
  "<h2>Per-trait curves</h2>",
  sprintf("<img src='%s' alt='signal-sweep curves'/>", basename(png_path)),
  "<h2>Means per (method, lambda, trait, metric)</h2>",
  tbl_html,
  "</body></html>"
)

writeLines(html, out_html_script)
writeLines(html, out_html_pkgdown)
cat("Wrote", out_html_script, "\n")
cat("Wrote", out_html_pkgdown, "\n")
