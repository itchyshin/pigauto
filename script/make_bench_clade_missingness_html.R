#!/usr/bin/env Rscript
# HTML generator for Phase 8.3 clade-correlated missingness.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto"
in_rds <- file.path(here, "script", "bench_clade_missingness.rds")
out_html_s <- file.path(here, "script", "bench_clade_missingness.html")
out_html_p <- file.path(here, "pkgdown", "assets", "dev",
                          "bench_clade_missingness.html")
dir.create(dirname(out_html_p), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(in_rds)) stop("run script/bench_clade_missingness.R first")

b    <- readRDS(in_rds)
res  <- b$results
cfg  <- b$config

plot_df <- aggregate(value ~ method + target_frac + trait + metric,
                      data = res, FUN = mean)

# Focus plot: focal trait across target_frac
focal_df <- plot_df[plot_df$trait == cfg$focal_trait, ]
p <- ggplot(focal_df, aes(x = target_frac, y = value,
                            colour = method, linetype = method)) +
  geom_line() + geom_point() +
  facet_wrap(~ metric, scales = "free_y") +
  scale_x_continuous(breaks = cfg$target_fracs) +
  labs(x = "target missing fraction (clade-correlated)",
       y = "metric value",
       title = sprintf("Phase 8.3: clade-correlated MAR on %s",
                        cfg$focal_trait),
       subtitle = sprintf("n=%d, %d reps, mixed-type data", cfg$n_species,
                          cfg$n_reps))
png_path <- file.path(here, "script", "bench_clade_missingness.png")
ggsave(png_path, p, width = 9, height = 5, dpi = 120)

tbl_html <- paste0(
  "<table border='1' cellpadding='4' style='border-collapse:collapse'>\n",
  "<tr><th>method</th><th>target_frac</th><th>trait</th><th>metric</th>",
  "<th>mean value</th></tr>\n",
  paste0(apply(plot_df, 1, function(r)
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%.3f</td></tr>",
            r[["method"]], r[["target_frac"]], r[["trait"]], r[["metric"]],
            as.numeric(r[["value"]]))),
    collapse = "\n"),
  "\n</table>")

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>Phase 8.3: clade-MAR</title>",
  "<style>body{font-family:system-ui,sans-serif;max-width:1100px;",
  "margin:2em auto;padding:0 1em;color:#222}h1{border-bottom:2px solid #ccc}",
  "img{max-width:100%}</style></head><body>",
  "<h1>Phase 8.3: clade-correlated missingness</h1>",
  sprintf("<p>Pick random internal nodes, union their descendants until a target fraction is covered, mask ONE focal trait (%s) in all covered tips. n=%d, %d target_fracs × %d reps × %d methods.</p>",
          cfg$focal_trait, cfg$n_species,
          length(cfg$target_fracs), cfg$n_reps, length(cfg$methods)),
  sprintf("<p>Total wall: %.1f min.</p>", b$script_wall / 60),
  sprintf("<h2>Focal trait (%s) curves</h2>", cfg$focal_trait),
  sprintf("<img src='%s' alt='clade-MAR focal curves'/>", basename(png_path)),
  "<h2>Full summary table</h2>", tbl_html,
  "</body></html>")

writeLines(html, out_html_s)
writeLines(html, out_html_p)
cat("Wrote ", out_html_s, " + ", out_html_p, "\n")
