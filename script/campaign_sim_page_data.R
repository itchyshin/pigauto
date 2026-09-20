#!/usr/bin/env Rscript
# Turn the aggregator's csv output into one compact JSON for the results page.
#
#   Rscript script/campaign_sim_page_data.R --agg /tmp/pig_pool2/agg --out /tmp/sim_data.json
#
# Reads <agg>_summary.csv and <agg>_paired.csv. Emits only what the page draws, so the page stays
# small and needs no library to parse it. Every number keeps its regime and its MCSE.

args <- commandArgs(trailingOnly = TRUE)
get <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1L] }
agg <- get("--agg"); out <- get("--out", "sim_data.json")
if (is.null(agg)) stop("usage: --agg <prefix> [--out file.json]")

s <- read.csv(paste0(agg, "_summary.csv"), stringsAsFactors = FALSE)
pth <- paste0(agg, "_paired.csv")
p <- if (file.exists(pth)) read.csv(pth, stringsAsFactors = FALSE) else NULL

# minimal JSON writer: no jsonlite dependency on a cluster login node
esc <- function(x) gsub('"', '\\\\"', as.character(x))
jval <- function(v) {
  if (is.numeric(v)) ifelse(is.finite(v), formatC(v, digits = 6, format = "g"), "null")
  else paste0('"', esc(v), '"')
}
jrows <- function(d) {
  if (is.null(d) || !nrow(d)) return("[]")
  cols <- names(d)
  body <- vapply(seq_len(nrow(d)), function(i)
    paste0("{", paste0('"', cols, '":', vapply(cols, function(k) jval(d[i, k]), ""), collapse = ","), "}"), "")
  paste0("[", paste(body, collapse = ","), "]")
}

keep_metrics <- c("zRMSE", "accuracy", "macroF1", "brier", "coverage", "width", "interval_score")
s <- s[s$metric %in% keep_metrics, ]

# trait type, so the page can pool continuous-family and discrete separately
cont <- c("c1", "c2", "c3", "c4", "cnt", "prp")
s$kind <- ifelse(s$trait %in% cont, "continuous", "discrete")

# pooled over traits within (cell, arm, metric): the primary contrast is stated pooled over types
pool <- do.call(rbind, lapply(
  split(s, list(s$n, s$lambda, s$rho, s$evo, s$miss, s$arm, s$metric, s$kind), drop = TRUE),
  function(d) data.frame(
    n = d$n[1], lambda = d$lambda[1], rho = d$rho[1], evo = d$evo[1], miss = d$miss[1],
    arm = d$arm[1], metric = d$metric[1], kind = d$kind[1],
    reps = round(mean(d$reps)), traits = nrow(d),
    mean = mean(d$mean, na.rm = TRUE),
    # MCSE of a mean over traits, treating traits as independent summaries: conservative, and the
    # page labels it as pooled. Per-trait numbers keep their own MCSE in `byTrait`.
    mcse = sqrt(sum(d$mcse^2, na.rm = TRUE)) / max(1, sum(is.finite(d$mcse))),
    stringsAsFactors = FALSE)))

meta <- list(generated = format(Sys.time(), "%Y-%m-%d %H:%M"),
             cells = length(unique(paste(s$n, s$lambda, s$rho, s$evo, s$miss))),
             arms = paste(sort(unique(s$arm)), collapse = ","),
             maxReps = max(s$reps, na.rm = TRUE))

json <- paste0(
  '{"meta":{"generated":"', meta$generated, '","cells":', meta$cells,
  ',"arms":"', meta$arms, '","maxReps":', meta$maxReps, '},',
  '"pooled":', jrows(pool), ',',
  '"byTrait":', jrows(s[, c("n","lambda","rho","evo","miss","arm","trait","kind","metric","reps","mean","mcse")]), ',',
  '"paired":', jrows(if (is.null(p)) NULL else p[, intersect(names(p),
      c("n","lambda","rho","evo","miss","trait","metric","arm","reference","n_common_seeds","mean_diff","mcse_diff"))]),
  '}')
writeLines(json, out)
cat("wrote", out, "-", nrow(pool), "pooled rows,", nrow(s), "per-trait rows,",
    if (is.null(p)) 0 else nrow(p), "paired rows\n")
