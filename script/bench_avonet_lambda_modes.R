# AVONET300: does pigauto's own Pagel-lambda mode close the BACE continuous-trait gap?
#
# Motivated by two same-day results:
#   docs/dev-log/2026-08-16-lambda-attribution-results.md  -- lambda_mode="bayes" gives
#     the best baseline repair on simulated data and eliminates the ML boundary collapse.
#   script/bench_avonet_bace.md -- BACE (MCMCglmm; Bayesian averaging over the
#     phylo/residual variance ratio, i.e. over what pigauto calls lambda) beats pigauto
#     on all four AVONET continuous traits while losing on categorical.
#
# Hypothesis: pigauto's continuous-trait deficit vs BACE is the lambda = 1 default, and
# lambda_mode = "bayes" recovers part of it. This is the real-data confirmation the
# simulated campaign explicitly said was still needed.
#
# Same data, seed, and mask as script/bench_avonet_bace.R so the numbers are comparable
# cell-for-cell with the BACE run recorded there.

SEED <- 2026L; MISS_FRAC <- 0.30
MODES <- c("fixed_1", "bayes")
out_rds <- "script/bench_avonet_lambda_modes.rds"
out_md  <- "script/bench_avonet_lambda_modes.md"

suppressPackageStartupMessages({ library(pigauto); library(ape) })
try(torch::torch_set_num_threads(1L), silent = TRUE)
t0 <- proc.time()[["elapsed"]]
log_line <- function(...) { cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
                                ..., "\n", sep = ""); flush.console() }

e <- new.env(parent = emptyenv())
utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300",   package = "pigauto", envir = e)
df <- e$avonet300; tree <- e$tree300
rownames(df) <- df$Species_Key; df$Species_Key <- NULL
stopifnot(all(rownames(df) == tree$tip.label))

# identical mask construction to bench_avonet_bace.R
set.seed(SEED)
df_truth <- df
mask_test <- matrix(FALSE, nrow(df), ncol(df), dimnames = list(rownames(df), names(df)))
for (v in names(df)) {
  obs_idx <- which(!is.na(df[[v]]))
  mask_test[sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC)), v] <- TRUE
}
df_miss <- df
for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA

safe_cor <- function(x, y) {
  i <- which(is.finite(x) & is.finite(y))
  if (length(i) < 2L || stats::sd(x[i]) == 0 || stats::sd(y[i]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[i], y[i]))
}
eval_completed <- function(completed, method, wall_s) {
  rows <- list()
  for (v in colnames(mask_test)) {
    idx <- which(mask_test[, v]); if (!length(idx)) next
    t_v <- df_truth[[v]][idx]; c_v <- completed[[v]][idx]
    if (is.factor(t_v) || is.ordered(t_v) || is.character(t_v)) {
      rows[[length(rows)+1L]] <- data.frame(method=method, trait=v, metric="accuracy",
        value=mean(as.character(c_v) == as.character(t_v), na.rm=TRUE),
        n_cells=length(idx), wall_s=wall_s)
    } else {
      rows[[length(rows)+1L]] <- data.frame(method=method, trait=v, metric="rmse",
        value=sqrt(mean((as.numeric(t_v)-as.numeric(c_v))^2, na.rm=TRUE)),
        n_cells=length(idx), wall_s=wall_s)
      rows[[length(rows)+1L]] <- data.frame(method=method, trait=v, metric="pearson_r",
        value=safe_cor(as.numeric(t_v), as.numeric(c_v)),
        n_cells=length(idx), wall_s=wall_s)
    }
  }
  do.call(rbind, rows)
}

results <- list()
for (md in MODES) {
  log_line("=== lambda_mode = ", md, " ===")
  w0 <- proc.time()[["elapsed"]]
  res <- tryCatch(pigauto::impute(df_miss, tree, lambda_mode = md, verbose = FALSE,
                                  seed = SEED, n_imputations = 20L),
                  error = function(e) e)
  w <- proc.time()[["elapsed"]] - w0
  if (inherits(res, "error")) {
    log_line("FAILED: ", conditionMessage(res)); next
  }
  results[[md]] <- eval_completed(res$completed, paste0("pigauto_", md), w)
  log_line(sprintf("done in %.1f s", w))
}
out <- do.call(rbind, results)
saveRDS(list(results = out, seed = SEED, miss_frac = MISS_FRAC), out_rds)

md_lines <- c("# AVONET300: pigauto lambda_mode fixed_1 vs bayes",
  "", sprintf("n = %d species x %d traits, seed = %d, MCAR %.0f%%, m = 20.",
              nrow(df), ncol(df), SEED, 100*MISS_FRAC),
  "Same data/seed/mask as `script/bench_avonet_bace.md`, so BACE's numbers there are",
  "directly comparable. Single seed -- a direction, not an interval.", "", "```",
  capture.output(print(out, row.names = FALSE)), "```")
writeLines(md_lines, out_md)
log_line("wrote ", out_rds, " and ", out_md)
print(out, row.names = FALSE)
