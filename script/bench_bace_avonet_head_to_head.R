#!/usr/bin/env Rscript
#
# script/bench_bace_avonet_head_to_head.R
#
# Phase 8 MVP: pigauto v0.9.1.9000 vs BACE on bundled avonet300/tree300.
# Identical splits, seeds, and metrics. Single-seed for the MVP; multi-seed
# is a Phase 8.x follow-up.
#
# Output
#   script/bench_bace_avonet_head_to_head.rds
#   script/bench_bace_avonet_head_to_head.md

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_bace_avonet_head_to_head.rds")
out_md  <- file.path(here, "script", "bench_bace_avonet_head_to_head.md")

# -------------------------------------------------------------------------
# Data
# -------------------------------------------------------------------------

e <- new.env()
utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300",   package = "pigauto", envir = e)
df <- e$avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL
tree <- e$tree300

cat("avonet300:", nrow(df), "species x", ncol(df), "traits\n")

# -------------------------------------------------------------------------
# Split (identical across methods)
# -------------------------------------------------------------------------

SEED      <- 2026L
MISS_FRAC <- 0.30
pd0 <- preprocess_traits(df, tree, log_transform = TRUE)
splits <- make_missing_splits(pd0$X_scaled, missing_frac = MISS_FRAC,
                                val_frac = 0.5, seed = SEED,
                                trait_map = pd0$trait_map)

# Hold-out mask in the USER's data.frame (so all 3 methods see the same NAs).
# We derive this from the splits' test_idx linear positions on X_scaled, then
# map back to trait columns via trait_map.
mask_test <- matrix(FALSE, nrow = nrow(pd0$X_scaled), ncol = ncol(pd0$X_scaled))
mask_test[splits$test_idx] <- TRUE
mask_val  <- matrix(FALSE, nrow = nrow(pd0$X_scaled), ncol = ncol(pd0$X_scaled))
mask_val[splits$val_idx]  <- TRUE
# We evaluate on the test mask; val is eaten by pigauto for calibration.

# Map latent-col hits back to user-col hits via trait_map.
user_mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                          dimnames = list(rownames(df), names(df)))
for (tm in pd0$trait_map) {
  user_col <- tm$name
  if (!user_col %in% colnames(user_mask_test)) next
  latent_cols <- tm$latent_cols
  row_hits <- apply(mask_test[, latent_cols, drop = FALSE], 1L, any)
  user_mask_test[row_hits, user_col] <- TRUE
}

# Build the masked data.frame (truth preserved as df; mask applies only to test).
df_miss <- df
for (v in colnames(user_mask_test)) {
  df_miss[[v]][user_mask_test[, v]] <- NA
}

cat("test cells per trait (user-scale):\n")
print(colSums(user_mask_test))

# -------------------------------------------------------------------------
# Methods
# -------------------------------------------------------------------------

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(val = val, wall = proc.time()[["elapsed"]] - t0)
}

# pigauto_default (em_iterations = 0L)
run_pigauto_default <- function() {
  pigauto::impute(df_miss, tree, log_transform = TRUE,
                    missing_frac = 0, verbose = FALSE, seed = SEED,
                    epochs = 500L, em_iterations = 0L)
}

# pigauto_em5 (Phase 6 EM, diagonal)
run_pigauto_em5 <- function() {
  pigauto::impute(df_miss, tree, log_transform = TRUE,
                    missing_frac = 0, verbose = FALSE, seed = SEED,
                    epochs = 500L, em_iterations = 5L)
}

# BACE (optional — skipped if BACE is not installed)
run_bace <- function() {
  if (!requireNamespace("BACE", quietly = TRUE)) {
    return(NULL)
  }
  # BACE's canonical API on AVONET: bace() with default prior + OVR. The
  # precise call varies by BACE version — below is the v0.9.0 comparison
  # setup. If this errors on a newer BACE, the try() returns NULL and the
  # pipeline continues.
  tryCatch({
    BACE::bace(
      data          = df_miss,
      tree          = tree,
      n_iter        = 2000L,
      burnin        = 500L,
      thin          = 5L,
      ovr           = TRUE,
      verbose       = FALSE
    )
  }, error = function(e) {
    message("BACE run failed: ", conditionMessage(e))
    NULL
  })
}

# -------------------------------------------------------------------------
# Evaluation
# -------------------------------------------------------------------------

eval_completed <- function(completed_df, truth_df, mask, method) {
  rows <- list()
  for (v in colnames(mask)) {
    idx <- which(mask[, v])
    if (!length(idx)) next
    t_v <- truth_df[[v]][idx]
    c_v <- completed_df[[v]][idx]
    if (is.factor(t_v) || is.ordered(t_v) || is.character(t_v)) {
      acc <- mean(as.character(c_v) == as.character(t_v), na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "accuracy", value = acc,
        n_cells = length(idx)
      )
    } else {
      rmse <- sqrt(mean((as.numeric(t_v) - as.numeric(c_v))^2, na.rm = TRUE))
      pear <- suppressWarnings(stats::cor(as.numeric(t_v), as.numeric(c_v),
                                            use = "complete.obs"))
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "rmse", value = rmse,
        n_cells = length(idx)
      )
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "pearson_r", value = pear,
        n_cells = length(idx)
      )
    }
  }
  do.call(rbind, rows)
}

# -------------------------------------------------------------------------
# Run
# -------------------------------------------------------------------------

cat("\n=== pigauto_default ===\n")
r_def <- timed(run_pigauto_default())
ev_def <- eval_completed(r_def$val$completed, df, user_mask_test,
                           method = "pigauto_default")
ev_def$wall_s <- r_def$wall

cat("\n=== pigauto_em5 ===\n")
r_em5 <- timed(run_pigauto_em5())
ev_em5 <- eval_completed(r_em5$val$completed, df, user_mask_test,
                           method = "pigauto_em5")
ev_em5$wall_s <- r_em5$wall

cat("\n=== BACE (may skip) ===\n")
r_bace <- timed(run_bace())
ev_bace <- if (!is.null(r_bace$val)) {
  # BACE::bace() output has $data or similar — depends on version. Best
  # effort: try common slot names; if none match, skip.
  completed_bace <- tryCatch({
    if ("completed" %in% names(r_bace$val)) r_bace$val$completed
    else if ("data"  %in% names(r_bace$val)) r_bace$val$data
    else if ("imputed_data" %in% names(r_bace$val)) r_bace$val$imputed_data
    else NULL
  }, error = function(e) NULL)
  if (is.null(completed_bace)) {
    message("BACE output shape not recognised; skipping evaluation.")
    NULL
  } else {
    out <- eval_completed(completed_bace, df, user_mask_test,
                           method = "bace_default")
    out$wall_s <- r_bace$wall
    out
  }
} else NULL

all_rows <- do.call(rbind, Filter(Negate(is.null),
                                     list(ev_def, ev_em5, ev_bace)))

saveRDS(list(results = all_rows,
              seed = SEED, miss_frac = MISS_FRAC,
              bace_ran = !is.null(ev_bace)),
         out_rds)

# -------------------------------------------------------------------------
# Markdown
# -------------------------------------------------------------------------

md <- c(
  "# Phase 8 MVP: AVONET 300 head-to-head (pigauto vs BACE)",
  "",
  sprintf("Seed = %d, miss_frac = %.2f, identical splits across methods.",
          SEED, MISS_FRAC),
  if (!is.null(ev_bace)) "" else "**BACE skipped** (not installed or failed).",
  "",
  "## Per-trait metrics",
  "",
  "```",
  capture.output(print(all_rows, row.names = FALSE, max = 500L)),
  "```"
)
writeLines(md, out_md)
cat("\n=== DONE ===\n")
cat("  rds :", out_rds, "\n")
cat("  md  :", out_md, "\n")
