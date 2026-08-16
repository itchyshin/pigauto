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
                    epochs = 500L, em_iterations = 0L,
                    n_imputations = 20L)
}

# pigauto_em5 (Phase 6 EM, diagonal)
run_pigauto_em5 <- function() {
  pigauto::impute(df_miss, tree, log_transform = TRUE,
                    missing_frac = 0, verbose = FALSE, seed = SEED,
                    epochs = 500L, em_iterations = 5L,
                    n_imputations = 20L)
}

# BACE (optional — skipped if BACE is not installed). Correct
# BACE::bace() call -- uses fixformula + ran_phylo_form + phylo + data
# (with a Species column), same pattern as bench_avonet_bace.R.
run_bace <- function() {
  if (!requireNamespace("BACE", quietly = TRUE)) {
    return(list(completed = NULL, error = "BACE not installed"))
  }
  # BACE (via MCMCglmm) refuses phylogenies with zero-length edges.
  tree_b <- tree
  if (any(tree_b$edge.length == 0, na.rm = TRUE)) {
    tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
  }

  df_b <- df_miss
  df_b$Species <- rownames(df_miss)
  all_traits <- setdiff(names(df_b), "Species")
  fixformula <- lapply(all_traits, function(v) {
    others <- setdiff(all_traits, v)
    paste0(v, " ~ ", paste(others, collapse = " + "))
  })

  bace_error <- NULL
  out <- tryCatch({
    BACE::bace(
      fixformula     = fixformula,
      ran_phylo_form = "~ 1 |Species",
      phylo          = tree_b,
      data           = df_b,
      nitt           = 2000L,
      burnin         = 500L,
      thin           = 5L,
      runs           = 2L,
      n_final        = 2L,
      verbose        = FALSE,
      skip_conv      = TRUE
    )
  }, error = function(e) {
    bace_error <<- conditionMessage(e)
    message("BACE run failed: ", bace_error)
    NULL
  })

  if (is.null(out)) return(list(completed = NULL, error = bace_error))

  # BACE returns a list with $imputed_datasets (list of M completed
  # data.frames). Pool by median / mode across draws for the final
  # completed dataset.
  imputed_sets <- tryCatch({
    if ("imputed_datasets" %in% names(out)) out$imputed_datasets
    else if ("imputed_data" %in% names(out)) out$imputed_data
    else if ("data" %in% names(out)) list(out$data)
    else NULL
  }, error = function(e) NULL)

  if (is.null(imputed_sets) || !length(imputed_sets)) {
    shape_msg <- sprintf("BACE output shape not recognised (keys: %s)",
                          paste(names(out), collapse = ", "))
    message(shape_msg)
    return(list(completed = NULL, error = shape_msg))
  }

  M <- length(imputed_sets)
  completed <- df_miss
  for (v in names(completed)) {
    if (!any(user_mask_test[, v])) next
    draws <- sapply(imputed_sets, function(d) d[[v]])
    if (!is.matrix(draws)) draws <- matrix(draws, ncol = M)
    if (is.factor(completed[[v]])) {
      for (i in which(user_mask_test[, v])) {
        vals <- as.character(draws[i, ])
        vals <- vals[!is.na(vals)]
        if (!length(vals)) next
        mode_val <- names(sort(table(vals), decreasing = TRUE))[1]
        completed[[v]][i] <- factor(mode_val, levels = levels(completed[[v]]),
                                     ordered = is.ordered(completed[[v]]))
      }
    } else {
      for (i in which(user_mask_test[, v])) {
        vals <- as.numeric(draws[i, ])
        completed[[v]][i] <- stats::median(vals, na.rm = TRUE)
      }
    }
  }

  list(completed = completed, error = NULL)
}

# -------------------------------------------------------------------------
# Evaluation
# -------------------------------------------------------------------------

eval_completed <- function(completed_df, truth_df, mask, method,
                            res_obj = NULL) {
  lo      <- if (!is.null(res_obj)) res_obj$prediction$conformal_lower else NULL
  hi      <- if (!is.null(res_obj)) res_obj$prediction$conformal_upper else NULL
  mi_list <- if (!is.null(res_obj)) res_obj$prediction$imputed_datasets else NULL

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
      t_num <- as.numeric(t_v)
      if (!is.null(lo) && !is.null(hi) && v %in% colnames(lo)) {
        lo_v <- lo[idx, v]; hi_v <- hi[idx, v]
        valid <- is.finite(lo_v) & is.finite(hi_v) & is.finite(t_num)
        if (any(valid)) {
          hits <- t_num[valid] >= lo_v[valid] & t_num[valid] <= hi_v[valid]
          rows[[length(rows) + 1L]] <- data.frame(
            method = method, trait = v, metric = "coverage95_conformal",
            value = mean(hits), n_cells = sum(valid))
        }
      }
      if (!is.null(mi_list) && length(mi_list) > 1L && v %in% names(mi_list[[1]])) {
        draws_mat <- vapply(mi_list, function(d) as.numeric(d[idx, v]),
                             numeric(length(idx)))
        if (!is.matrix(draws_mat))
          draws_mat <- matrix(draws_mat, ncol = length(mi_list))
        if (nrow(draws_mat) == length(idx) && ncol(draws_mat) > 1L) {
          q_lo <- apply(draws_mat, 1L, stats::quantile, probs = 0.025, na.rm = TRUE)
          q_hi <- apply(draws_mat, 1L, stats::quantile, probs = 0.975, na.rm = TRUE)
          valid <- is.finite(q_lo) & is.finite(q_hi) & is.finite(t_num)
          if (any(valid)) {
            hits <- t_num[valid] >= q_lo[valid] & t_num[valid] <= q_hi[valid]
            rows[[length(rows) + 1L]] <- data.frame(
              method = method, trait = v, metric = "coverage95_mcdropout",
              value = mean(hits), n_cells = sum(valid))
          }
        }
      }
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
                           method = "pigauto_default",
                           res_obj = r_def$val)
ev_def$wall_s <- r_def$wall

cat("\n=== pigauto_em5 ===\n")
r_em5 <- timed(run_pigauto_em5())
ev_em5 <- eval_completed(r_em5$val$completed, df, user_mask_test,
                           method = "pigauto_em5",
                           res_obj = r_em5$val)
ev_em5$wall_s <- r_em5$wall

cat("\n=== BACE (may skip) ===\n")
r_bace <- timed(run_bace())
ev_bace <- if (!is.null(r_bace$val$completed)) {
  out <- eval_completed(r_bace$val$completed, df, user_mask_test,
                         method = "bace_default")
  out$wall_s <- r_bace$wall
  out
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

bace_skip_msg <- if (!is.null(ev_bace)) "" else {
  err <- r_bace$val$error
  if (is.null(err) || !nzchar(err)) err <- "not installed or failed (no error captured)"
  sprintf("**BACE skipped**: %s", err)
}

md <- c(
  "# Phase 8 MVP: AVONET 300 head-to-head (pigauto vs BACE)",
  "",
  sprintf("Seed = %d, miss_frac = %.2f, identical splits across methods.",
          SEED, MISS_FRAC),
  bace_skip_msg,
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
