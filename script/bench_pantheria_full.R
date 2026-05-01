#!/usr/bin/env Rscript
#
# script/bench_pantheria_full.R
#
# Full-scale PanTHERIA bench: pigauto_default vs pigauto_em5 vs
# mean_baseline on the aligned PanTHERIA / mammal-tree intersection.
# Single seed (2026); 30% MCAR per trait; 8 canonical traits.
#
# Run after: Rscript script/fetch_pantheria_and_tree.R
#
# Output
#   script/bench_pantheria_full.{rds,md}

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_pantheria_full.rds")
out_md  <- file.path(here, "script", "bench_pantheria_full.md")

cache_dir <- file.path(here, "script", "data-cache")
pantheria_local <- file.path(cache_dir, "pantheria.txt")
tree_local      <- file.path(cache_dir, "mammal_tree.tre")

if (!file.exists(pantheria_local) || !file.exists(tree_local)) {
  stop("Data cache missing. Run script/fetch_pantheria_and_tree.R first.",
       call. = FALSE)
}

# -------------------------------------------------------------------------
# Load + align
# -------------------------------------------------------------------------

cat("Loading PanTHERIA + tree...\n")
pan <- utils::read.table(pantheria_local, header = TRUE, sep = "\t",
                          na.strings = "-999",
                          stringsAsFactors = FALSE, quote = "",
                          comment.char = "")
pan$species_key <- paste(pan$MSW93_Genus, pan$MSW93_Species, sep = "_")
pan$species_key <- gsub("[^A-Za-z0-9_]", "", pan$species_key)

tree <- ape::read.tree(tree_local)
tree$tip.label <- gsub("[^A-Za-z0-9_]", "", tree$tip.label)

overlap <- intersect(tree$tip.label, pan$species_key)
cat(sprintf("Aligned: %d species (of %d PanTHERIA / %d tree)\n",
            length(overlap), nrow(pan), ape::Ntip(tree)))
if (length(overlap) < 500) {
  stop("Fewer than 500 aligned species — alignment quality too low.",
       call. = FALSE)
}

# Prune tree + subset PanTHERIA to overlap
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, overlap))
pan  <- pan[pan$species_key %in% overlap, ]
pan  <- pan[match(tree$tip.label, pan$species_key), ]
rownames(pan) <- pan$species_key

# -------------------------------------------------------------------------
# Trait selection (8 canonical mixed types)
# -------------------------------------------------------------------------

trait_cols <- list(
  body_mass_g         = list(src = "X5.1_AdultBodyMass_g",         type = "cont_log"),
  head_body_length_mm = list(src = "X13.1_AdultHeadBodyLen_mm",    type = "cont_log"),
  gestation_d         = list(src = "X9.1_GestationLen_d",          type = "cont_log"),
  litter_size         = list(src = "X15.1_LitterSize",              type = "count"),
  max_longevity_m     = list(src = "X17.1_MaxLongevity_m",         type = "cont_log"),
  diet_breadth        = list(src = "X6.1_DietBreadth",              type = "ordinal"),
  habitat_breadth     = list(src = "X12.1_HabitatBreadth",         type = "ordinal"),
  terrestriality      = list(src = "X12.2_Terrestriality",          type = "categorical")
)

# Verify columns exist; missing ones are dropped with a notice
keep <- names(trait_cols)[vapply(names(trait_cols), function(nm) {
  trait_cols[[nm]]$src %in% names(pan)
}, logical(1))]
missing_cols <- setdiff(names(trait_cols), keep)
if (length(missing_cols)) {
  cat("Missing PanTHERIA columns (dropped):\n")
  for (m in missing_cols) cat("  ", m, "=>", trait_cols[[m]]$src, "\n")
}
trait_cols <- trait_cols[keep]

df <- data.frame(row.names = pan$species_key)
for (nm in names(trait_cols)) {
  tt <- trait_cols[[nm]]
  x <- pan[[tt$src]]
  if (tt$type == "cont_log") {
    x[x <= 0 | !is.finite(x)] <- NA
    df[[nm]] <- log(x)
  } else if (tt$type == "count") {
    x[!is.finite(x)] <- NA
    df[[nm]] <- as.integer(round(x))
  } else if (tt$type == "ordinal") {
    x[!is.finite(x)] <- NA
    lv <- sort(unique(stats::na.omit(x)))
    df[[nm]] <- factor(as.character(x), levels = as.character(lv),
                        ordered = TRUE)
  } else if (tt$type == "categorical") {
    x[!is.finite(x)] <- NA
    lv <- sort(unique(stats::na.omit(x)))
    df[[nm]] <- factor(as.character(x), levels = as.character(lv))
  }
}

# Drop species with EVERYTHING NA (nothing to impute)
nna <- rowSums(!is.na(df))
keep_rows <- nna > 0L
df   <- df[keep_rows, , drop = FALSE]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(df)))
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(df)))
df   <- df[match(tree$tip.label, rownames(df)), , drop = FALSE]

cat(sprintf("\nFinal: n = %d species × %d traits\n", nrow(df), ncol(df)))
cat("Per-trait observed counts:\n")
print(vapply(df, function(v) sum(!is.na(v)), integer(1)))

# -------------------------------------------------------------------------
# MCAR 30% on top of native NA
# -------------------------------------------------------------------------

SEED      <- 2026L
MISS_FRAC <- 0.30
set.seed(SEED)
df_truth <- df
mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                     dimnames = list(rownames(df), names(df)))
for (v in names(df)) {
  obs_idx <- which(!is.na(df[[v]]))
  to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
  mask_test[to_hide, v] <- TRUE
}
df_miss <- df
for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA

cat("\nHeld-out test cells per trait:\n")
print(colSums(mask_test))

# -------------------------------------------------------------------------
# Methods
# -------------------------------------------------------------------------

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(val = val, wall = proc.time()[["elapsed"]] - t0)
}

method_mean <- function(df_miss, df_truth, mask) {
  out <- df_miss
  for (v in names(out)) {
    if (is.factor(out[[v]])) {
      mode <- names(sort(table(out[[v]], useNA = "no"), decreasing = TRUE))[1]
      out[[v]][mask[, v]] <- factor(mode, levels = levels(out[[v]]),
                                      ordered = is.ordered(out[[v]]))
    } else {
      out[[v]][mask[, v]] <- mean(out[[v]], na.rm = TRUE)
    }
  }
  out
}

method_pigauto <- function(df_miss, tree, em_iter) {
  # n_imputations = 20 captures MC-dropout draws in res$prediction$imputed_datasets
  # so downstream eval can compute both conformal coverage AND MC-dropout
  # quantile coverage on held-out cells. Adds ~15 % wall-clock vs n=1.
  res <- pigauto::impute(df_miss, tree, log_transform = FALSE,
                           missing_frac = 0.20,
                           verbose = FALSE, seed = SEED,
                           epochs = 500L, em_iterations = em_iter,
                           n_imputations = 20L)
  res
}

# -------------------------------------------------------------------------
# Evaluate
# -------------------------------------------------------------------------

safe_cor <- function(x, y) {
  idx <- which(is.finite(x) & is.finite(y))
  if (length(idx) < 2L) return(NA_real_)
  if (stats::sd(x[idx]) == 0 || stats::sd(y[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[idx], y[idx]))
}

# Evaluate a completed data.frame. For pigauto, pass res_obj = the full
# impute() result so conformal AND MC-dropout coverage can also be computed.
# For mean_baseline pass res_obj = NULL (no intervals).
eval_completed <- function(completed, truth, mask, method_name, wall_s,
                            res_obj = NULL) {
  rows <- list()
  # Conformal 95% interval per cell (stored by pigauto's predict path)
  lo <- if (!is.null(res_obj)) res_obj$prediction$conformal_lower else NULL
  hi <- if (!is.null(res_obj)) res_obj$prediction$conformal_upper else NULL
  # MC-dropout draws: list of n_imputations data.frames; per cell take
  # the [q_0.025, q_0.975] quantile of the M draws.
  mi_list <- if (!is.null(res_obj)) res_obj$prediction$imputed_datasets else NULL

  for (v in names(truth)) {
    idx <- which(mask[, v])
    if (!length(idx)) next
    t_v <- truth[[v]][idx]; c_v <- completed[[v]][idx]
    if (is.factor(t_v)) {
      acc <- mean(as.character(c_v) == as.character(t_v), na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method_name, trait = v, metric = "accuracy",
        value = acc, n_cells = length(idx), wall_s = wall_s
      )
    } else {
      rmse <- sqrt(mean((as.numeric(t_v) - as.numeric(c_v))^2, na.rm = TRUE))
      pear <- safe_cor(as.numeric(t_v), as.numeric(c_v))
      rows[[length(rows) + 1L]] <- data.frame(
        method = method_name, trait = v, metric = "rmse",
        value = rmse, n_cells = length(idx), wall_s = wall_s
      )
      rows[[length(rows) + 1L]] <- data.frame(
        method = method_name, trait = v, metric = "pearson_r",
        value = pear, n_cells = length(idx), wall_s = wall_s
      )
      t_num <- as.numeric(t_v)
      # ---- conformal coverage ---------------------------------------
      if (!is.null(lo) && !is.null(hi) && v %in% colnames(lo)) {
        lo_v <- lo[idx, v]; hi_v <- hi[idx, v]
        valid <- is.finite(lo_v) & is.finite(hi_v) & is.finite(t_num)
        if (any(valid)) {
          hits <- t_num[valid] >= lo_v[valid] & t_num[valid] <= hi_v[valid]
          rows[[length(rows) + 1L]] <- data.frame(
            method = method_name, trait = v, metric = "coverage95_conformal",
            value = mean(hits), n_cells = sum(valid), wall_s = wall_s
          )
        }
      }
      # ---- MC-dropout coverage --------------------------------------
      if (!is.null(mi_list) && length(mi_list) > 1L && v %in% names(mi_list[[1]])) {
        # draws_mat: M x length(idx); row m = draw-m values at masked cells
        draws_mat <- vapply(mi_list, function(d) {
          as.numeric(d[idx, v])
        }, numeric(length(idx)))
        if (!is.matrix(draws_mat)) {
          draws_mat <- matrix(draws_mat, ncol = length(mi_list))
        }
        # draws_mat is length(idx) x M after vapply above; check dims
        if (nrow(draws_mat) == length(idx) && ncol(draws_mat) > 1L) {
          q_lo <- apply(draws_mat, 1L, stats::quantile, probs = 0.025,
                         na.rm = TRUE)
          q_hi <- apply(draws_mat, 1L, stats::quantile, probs = 0.975,
                         na.rm = TRUE)
          valid <- is.finite(q_lo) & is.finite(q_hi) & is.finite(t_num)
          if (any(valid)) {
            hits <- t_num[valid] >= q_lo[valid] & t_num[valid] <= q_hi[valid]
            rows[[length(rows) + 1L]] <- data.frame(
              method = method_name, trait = v, metric = "coverage95_mcdropout",
              value = mean(hits), n_cells = sum(valid), wall_s = wall_s
            )
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

cat("\n=== mean_baseline ===\n")
r_mean <- timed(method_mean(df_miss, df_truth, mask_test))
ev_mean <- eval_completed(r_mean$val, df_truth, mask_test,
                            "mean_baseline", r_mean$wall)

cat("=== pigauto_default ===\n")
r_def <- timed(method_pigauto(df_miss, tree, em_iter = 0L))
ev_def <- eval_completed(r_def$val$completed, df_truth, mask_test,
                           "pigauto_default", r_def$wall,
                           res_obj = r_def$val)

cat("=== pigauto_em5 ===\n")
r_em5 <- timed(method_pigauto(df_miss, tree, em_iter = 5L))
ev_em5 <- eval_completed(r_em5$val$completed, df_truth, mask_test,
                           "pigauto_em5", r_em5$wall,
                           res_obj = r_em5$val)

all_rows <- rbind(ev_mean, ev_def, ev_em5)
saveRDS(list(results = all_rows, n_species = nrow(df), n_traits = ncol(df),
              seed = SEED, miss_frac = MISS_FRAC,
              overlap = length(overlap)), out_rds)

md <- c(
  "# PanTHERIA full-scale bench",
  "",
  sprintf("n_species = %d (PanTHERIA ∩ mammal tree), n_traits = %d, seed = %d, miss_frac = %.2f",
          nrow(df), ncol(df), SEED, MISS_FRAC),
  "",
  "## Per-trait metrics",
  "",
  "```",
  capture.output(print(all_rows, row.names = FALSE, max = 500)),
  "```"
)
writeLines(md, out_md)
cat("\n=== DONE ===  ", out_rds, " ", out_md, "\n")
