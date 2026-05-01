#!/usr/bin/env Rscript
#
# script/bench_pantheria_bace_head_to_head.R
#
# PanTHERIA head-to-head: pigauto vs BACE on a 500-species random subset
# of the aligned PanTHERIA + mammal tree. Same splits, same seed, same
# metrics. BACE is optional (requireNamespace guard).
#
# Output
#   script/bench_pantheria_bace_head_to_head.{rds,md}

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_pantheria_bace_head_to_head.rds")
out_md  <- file.path(here, "script", "bench_pantheria_bace_head_to_head.md")

cache_dir <- file.path(here, "script", "data-cache")
pantheria_local <- file.path(cache_dir, "pantheria.txt")
tree_local      <- file.path(cache_dir, "mammal_tree.tre")

if (!file.exists(pantheria_local) || !file.exists(tree_local)) {
  stop("Data cache missing. Run script/fetch_pantheria_and_tree.R first.",
       call. = FALSE)
}

SEED      <- 2026L
MISS_FRAC <- 0.30
N_SUBSET  <- 500L

# -------------------------------------------------------------------------
# Same loader / trait selection as bench_pantheria_full.R — duplicated
# inline to keep this bench self-contained.
# -------------------------------------------------------------------------

pan <- utils::read.table(pantheria_local, header = TRUE, sep = "\t",
                          na.strings = "-999",
                          stringsAsFactors = FALSE, quote = "",
                          comment.char = "")
pan$species_key <- paste(pan$MSW93_Genus, pan$MSW93_Species, sep = "_")
pan$species_key <- gsub("[^A-Za-z0-9_]", "", pan$species_key)

tree <- ape::read.tree(tree_local)
tree$tip.label <- gsub("[^A-Za-z0-9_]", "", tree$tip.label)
overlap <- intersect(tree$tip.label, pan$species_key)
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, overlap))
pan  <- pan[pan$species_key %in% overlap, ]
pan  <- pan[match(tree$tip.label, pan$species_key), ]
rownames(pan) <- pan$species_key

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
keep <- names(trait_cols)[vapply(names(trait_cols), function(nm) {
  trait_cols[[nm]]$src %in% names(pan)
}, logical(1))]
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

# Drop species with 0 observed traits
keep_rows <- rowSums(!is.na(df)) > 0L
df   <- df[keep_rows, , drop = FALSE]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(df)))
df   <- df[match(tree$tip.label, rownames(df)), , drop = FALSE]

# -------------------------------------------------------------------------
# Stratified random subset
# -------------------------------------------------------------------------

set.seed(SEED)
n_full <- nrow(df)
sub_idx <- sample.int(n_full, min(N_SUBSET, n_full))
df_sub   <- df[sub_idx, , drop = FALSE]
tree_sub <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(df_sub)))
df_sub   <- df_sub[match(tree_sub$tip.label, rownames(df_sub)), , drop = FALSE]
cat(sprintf("Subset: n = %d (of %d)\n", nrow(df_sub), n_full))

# -------------------------------------------------------------------------
# MCAR 30% test hold-out
# -------------------------------------------------------------------------

set.seed(SEED)
mask_test <- matrix(FALSE, nrow = nrow(df_sub), ncol = ncol(df_sub),
                     dimnames = list(rownames(df_sub), names(df_sub)))
for (v in names(df_sub)) {
  obs_idx <- which(!is.na(df_sub[[v]]))
  to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
  mask_test[to_hide, v] <- TRUE
}
df_miss <- df_sub
for (v in names(df_sub)) df_miss[[v]][mask_test[, v]] <- NA

# -------------------------------------------------------------------------
# Methods
# -------------------------------------------------------------------------

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(val = val, wall = proc.time()[["elapsed"]] - t0)
}

safe_cor <- function(x, y) {
  idx <- which(is.finite(x) & is.finite(y))
  if (length(idx) < 2L) return(NA_real_)
  if (stats::sd(x[idx]) == 0 || stats::sd(y[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[idx], y[idx]))
}

run_pigauto <- function(em_iter) {
  # n_imputations = 20 gives MC-dropout draws for the second coverage type
  pigauto::impute(df_miss, tree_sub, log_transform = FALSE,
                    missing_frac = 0.20, verbose = FALSE, seed = SEED,
                    epochs = 300L, em_iterations = em_iter,
                    n_imputations = 20L)
}

run_bace <- function() {
  if (!requireNamespace("BACE", quietly = TRUE)) return(NULL)
  tryCatch({
    res <- BACE::bace(data = df_miss, tree = tree_sub,
                        n_iter = 2000L, burnin = 500L, thin = 5L,
                        ovr = TRUE, verbose = FALSE)
    # BACE return shape varies; try common slots
    if ("completed" %in% names(res))     res$completed
    else if ("data"  %in% names(res))    res$data
    else if ("imputed_data" %in% names(res)) res$imputed_data
    else NULL
  }, error = function(e) { message("BACE run failed: ", conditionMessage(e)); NULL })
}

eval_completed <- function(completed, truth, mask, method_name, wall_s,
                            res_obj = NULL) {
  if (is.null(completed)) return(NULL)
  rows <- list()
  lo <- if (!is.null(res_obj)) res_obj$prediction$conformal_lower else NULL
  hi <- if (!is.null(res_obj)) res_obj$prediction$conformal_upper else NULL
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
      # Conformal coverage
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
      # MC-dropout coverage
      if (!is.null(mi_list) && length(mi_list) > 1L && v %in% names(mi_list[[1]])) {
        draws_mat <- vapply(mi_list, function(d) as.numeric(d[idx, v]),
                             numeric(length(idx)))
        if (!is.matrix(draws_mat)) {
          draws_mat <- matrix(draws_mat, ncol = length(mi_list))
        }
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
  if (!length(rows)) NULL else do.call(rbind, rows)
}

cat("\n=== pigauto_default ===\n")
r_def <- timed(run_pigauto(0L))
ev_def <- eval_completed(r_def$val$completed, df_sub, mask_test,
                           "pigauto_default", r_def$wall,
                           res_obj = r_def$val)

cat("=== pigauto_em5 ===\n")
r_em5 <- timed(run_pigauto(5L))
ev_em5 <- eval_completed(r_em5$val$completed, df_sub, mask_test,
                           "pigauto_em5", r_em5$wall,
                           res_obj = r_em5$val)

cat("=== BACE (may skip) ===\n")
r_bace <- timed(run_bace())
# BACE doesn't expose conformal intervals; pass NULL for res_obj
ev_bace <- eval_completed(r_bace$val, df_sub, mask_test,
                            "bace_default", r_bace$wall, res_obj = NULL)

all_rows <- do.call(rbind, Filter(Negate(is.null), list(ev_def, ev_em5, ev_bace)))

saveRDS(list(results = all_rows,
              seed = SEED, miss_frac = MISS_FRAC,
              n_subset = nrow(df_sub),
              bace_ran = !is.null(ev_bace)),
         out_rds)

md <- c(
  "# PanTHERIA head-to-head (pigauto vs BACE)",
  "",
  sprintf("Subset: n = %d, seed = %d, miss_frac = %.2f",
          nrow(df_sub), SEED, MISS_FRAC),
  if (!is.null(ev_bace)) "" else "**BACE skipped** (not installed or failed).",
  "",
  "## Per-trait metrics",
  "",
  "```",
  capture.output(print(all_rows, row.names = FALSE, max = 500)),
  "```"
)
writeLines(md, out_md)
cat("\n=== DONE ===  ", out_rds, " ", out_md, "\n")
