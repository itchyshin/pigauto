# script/bench_avonet_phase6.R
# Minimal AVONET benchmark post-Phase-6 to compare against BACE's OVR numbers.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-6",
    quiet = TRUE
  )
})

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-6"
out_md <- file.path(here, "script", "bench_avonet_phase6.md")

# Load AVONET 300 (bundled)
data(avonet300, tree300, package = "pigauto")
traits <- avonet300
rownames(traits) <- traits$Species_Key
traits$Species_Key <- NULL
tree <- tree300

set.seed(42)
pd     <- preprocess_traits(traits, tree)
splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                               seed = 42, trait_map = pd$trait_map)

# Phase 6 path
bl_p6 <- tryCatch(fit_baseline(pd, tree, splits = splits),
                   error = function(e) { cat("P6 failed:", conditionMessage(e), "\n"); NULL })

# LP path
orig <- pigauto:::joint_mvn_available
assignInNamespace("joint_mvn_available", function() FALSE, ns = "pigauto")
bl_lp <- fit_baseline(pd, tree, splits = splits)
assignInNamespace("joint_mvn_available", orig, ns = "pigauto")

if (is.null(bl_p6)) {
  writeLines(c("# AVONET Phase 6 bench", "",
                "Phase 6 fit failed on AVONET — see log"), out_md)
  quit(save = "no")
}

# Evaluate on test set
n <- nrow(pd$X_scaled); p <- ncol(pd$X_scaled)
row_i <- ((splits$test_idx - 1L) %% n) + 1L
col_j <- ((splits$test_idx - 1L) %/% n) + 1L

# Metrics per trait type
evaluate_col <- function(bl, tm, col_idx) {
  test_mask <- col_j == col_idx
  if (sum(test_mask) == 0L) return(NULL)
  rows <- row_i[test_mask]
  truth <- pd$X_scaled[rows, col_idx]
  pred  <- bl$mu[rows, col_idx]
  if (tm$type %in% c("continuous", "count", "ordinal", "proportion")) {
    data.frame(metric = "RMSE", value = sqrt(mean((truth - pred)^2)))
  } else if (tm$type == "binary") {
    p <- plogis(pred)
    data.frame(metric = "accuracy",
                value = mean(as.integer(p > 0.5) == truth))
  } else {
    NULL   # categorical handled per-trait, not per-col
  }
}

# Scalar-output traits (continuous + binary)
rows_list <- list()
for (tm in pd$trait_map) {
  if (tm$type == "categorical") next
  for (lc in tm$latent_cols) {
    m_p6 <- evaluate_col(bl_p6, tm, lc)
    m_lp <- evaluate_col(bl_lp, tm, lc)
    if (!is.null(m_p6) && !is.null(m_lp)) {
      rows_list[[length(rows_list) + 1L]] <- data.frame(
        trait = colnames(pd$X_scaled)[lc],
        type  = tm$type,
        metric = m_p6$metric,
        levelC = m_p6$value,
        lp     = m_lp$value
      )
    }
  }
}

# Categorical traits: argmax accuracy across held-out rows
for (tm in pd$trait_map) {
  if (tm$type != "categorical") next
  k_cols   <- tm$latent_cols
  z_rows   <- unique(row_i[col_j %in% k_cols])
  if (length(z_rows) == 0L) next
  z_cols_nm <- colnames(pd$X_scaled)[k_cols]
  truth <- apply(pd$X_scaled[z_rows, z_cols_nm, drop = FALSE], 1, which.max)
  acc_p6 <- mean(apply(bl_p6$mu[z_rows, z_cols_nm, drop = FALSE], 1, which.max) == truth)
  acc_lp <- mean(apply(bl_lp$mu[z_rows, z_cols_nm, drop = FALSE], 1, which.max) == truth)
  cat_name <- sub("=.*$", "", z_cols_nm[1])
  rows_list[[length(rows_list) + 1L]] <- data.frame(
    trait = cat_name, type = "categorical", metric = "accuracy",
    levelC = acc_p6, lp = acc_lp
  )
}

res <- do.call(rbind, rows_list)
res$lift <- ifelse(res$metric == "RMSE", res$lp - res$levelC,
                   res$levelC - res$lp)

md <- c("# AVONET 300 post-Phase-6 benchmark",
        "",
        sprintf("Run on: %s", format(Sys.time())),
        "Level-C baseline (Phase 6 active) vs LP baseline.",
        "",
        "RMSE: lift = lp - levelC (positive = levelC better).",
        "Accuracy: lift = levelC - lp (positive = levelC better).",
        "",
        "```",
        capture.output(print(res, row.names = FALSE)),
        "```")

writeLines(md, out_md)
cat(paste(md, collapse = "\n"), "\n")
