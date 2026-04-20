#!/usr/bin/env Rscript
#
# script/bench_clade_missingness.R
#
# Phase 8.3: clade-correlated missingness (realistic MAR pattern). Real
# databases under-sample whole clades at a time, not random species.
# Tests whether pigauto's phylo prior actually helps when the missing
# cells cluster on the tree.
#
# Design
#   - Tree:    ape::rcoal(300).
#   - Traits:  3 continuous + 1 binary + 1 categorical K=3 (mixed type).
#   - Missingness: pick random internal nodes such that the union of
#                  their descendants covers approximately `target_frac`
#                  of tips; mask ONE focal trait in all covered tips.
#                  Non-focal traits masked MCAR at 0.10 for realism.
#   - target_frac: {0.10, 0.25, 0.40} (depth of clade hole).
#   - Reps:    3.
#   - Methods: mean_baseline, pigauto_default, pigauto_em5.
#
# Output
#   script/bench_clade_missingness.{rds,md}
#
# Runtime: ~20 min.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_clade_missingness.rds")
out_md  <- file.path(here, "script", "bench_clade_missingness.md")

CONFIG <- list(
  target_fracs = c(0.10, 0.25, 0.40),
  n_species    = 300L,
  n_reps       = 3L,
  methods      = c("mean_baseline", "pigauto_default", "pigauto_em5"),
  epochs       = 200L,
  focal_trait  = "c1"   # continuous; clade-masked
)

# -------------------------------------------------------------------------
# Mixed-type simulator (3 cont + 1 bin + 1 cat K=3) under strong BM
# -------------------------------------------------------------------------

sim_mixed_bm <- function(tree, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  L <- t(chol(ape::vcv(tree) + diag(1e-8, n)))
  latent <- matrix(stats::rnorm(n * 6), n, 6)
  latent <- L %*% latent  # apply tree-side correlation
  c1 <- as.numeric(latent[, 1])
  c2 <- as.numeric(latent[, 2])
  c3 <- as.numeric(latent[, 3])
  b1 <- factor(ifelse(latent[, 4] > 0, "1", "0"), levels = c("0", "1"))
  cat_scores <- latent[, 5:6]
  cat_scores <- cbind(cat_scores, 0)   # 3 classes, third = 0
  cat_labels <- apply(cat_scores, 1L, function(r) c("A","B","C")[which.max(r)])
  cat1 <- factor(cat_labels, levels = c("A", "B", "C"))
  df <- data.frame(c1 = c1, c2 = c2, c3 = c3, b1 = b1, cat1 = cat1,
                    row.names = tree$tip.label)
  df
}

# -------------------------------------------------------------------------
# Clade-correlated missingness on focal trait
# -------------------------------------------------------------------------
#
# Pick random internal nodes, collect their tip descendants, union until
# the union covers ~target_frac of the tree. Mask the focal trait in all
# covered tips. Non-focal traits get a light MCAR pattern.
pick_clade_tips <- function(tree, target_frac, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  n_int <- tree$Nnode
  internal <- seq(n + 1, n + n_int)
  target_n <- max(1L, min(n - 2L, ceiling(n * target_frac)))
  # Prefer smaller clades — shuffle internal nodes and order by descendant
  # count (ascending) to add small clades first and avoid overshooting the
  # target by a huge margin.
  order_random <- sample(internal)
  sizes <- vapply(order_random, function(nd) {
    length(phangorn::Descendants(tree, nd, type = "tips")[[1]])
  }, integer(1))
  order_random <- order_random[order(sizes)]
  covered <- integer(0)
  for (nd in order_random) {
    desc <- phangorn::Descendants(tree, nd, type = "tips")[[1]]
    covered <- union(covered, desc)
    if (length(covered) >= target_n) break
  }
  # Cap at target_n so we never mask 100% of the tree
  if (length(covered) > target_n) covered <- covered[seq_len(target_n)]
  covered
}

# Fallback if phangorn not available
pick_clade_tips_fallback <- function(tree, target_frac, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  target_n <- max(1L, min(n - 2L, ceiling(n * target_frac)))
  # Use ape::node.depth + descendants: pick random tips and add their
  # nearest-10% neighbourhood via cophenetic distance.
  D <- ape::cophenetic.phylo(tree)
  covered <- integer(0)
  order_random <- sample(seq_len(n))
  thresh <- stats::quantile(D[lower.tri(D)], 0.10)
  for (root_tip in order_random) {
    neighbours <- which(D[root_tip, ] <= thresh)
    covered <- union(covered, neighbours)
    if (length(covered) >= target_n) break
  }
  if (length(covered) > target_n) covered <- covered[seq_len(target_n)]
  covered
}

inject_clade_miss <- function(df, tree, focal, target_frac, seed) {
  tips_masked <- if (requireNamespace("phangorn", quietly = TRUE)) {
    pick_clade_tips(tree, target_frac, seed)
  } else {
    pick_clade_tips_fallback(tree, target_frac, seed)
  }
  tip_names_masked <- tree$tip.label[tips_masked]
  miss_idx <- which(rownames(df) %in% tip_names_masked)

  mask <- vector("list", length(names(df)))
  names(mask) <- names(df)
  for (v in names(df)) mask[[v]] <- integer(0)
  mask[[focal]] <- miss_idx
  df[[focal]][miss_idx] <- NA

  # Light MCAR 0.10 on other traits
  set.seed(seed + 7)
  for (v in setdiff(names(df), focal)) {
    idx <- sample.int(nrow(df), ceiling(nrow(df) * 0.10))
    mask[[v]] <- idx
    df[idx, v] <- NA
  }

  list(data = df, mask = mask, clade_size = length(miss_idx))
}

method_mean <- function(df_miss, tree, mask) {
  out <- df_miss
  for (v in names(out)) {
    if (is.factor(out[[v]])) {
      mode <- names(sort(table(out[[v]], useNA = "no"), decreasing = TRUE))[1]
      out[[v]][mask[[v]]] <- factor(mode, levels = levels(out[[v]]))
    } else {
      out[[v]][mask[[v]]] <- mean(out[[v]], na.rm = TRUE)
    }
  }
  out
}

method_pigauto <- function(df_miss, tree, em_iter, seed) {
  res <- pigauto::impute(df_miss, tree, log_transform = FALSE,
                           missing_frac = 0.25, n_imputations = 1L,
                           epochs = CONFIG$epochs, verbose = FALSE,
                           seed = seed, em_iterations = em_iter)
  res$completed
}

safe_cor <- function(x, y) {
  # Returns NA_real_ if either vector is constant or has < 2 complete pairs.
  idx <- which(is.finite(x) & is.finite(y))
  if (length(idx) < 2L) return(NA_real_)
  if (stats::sd(x[idx]) == 0 || stats::sd(y[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[idx], y[idx]))
}

eval_cell <- function(truth, completed, mask) {
  rows <- list()
  for (v in names(truth)) {
    idx <- mask[[v]]
    if (length(idx) == 0L) next
    t_v <- truth[[v]][idx]; c_v <- completed[[v]][idx]
    if (is.factor(t_v)) {
      acc <- mean(as.character(c_v) == as.character(t_v), na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(trait = v, metric = "accuracy",
                                                value = acc)
    } else {
      rmse <- sqrt(mean((t_v - c_v)^2, na.rm = TRUE))
      pear <- safe_cor(as.numeric(t_v), as.numeric(c_v))
      rows[[length(rows) + 1L]] <- data.frame(trait = v,
        metric = c("rmse", "pearson_r"), value = c(rmse, pear))
    }
  }
  do.call(rbind, rows)
}

cat(sprintf("Phase 8.3 clade-MAR sweep: %d fracs × %d reps × %d methods\n",
            length(CONFIG$target_fracs), CONFIG$n_reps,
            length(CONFIG$methods)))
results <- list()
t_start <- proc.time()[["elapsed"]]

for (rep_i in seq_len(CONFIG$n_reps)) {
  set.seed(rep_i)
  tree <- ape::rcoal(CONFIG$n_species,
                      tip.label = paste0("sp", seq_len(CONFIG$n_species)))
  truth <- sim_mixed_bm(tree, seed = rep_i * 7L)

  for (frac in CONFIG$target_fracs) {
    sim_seed <- 1000L * rep_i + as.integer(frac * 100)
    miss <- inject_clade_miss(truth, tree, CONFIG$focal_trait, frac, sim_seed)
    cat(sprintf("  rep=%d frac=%.2f focal_miss=%d tips\n",
                rep_i, frac, miss$clade_size))

    for (meth in CONFIG$methods) {
      t0 <- proc.time()[["elapsed"]]
      completed <- tryCatch(switch(meth,
        mean_baseline   = method_mean(miss$data, tree, miss$mask),
        pigauto_default = method_pigauto(miss$data, tree, 0L, sim_seed + 2),
        pigauto_em5     = method_pigauto(miss$data, tree, 5L, sim_seed + 2)
      ), error = function(e) { message("  ", meth, " @ frac=", frac, ": ",
                                          conditionMessage(e)); NULL })
      wall <- proc.time()[["elapsed"]] - t0
      if (is.null(completed)) next

      ev <- eval_cell(truth, completed, miss$mask)
      ev$method <- meth; ev$target_frac <- frac
      ev$focal_miss <- miss$clade_size; ev$rep <- rep_i; ev$wall_s <- wall
      results[[length(results) + 1L]] <- ev
      cat(sprintf("    %-18s %.1fs\n", meth, wall))
      saveRDS(list(results = do.call(rbind, results),
                    config = CONFIG,
                    script_wall = proc.time()[["elapsed"]] - t_start), out_rds)
    }
  }
}

all_res <- do.call(rbind, results)
saveRDS(list(results = all_res, config = CONFIG,
              script_wall = proc.time()[["elapsed"]] - t_start), out_rds)

summary_tbl <- aggregate(value ~ method + target_frac + trait + metric,
                          data = all_res, FUN = mean)
md <- c(
  "# Phase 8.3: clade-correlated missingness (realistic MAR)",
  "",
  sprintf("n=%d, mixed types, n_reps=%d, focal trait=%s",
          CONFIG$n_species, CONFIG$n_reps, CONFIG$focal_trait),
  sprintf("Total wall: %.1f min",
          (proc.time()[["elapsed"]] - t_start) / 60),
  "",
  "## Focal-trait metrics per (method, target_frac)",
  "",
  "```",
  capture.output(print(summary_tbl[summary_tbl$trait == CONFIG$focal_trait, ],
                        row.names = FALSE)),
  "```",
  "",
  "## All traits summary",
  "",
  "```",
  capture.output(print(summary_tbl, row.names = FALSE, max = 500)),
  "```"
)
writeLines(md, out_md)
cat("\n=== DONE === ", out_rds, " ", out_md, "\n")
