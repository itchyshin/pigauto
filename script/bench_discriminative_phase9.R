#!/usr/bin/env Rscript
# script/bench_discriminative_phase9.R
# Phase 9 A/B: LP vs Level-C legacy GNN vs Level-C transformer GNN

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/phase-9-transformer",
    quiet = TRUE
  )
})

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/phase-9-transformer"
out_md <- file.path(here, "script", "bench_discriminative_phase9.md")

n_species <- 200L
reps      <- 2L
epochs    <- 100L

scenarios <- list(
  high_signal_correlated = list(lambda = 1.0, rho = 0.7),
  mod_signal_correlated  = list(lambda = 0.6, rho = 0.6),
  low_signal_correlated  = list(lambda = 0.2, rho = 0.6),
  high_signal_indep      = list(lambda = 1.0, rho = 0.0)
)

tree_with_lambda <- function(tree, lambda) {
  if (abs(lambda - 1) < 1e-6) return(tree)
  t2 <- tree
  is_internal <- t2$edge[, 2] > length(t2$tip.label)
  t2$edge.length[is_internal] <- t2$edge.length[is_internal] * lambda
  t2
}

simulate_mixed <- function(tree, rho, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  V <- ape::vcv(tree)
  L <- chol(V)
  z_x1 <- as.numeric(t(L) %*% rnorm(n))
  z_x2 <- rho * z_x1 + sqrt(1 - rho^2) * as.numeric(t(L) %*% rnorm(n))
  z_y  <- rho * z_x1 + sqrt(1 - rho^2) * as.numeric(t(L) %*% rnorm(n))
  z_liab <- matrix(0, nrow = n, ncol = 4)
  for (k in seq_len(4)) {
    z_k <- as.numeric(t(L) %*% rnorm(n))
    z_liab[, k] <- rho * z_x1 + sqrt(1 - rho^2) * z_k
  }
  z_class <- apply(z_liab, 1, which.max)
  data.frame(
    x1 = z_x1, x2 = z_x2,
    y  = factor(ifelse(z_y > 0, "B", "A"), levels = c("A", "B")),
    z  = factor(LETTERS[z_class], levels = LETTERS[1:4]),
    row.names = tree$tip.label,
    stringsAsFactors = FALSE
  )
}

# Temporarily force joint_mvn_available() -> FALSE so fit_baseline uses LP
with_lp_forced <- function(fn) {
  orig <- get("joint_mvn_available", envir = asNamespace("pigauto"))
  assignInNamespace("joint_mvn_available", function() FALSE, ns = "pigauto")
  on.exit(
    assignInNamespace("joint_mvn_available", orig, ns = "pigauto"),
    add = TRUE
  )
  fn()
}

evaluate_one_method <- function(pd, tree, splits, graph, baseline,
                                use_transformer, epochs, seed) {
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    fit_pigauto(
      data = pd, tree = tree, splits = splits, graph = graph,
      baseline = baseline,
      use_transformer_blocks = use_transformer,
      epochs    = as.integer(epochs),
      eval_every = 10L,
      verbose   = FALSE,
      seed      = as.integer(seed)
    ),
    error = function(e) {
      cat("fit failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  wall <- proc.time()[["elapsed"]] - t0
  if (is.null(fit)) return(list(pred = NULL, wall = wall))
  pred <- tryCatch(predict(fit, return_se = FALSE), error = function(e) NULL)
  list(pred = pred, wall = wall)
}

score_pred <- function(pred, pd, splits) {
  if (is.null(pred) || is.null(pred$imputed_latent)) {
    return(data.frame(rmse_cont = NA_real_, acc_bin = NA_real_,
                      acc_cat = NA_real_))
  }
  n <- nrow(pd$X_scaled)
  row_i <- ((splits$test_idx - 1L) %% n) + 1L
  col_j <- ((splits$test_idx - 1L) %/% n) + 1L
  imp <- as.matrix(pred$imputed_latent)

  cont_cols <- which(colnames(pd$X_scaled) %in% c("x1", "x2"))
  cont_mask <- col_j %in% cont_cols
  rmse_c <- if (any(cont_mask)) {
    sqrt(mean((pd$X_scaled[cbind(row_i[cont_mask], col_j[cont_mask])] -
                 imp[cbind(row_i[cont_mask], col_j[cont_mask])])^2))
  } else NA_real_

  y_col  <- which(colnames(pd$X_scaled) == "y")
  y_rows <- row_i[col_j == y_col]
  acc_b <- if (length(y_rows) > 0L) {
    p_hat <- plogis(imp[y_rows, y_col])
    mean(as.integer(p_hat > 0.5) == pd$X_scaled[y_rows, y_col])
  } else NA_real_

  z_cols    <- grep("^z=", colnames(pd$X_scaled), value = TRUE)
  z_col_idx <- which(colnames(pd$X_scaled) %in% z_cols)
  z_rows    <- unique(row_i[col_j %in% z_col_idx])
  acc_c <- if (length(z_rows) > 0L && length(z_cols) > 0L) {
    truth  <- apply(pd$X_scaled[z_rows, z_cols, drop = FALSE], 1, which.max)
    pred_c <- apply(imp[z_rows, z_cols, drop = FALSE], 1, which.max)
    mean(pred_c == truth)
  } else NA_real_

  data.frame(rmse_cont = rmse_c, acc_bin = acc_b, acc_cat = acc_c)
}

run_one <- function(sc_name, sc, rep_id) {
  seed <- rep_id * 1000L + nchar(sc_name)
  set.seed(seed)
  tree0 <- ape::rtree(n_species)
  tree  <- tree_with_lambda(tree0, sc$lambda)
  df    <- simulate_mixed(tree, rho = sc$rho, seed = seed)
  # Hide 30% of each trait
  for (col in c("x1", "x2", "y", "z")) {
    df[sample(n_species, n_species %/% 3L), col] <- NA
  }
  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                seed = rep_id, trait_map = pd$trait_map)
  graph  <- build_phylo_graph(tree, k_eigen = "auto")

  bl_levelC <- fit_baseline(pd, tree, splits = splits, graph = graph)
  bl_lp     <- with_lp_forced(function()
    fit_baseline(pd, tree, splits = splits, graph = graph))

  # Method 1: LP baseline + legacy GNN
  cat(sprintf("  lp ...\n"))
  m_lp <- with_lp_forced(function()
    evaluate_one_method(pd, tree, splits, graph, bl_lp,
                        use_transformer = FALSE, epochs = epochs, seed = seed))

  # Method 2: Level-C baseline + legacy GNN
  cat(sprintf("  levelC_legacy ...\n"))
  m_leg <- evaluate_one_method(pd, tree, splits, graph, bl_levelC,
                               use_transformer = FALSE, epochs = epochs,
                               seed = seed)

  # Method 3: Level-C baseline + transformer GNN
  cat(sprintf("  levelC_transformer ...\n"))
  m_tr <- evaluate_one_method(pd, tree, splits, graph, bl_levelC,
                              use_transformer = TRUE, epochs = epochs,
                              seed = seed)

  rbind(
    cbind(scenario = sc_name, rep = rep_id, method = "lp",
          score_pred(m_lp$pred, pd, splits), wall = round(m_lp$wall, 1)),
    cbind(scenario = sc_name, rep = rep_id, method = "levelC_legacy",
          score_pred(m_leg$pred, pd, splits), wall = round(m_leg$wall, 1)),
    cbind(scenario = sc_name, rep = rep_id, method = "levelC_transformer",
          score_pred(m_tr$pred, pd, splits), wall = round(m_tr$wall, 1))
  )
}

all_results <- list()
t_total_start <- proc.time()[["elapsed"]]
for (sc_name in names(scenarios)) {
  for (r in seq_len(reps)) {
    cat(sprintf("[%s] running %s rep %d ...\n",
                format(Sys.time(), "%H:%M:%S"), sc_name, r))
    res <- tryCatch(
      run_one(sc_name, scenarios[[sc_name]], r),
      error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(res)) all_results[[length(all_results) + 1L]] <- res
  }
}
total_wall <- round(proc.time()[["elapsed"]] - t_total_start, 0)
all_results <- do.call(rbind, all_results)

# Aggregate per scenario x method
agg <- aggregate(cbind(rmse_cont, acc_bin, acc_cat, wall) ~ scenario + method,
                 data = all_results,
                 FUN = function(v) round(mean(v, na.rm = TRUE), 4))

# Wide pivot: one row per scenario
library_safe_pivot <- function(agg) {
  methods <- c("lp", "levelC_legacy", "levelC_transformer")
  out_list <- list()
  for (sc in unique(agg$scenario)) {
    sub <- agg[agg$scenario == sc, , drop = FALSE]
    row <- list(scenario = sc)
    for (m in methods) {
      s <- sub[sub$method == m, , drop = FALSE]
      if (nrow(s) == 0L) {
        row[[paste0(m, "_rmse")]] <- NA
        row[[paste0(m, "_accb")]] <- NA
        row[[paste0(m, "_accc")]] <- NA
        row[[paste0(m, "_wall")]] <- NA
      } else {
        row[[paste0(m, "_rmse")]] <- s$rmse_cont
        row[[paste0(m, "_accb")]] <- s$acc_bin
        row[[paste0(m, "_accc")]] <- s$acc_cat
        row[[paste0(m, "_wall")]] <- s$wall
      }
    }
    out_list[[length(out_list) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  do.call(rbind, out_list)
}
wide <- library_safe_pivot(agg)

md <- c(
  "# Phase 9 A/B: Level-C baseline + transformer GNN vs legacy GNN",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Species per scenario: %d, reps: %d, epochs: %d, total wall: %d s",
          n_species, reps, epochs, total_wall),
  "",
  "Three methods:",
  "- **lp**: LP baseline + legacy GNN (use_transformer_blocks = FALSE)",
  "- **levelC_legacy**: Phase 6 joint-threshold baseline + legacy GNN",
  "- **levelC_transformer**: Phase 6 joint-threshold baseline + Phase 9 transformer GNN",
  "",
  "Metrics: rmse = continuous RMSE (lower=better), accb = binary accuracy,",
  "accc = categorical accuracy (both higher=better), wall = training seconds.",
  "",
  "Wide summary (mean over reps):",
  "",
  "```",
  capture.output(print(wide, row.names = FALSE)),
  "```",
  "",
  "Raw per-rep results:",
  "",
  "```",
  capture.output(print(all_results, row.names = FALSE)),
  "```"
)

writeLines(md, out_md)
cat("\n--- Summary ---\n")
cat(paste(md, collapse = "\n"), "\n")
