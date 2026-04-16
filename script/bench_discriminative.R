#!/usr/bin/env Rscript
# script/bench_discriminative.R
#
# Phase 8 discriminative benchmark: mixed-type, mixed-signal scenarios
# designed to show where the Level-C joint baseline wins and where it ties.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-6",
    quiet = TRUE
  )
})

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-6"
out_md <- file.path(here, "script", "bench_discriminative.md")

n_species <- 200L
reps      <- 3L

# Scenarios: each maps an (evo_model, signal, correlation) tuple to a label
scenarios <- list(
  "high_signal_correlated" = list(lambda = 1.0, rho = 0.7,
    desc = "Strong phylo signal + strong cross-trait correlation"),
  "moderate_signal_correlated" = list(lambda = 0.6, rho = 0.6,
    desc = "Moderate phylo signal + moderate cross-trait correlation"),
  "low_signal_correlated" = list(lambda = 0.2, rho = 0.6,
    desc = "Weak phylo signal + moderate cross-trait correlation"),
  "high_signal_independent" = list(lambda = 1.0, rho = 0.0,
    desc = "Strong phylo signal, no cross-trait correlation")
)

# Build a tree with adjustable Pagel's lambda (scale internal branches)
tree_with_lambda <- function(tree, lambda) {
  if (abs(lambda - 1) < 1e-6) return(tree)
  # Scale internal branches by lambda (tip branches untouched)
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
  # 2 continuous + 1 binary + 1 categorical (K=4)
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
    y = factor(ifelse(z_y > 0, "B", "A"), levels = c("A", "B")),
    z = factor(LETTERS[z_class], levels = LETTERS[1:4]),
    row.names = tree$tip.label,
    stringsAsFactors = FALSE
  )
}

force_lp <- function(fn) {
  orig <- pigauto:::joint_mvn_available
  assignInNamespace("joint_mvn_available", function() FALSE, ns = "pigauto")
  on.exit(assignInNamespace("joint_mvn_available", orig, ns = "pigauto"),
          add = TRUE)
  res <- fn()
  assignInNamespace("joint_mvn_available", orig, ns = "pigauto")
  res
}

run_one <- function(sc_name, sc, rep_id) {
  seed <- rep_id * 1000L + nchar(sc_name)
  set.seed(seed)
  tree0 <- ape::rtree(n_species)
  tree  <- tree_with_lambda(tree0, sc$lambda)
  df    <- simulate_mixed(tree, rho = sc$rho, seed = seed)
  # Impute 30% of EACH trait type
  for (col in c("x1", "x2", "y", "z")) {
    n_hide <- n_species %/% 3
    hide_rows <- sample(n_species, n_hide)
    df[hide_rows, col] <- NA
  }
  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                 seed = rep_id, trait_map = pd$trait_map)
  bl_levelC <- tryCatch(fit_baseline(pd, tree, splits = splits),
                         error = function(e) NULL)
  bl_lp     <- force_lp(function() fit_baseline(pd, tree, splits = splits))
  if (is.null(bl_levelC)) return(NULL)

  n <- nrow(pd$X_scaled); p <- ncol(pd$X_scaled)
  hold_row_i <- ((splits$test_idx - 1L) %% n) + 1L
  hold_col_j <- ((splits$test_idx - 1L) %/% n) + 1L

  # Continuous RMSE (x1, x2)
  cont_cols <- which(colnames(pd$X_scaled) %in% c("x1", "x2"))
  cont_mask <- hold_col_j %in% cont_cols
  if (any(cont_mask)) {
    truth <- pd$X_scaled[cbind(hold_row_i[cont_mask], hold_col_j[cont_mask])]
    p_lc  <- bl_levelC$mu[cbind(hold_row_i[cont_mask], hold_col_j[cont_mask])]
    p_lp  <- bl_lp$mu[cbind(hold_row_i[cont_mask], hold_col_j[cont_mask])]
    rmse_lc <- sqrt(mean((truth - p_lc)^2))
    rmse_lp <- sqrt(mean((truth - p_lp)^2))
  } else {
    rmse_lc <- NA; rmse_lp <- NA
  }

  # Binary accuracy (y)
  y_col <- which(colnames(pd$X_scaled) == "y")
  y_rows <- hold_row_i[hold_col_j == y_col]
  if (length(y_rows) > 0L) {
    truth_y <- pd$X_scaled[y_rows, y_col]
    p_y_lc  <- as.integer(plogis(bl_levelC$mu[y_rows, y_col]) > 0.5)
    p_y_lp  <- as.integer(plogis(bl_lp$mu[y_rows, y_col])    > 0.5)
    acc_y_lc <- mean(p_y_lc == truth_y)
    acc_y_lp <- mean(p_y_lp == truth_y)
  } else {
    acc_y_lc <- NA; acc_y_lp <- NA
  }

  # Categorical accuracy (z)
  z_cols <- grep("^z=", colnames(pd$X_scaled), value = TRUE)
  z_col_idx <- which(colnames(pd$X_scaled) %in% z_cols)
  z_rows <- unique(hold_row_i[hold_col_j %in% z_col_idx])
  if (length(z_rows) > 0L) {
    truth_z <- apply(pd$X_scaled[z_rows, z_cols, drop = FALSE], 1, which.max)
    p_z_lc  <- apply(bl_levelC$mu[z_rows, z_cols, drop = FALSE], 1, which.max)
    p_z_lp  <- apply(bl_lp$mu[z_rows, z_cols, drop = FALSE], 1, which.max)
    acc_z_lc <- mean(p_z_lc == truth_z)
    acc_z_lp <- mean(p_z_lp == truth_z)
  } else {
    acc_z_lc <- NA; acc_z_lp <- NA
  }

  data.frame(
    scenario = sc_name, rep = rep_id,
    rmse_lc = rmse_lc, rmse_lp = rmse_lp, rmse_lift = rmse_lp - rmse_lc,
    acc_y_lc = acc_y_lc, acc_y_lp = acc_y_lp, acc_y_lift = acc_y_lc - acc_y_lp,
    acc_z_lc = acc_z_lc, acc_z_lp = acc_z_lp, acc_z_lift = acc_z_lc - acc_z_lp
  )
}

all_results <- list()
for (sc_name in names(scenarios)) {
  for (r in seq_len(reps)) {
    cat(sprintf("Running %s rep %d ...\n", sc_name, r))
    res <- run_one(sc_name, scenarios[[sc_name]], r)
    if (!is.null(res)) all_results[[length(all_results) + 1L]] <- res
  }
}
all_results <- do.call(rbind, all_results)

agg <- aggregate(cbind(rmse_lc, rmse_lp, rmse_lift,
                         acc_y_lc, acc_y_lp, acc_y_lift,
                         acc_z_lc, acc_z_lp, acc_z_lift) ~ scenario,
                  data = all_results, FUN = mean)

md <- c("# Phase 8 discriminative benchmark",
        "",
        sprintf("Run on: %s", format(Sys.time())),
        sprintf("Species per scenario: %d, reps per scenario: %d", n_species, reps),
        "",
        "Level-C baseline (Phase 6 active) vs LP (force_lp).",
        "",
        "Continuous RMSE (lower better; lift = LP - LevelC):",
        "",
        "```",
        capture.output(print(agg[, c("scenario", "rmse_lc", "rmse_lp", "rmse_lift")],
                              row.names = FALSE)),
        "```",
        "",
        "Binary accuracy (higher better; lift = LevelC - LP):",
        "",
        "```",
        capture.output(print(agg[, c("scenario", "acc_y_lc", "acc_y_lp", "acc_y_lift")],
                              row.names = FALSE)),
        "```",
        "",
        "Categorical (K=4) accuracy (higher better; lift = LevelC - LP):",
        "",
        "```",
        capture.output(print(agg[, c("scenario", "acc_z_lc", "acc_z_lp", "acc_z_lift")],
                              row.names = FALSE)),
        "```",
        "")

writeLines(md, out_md)
cat(paste(md, collapse = "\n"), "\n")
