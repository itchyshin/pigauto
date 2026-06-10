#!/usr/bin/env Rscript
#
# script/bench_tabpfn_phylo_features.R
#
# Experimental benchmark branch for the Russell-style TabPFN comparison:
# TabPFN regression on AVONET continuous targets with phylogenetic features.
#
# This is not the retired v0.4 TabPFN baseline. That old path was tree-blind.
# This script tests feature-engineered TabPFN variants:
#   plain          same-row trait features only
#   lappe          plain + Laplacian phylogenetic eigenvectors
#   lappe_nfa      lappe + nearest-training-target aggregate features
#   knn            plain + nearest-training-target aggregate features
#
# The driver reuses pigauto's preprocess_traits(), make_missing_splits(), and
# evaluate_imputation() machinery so TabPFN, BM baseline, and pigauto can be
# compared on the same held-out cells.
#
# Smoke check without TabPFN installed:
#   PIGAUTO_TABPFN_DRY_RUN=true PIGAUTO_TABPFN_SCALES=50 \
#     PIGAUTO_TABPFN_REPS=1 Rscript script/bench_tabpfn_phylo_features.R
#
# Small actual run:
#   PIGAUTO_TABPFN_SCALES=50 PIGAUTO_TABPFN_REPS=1 \
#     PIGAUTO_TABPFN_VARIANTS=lappe_nfa \
#     Rscript script/bench_tabpfn_phylo_features.R
#
# Larger Russell-like run, if a GPU TabPFN environment is available:
#   PIGAUTO_TABPFN_SCALES=50,75,300,2000,9993 PIGAUTO_TABPFN_REPS=3 \
#     PIGAUTO_TABPFN_RUN_PIGAUTO=true \
#     PIGAUTO_TABPFN_PYTHON=/path/to/python \
#     Rscript script/bench_tabpfn_phylo_features.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(".", quiet = TRUE)
})

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ...,
      "\n", sep = "")
  flush.console()
}

env_chr <- function(name, default) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) default else x
}

env_int <- function(name, default) {
  as.integer(env_chr(name, as.character(default)))
}

env_num <- function(name, default) {
  as.numeric(env_chr(name, as.character(default)))
}

env_lgl <- function(name, default) {
  x <- tolower(env_chr(name, if (isTRUE(default)) "true" else "false"))
  x %in% c("1", "true", "t", "yes", "y")
}

env_vec <- function(name, default) {
  x <- env_chr(name, paste(default, collapse = ","))
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_rds <- file.path(repo_root, "script", "bench_tabpfn_phylo_features.rds")
out_md <- file.path(repo_root, "script", "bench_tabpfn_phylo_features.md")
runner_py <- file.path(repo_root, "script", "run_tabpfn_phylo.py")

scales <- as.integer(env_vec("PIGAUTO_TABPFN_SCALES", c("50", "75", "300")))
variants <- env_vec("PIGAUTO_TABPFN_VARIANTS",
                    c("plain", "lappe", "lappe_nfa", "knn"))
n_reps <- env_int("PIGAUTO_TABPFN_REPS", 1L)
seed0 <- env_int("PIGAUTO_TABPFN_SEED", 20260610L)
missing_frac <- env_num("PIGAUTO_TABPFN_MISSING_FRAC", 0.25)
val_frac <- env_num("PIGAUTO_TABPFN_VAL_FRAC", 0.5)
epochs <- env_int("PIGAUTO_TABPFN_EPOCHS", 500L)
nfa_k <- env_int("PIGAUTO_TABPFN_NFA_K", 15L)
dry_run <- env_lgl("PIGAUTO_TABPFN_DRY_RUN", FALSE)
run_pigauto <- env_lgl("PIGAUTO_TABPFN_RUN_PIGAUTO", !dry_run)
python <- env_chr("PIGAUTO_TABPFN_PYTHON", "python3")
device <- env_chr("PIGAUTO_TABPFN_DEVICE", "auto")
target_env <- env_vec("PIGAUTO_TABPFN_TARGETS", character())

valid_variants <- c("plain", "lappe", "lappe_nfa", "knn")
if (!all(variants %in% valid_variants)) {
  stop("Unknown variants: ",
       paste(setdiff(variants, valid_variants), collapse = ", "))
}

metadata <- list(
  scales = scales,
  variants = variants,
  n_reps = n_reps,
  seed0 = seed0,
  missing_frac = missing_frac,
  val_frac = val_frac,
  epochs = epochs,
  nfa_k = nfa_k,
  dry_run = dry_run,
  run_pigauto = run_pigauto,
  python = python,
  device = device,
  commit = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                    error = function(e) "unknown")
)

rmse_vec <- function(truth, pred) {
  sqrt(mean((truth - pred)^2, na.rm = TRUE))
}

cor_vec <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (sum(ok) < 2L) return(NA_real_)
  stats::cor(truth[ok], pred[ok])
}

conformal_q <- function(scores, alpha = 0.05) {
  scores <- sort(scores[is.finite(scores)])
  if (!length(scores)) return(NA_real_)
  idx <- ceiling((length(scores) + 1L) * (1 - alpha))
  scores[min(idx, length(scores))]
}

load_avonet_subset <- function(n, seed) {
  env <- new.env(parent = emptyenv())
  if (n <= 300L) {
    data("avonet300", package = "pigauto", envir = env)
    data("tree300", package = "pigauto", envir = env)
    df <- env$avonet300
    tree <- env$tree300
  } else {
    data("avonet_full", package = "pigauto", envir = env)
    data("tree_full", package = "pigauto", envir = env)
    df <- env$avonet_full
    tree <- env$tree_full
  }

  if (n > nrow(df)) {
    stop("Requested n = ", n, " but AVONET source has only ", nrow(df), " rows.")
  }

  if (n < nrow(df)) {
    set.seed(seed)
    species <- sample(df$Species_Key, n)
    tree <- ape::keep.tip(tree, species)
    df <- df[match(tree$tip.label, df$Species_Key), , drop = FALSE]
  } else {
    df <- df[match(tree$tip.label, df$Species_Key), , drop = FALSE]
  }

  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  list(df = df, tree = tree)
}

target_rows_from_split <- function(splits, target_col, n) {
  split_to_rows <- function(idx) {
    row_i <- ((idx - 1L) %% n) + 1L
    col_j <- ((idx - 1L) %/% n) + 1L
    row_i[col_j == target_col]
  }
  list(
    val = split_to_rows(splits$val_idx),
    test = split_to_rows(splits$test_idx)
  )
}

neighbor_features <- function(y, D, train_rows, k) {
  n <- length(y)
  out <- matrix(NA_real_, nrow = n, ncol = 4L)
  colnames(out) <- c("nfa_mean", "nfa_sd", "nfa_min_dist", "nfa_n")

  train_rows <- train_rows[is.finite(y[train_rows])]
  if (length(train_rows) < 2L) return(out)

  for (i in seq_len(n)) {
    candidates <- setdiff(train_rows, i)
    if (!length(candidates)) next
    d <- D[i, candidates]
    ord <- order(d, decreasing = FALSE, na.last = NA)
    if (!length(ord)) next
    keep <- candidates[ord[seq_len(min(k, length(ord)))]]
    vals <- y[keep]
    out[i, "nfa_mean"] <- mean(vals, na.rm = TRUE)
    out[i, "nfa_sd"] <- if (length(vals) > 1L) stats::sd(vals, na.rm = TRUE) else 0
    out[i, "nfa_min_dist"] <- min(D[i, keep], na.rm = TRUE)
    out[i, "nfa_n"] <- length(vals)
  }

  out
}

feature_matrix <- function(pd, graph, splits, target_col, y, train_rows,
                           variant, nfa_k) {
  X <- pd$X_scaled
  X[c(splits$val_idx, splits$test_idx)] <- NA_real_
  X <- X[, setdiff(seq_len(ncol(X)), target_col), drop = FALSE]

  parts <- list(traits = X)
  if (variant %in% c("lappe", "lappe_nfa")) {
    coords <- graph$coords
    colnames(coords) <- paste0("lappe_", seq_len(ncol(coords)))
    parts$lappe <- coords
  }
  if (variant %in% c("lappe_nfa", "knn")) {
    parts$nfa <- neighbor_features(y, graph$D, train_rows, nfa_k)
  }

  out <- do.call(cbind, parts)
  storage.mode(out) <- "double"
  out
}

run_tabpfn_predict <- function(train_x, train_y, pred_x, seed) {
  tmpdir <- tempfile("tabpfn_phylo_")
  dir.create(tmpdir, recursive = TRUE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  train_path <- file.path(tmpdir, "train.csv")
  pred_path <- file.path(tmpdir, "predict.csv")
  out_path <- file.path(tmpdir, "pred.csv")
  meta_path <- file.path(tmpdir, "metadata.json")

  train_frame <- data.frame(.target = train_y, train_x, check.names = FALSE)
  pred_frame <- data.frame(pred_x, check.names = FALSE)
  utils::write.csv(train_frame, train_path, row.names = FALSE, na = "")
  utils::write.csv(pred_frame, pred_path, row.names = FALSE, na = "")

  cmd_args <- c(
    shQuote(runner_py),
    "--train", shQuote(train_path),
    "--predict", shQuote(pred_path),
    "--target", ".target",
    "--out", shQuote(out_path),
    "--metadata", shQuote(meta_path),
    "--device", device,
    "--seed", as.character(seed)
  )

  cmd_out <- system2(
    python,
    cmd_args,
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(cmd_out, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(cmd_out, collapse = "\n"))
  }

  utils::read.csv(out_path)$prediction
}

tabpfn_one_target <- function(pd, graph, splits, target, variant, seed,
                              scale_n, rep_id) {
  tm <- pd$trait_map[[target]]
  target_col <- tm$latent_cols[1L]
  n <- nrow(pd$X_scaled)
  y <- pd$X_scaled[, target_col]
  rows <- target_rows_from_split(splits, target_col, n)
  val_rows <- rows$val[is.finite(y[rows$val])]
  test_rows <- rows$test[is.finite(y[rows$test])]
  train_rows <- setdiff(which(is.finite(y)), c(val_rows, test_rows))

  x <- feature_matrix(pd, graph, splits, target_col, y, train_rows,
                      variant, nfa_k)

  base_row <- data.frame(
    method = paste0("tabpfn_", variant),
    scale_n = scale_n,
    rep = rep_id,
    trait = target,
    split = "test",
    n_train = length(train_rows),
    n_val = length(val_rows),
    n_test = length(test_rows),
    n_features = ncol(x),
    rmse = NA_real_,
    pearson_r = NA_real_,
    coverage_95 = NA_real_,
    width_95 = NA_real_,
    qhat = NA_real_,
    wall_sec = NA_real_,
    status = if (dry_run) "dry_run" else "pending",
    stringsAsFactors = FALSE
  )

  if (length(train_rows) < 5L || length(val_rows) < 2L ||
      length(test_rows) < 2L) {
    base_row$status <- "skipped_too_few_cells"
    return(base_row)
  }
  if (dry_run) return(base_row)

  t0 <- proc.time()[["elapsed"]]
  fit_result <- tryCatch({
    pred_rows <- c(val_rows, test_rows)
    pred_all <- run_tabpfn_predict(x[train_rows, , drop = FALSE],
                                   y[train_rows],
                                   x[pred_rows, , drop = FALSE],
                                   seed)
    val_pred <- pred_all[seq_along(val_rows)]
    test_pred <- pred_all[length(val_rows) + seq_along(test_rows)]
    qhat <- conformal_q(abs(y[val_rows] - val_pred), alpha = 0.05)
    list(ok = TRUE, val_pred = val_pred, test_pred = test_pred, qhat = qhat)
  }, error = function(e) {
    list(ok = FALSE, message = conditionMessage(e))
  })

  base_row$wall_sec <- proc.time()[["elapsed"]] - t0
  if (!isTRUE(fit_result$ok)) {
    base_row$status <- paste0("error: ", fit_result$message)
    return(base_row)
  }

  qhat <- fit_result$qhat
  base_row$rmse <- rmse_vec(y[test_rows], fit_result$test_pred)
  base_row$pearson_r <- cor_vec(y[test_rows], fit_result$test_pred)
  base_row$qhat <- qhat
  base_row$coverage_95 <- if (is.finite(qhat)) {
    mean(abs(y[test_rows] - fit_result$test_pred) <= qhat)
  } else {
    NA_real_
  }
  base_row$width_95 <- 2 * qhat
  base_row$status <- "ok"
  base_row
}

standardize_eval <- function(ev, method, scale_n, rep_id, targets, wall_sec) {
  if (is.null(ev) || !nrow(ev)) return(NULL)
  keep <- ev$split == "test" & ev$trait %in% targets
  ev <- ev[keep, , drop = FALSE]
  if (!nrow(ev)) return(NULL)
  data.frame(
    method = method,
    scale_n = scale_n,
    rep = rep_id,
    trait = ev$trait,
    split = ev$split,
    n_train = NA_integer_,
    n_val = NA_integer_,
    n_test = ev$n,
    n_features = NA_integer_,
    rmse = ev$rmse,
    pearson_r = ev$pearson_r,
    coverage_95 = ev$coverage_95,
    width_95 = NA_real_,
    qhat = NA_real_,
    wall_sec = wall_sec,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

pigauto_methods <- function(pd, tree, graph, splits, targets, seed,
                            scale_n, rep_id) {
  if (!run_pigauto) return(NULL)

  t0 <- proc.time()[["elapsed"]]
  baseline <- fit_baseline(pd, tree, splits = splits, graph = graph)
  ev_bl <- evaluate_imputation(baseline$mu, pd$X_scaled, splits,
                               trait_map = pd$trait_map)
  baseline_wall <- proc.time()[["elapsed"]] - t0

  graph_fit <- graph
  graph_fit$D <- NULL
  invisible(gc(full = TRUE, verbose = FALSE))

  t1 <- proc.time()[["elapsed"]]
  fit <- fit_pigauto(pd, tree, splits = splits, baseline = baseline,
                     graph = graph_fit, epochs = epochs, verbose = FALSE,
                     seed = seed)
  pred <- stats::predict(fit, return_se = TRUE, n_imputations = 1L)
  ev_pg <- evaluate_imputation(pred, pd$X_scaled, splits)
  pigauto_wall <- proc.time()[["elapsed"]] - t1

  rbind(
    standardize_eval(ev_bl, "baseline", scale_n, rep_id, targets,
                     baseline_wall),
    standardize_eval(ev_pg, "pigauto", scale_n, rep_id, targets,
                     pigauto_wall)
  )
}

write_summary <- function(results, metadata, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("# TabPFN Phylo-Feature Benchmark\n", con)
  writeLines(sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("- Git commit: `%s`", metadata$commit), con)
  writeLines(sprintf("- Scales: `%s`", paste(metadata$scales, collapse = ", ")), con)
  writeLines(sprintf("- Variants: `%s`", paste(metadata$variants, collapse = ", ")), con)
  writeLines(sprintf("- Replicates: `%d`", metadata$n_reps), con)
  writeLines(sprintf("- Dry run: `%s`", metadata$dry_run), con)
  writeLines("", con)

  writeLines("## Design", con)
  writeLines("", con)
  writeLines(paste(
    "Continuous AVONET targets are held out with pigauto's",
    "`make_missing_splits()` and scored on the same test cells used for",
    "the BM baseline and pigauto. TabPFN feature sets are built from same-row",
    "traits, optional Laplacian phylogenetic eigenvectors, and optional",
    "nearest-training-target aggregate features from cophenetic distances.",
    "Split conformal intervals use validation residuals from the same split."
  ), con)
  writeLines("", con)

  writeLines("## Status Counts", con)
  writeLines("", con)
  status_counts <- as.data.frame(table(
    method = results$method,
    status = results$status,
    useNA = "ifany"
  ))
  status_counts <- status_counts[status_counts$Freq > 0L, , drop = FALSE]
  writeLines(capture.output(print(status_counts, row.names = FALSE)), con)
  writeLines("", con)

  scored <- results[results$status == "ok", , drop = FALSE]
  if (nrow(scored)) {
    writeLines("## Test Summary", con)
    writeLines("", con)
    agg <- stats::aggregate(
      cbind(rmse, pearson_r, coverage_95) ~ method + scale_n + trait,
      data = scored,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    agg <- agg[order(agg$scale_n, agg$trait, agg$method), , drop = FALSE]
    writeLines(capture.output(print(agg, row.names = FALSE, digits = 4)), con)
  } else {
    writeLines("## Test Summary", con)
    writeLines("", con)
    if (isTRUE(metadata$dry_run)) {
      writeLines("No scored rows yet. This is expected for `PIGAUTO_TABPFN_DRY_RUN=true`.", con)
    } else {
      writeLines("No scored rows were produced. Inspect `Status Counts` above for the failing method status.", con)
    }
  }
}

all_results <- list()

for (scale_n in scales) {
  for (rep_id in seq_len(n_reps)) {
    seed <- seed0 + scale_n * 100L + rep_id
    log_line("Preparing AVONET n=", scale_n, " rep=", rep_id)
    dat <- load_avonet_subset(scale_n, seed)
    pd <- preprocess_traits(dat$df, dat$tree, log_transform = TRUE)
    graph <- build_phylo_graph(dat$tree, k_eigen = "auto")
    splits <- make_missing_splits(pd$X_scaled, missing_frac = missing_frac,
                                  val_frac = val_frac, seed = seed,
                                  trait_map = pd$trait_map)

    continuous_targets <- names(Filter(function(tm) {
      identical(tm$type, "continuous")
    }, pd$trait_map))
    if (length(target_env)) {
      continuous_targets <- intersect(continuous_targets, target_env)
    }
    if (!length(continuous_targets)) {
      stop("No continuous targets selected.")
    }

    log_line("Targets: ", paste(continuous_targets, collapse = ", "))

    all_results[[length(all_results) + 1L]] <-
      pigauto_methods(pd, dat$tree, graph, splits, continuous_targets,
                      seed, scale_n, rep_id)

    for (variant in variants) {
      for (target in continuous_targets) {
        log_line("TabPFN ", variant, " target=", target,
                 " n=", scale_n, " rep=", rep_id)
        row <- tryCatch(
          tabpfn_one_target(pd, graph, splits, target, variant, seed,
                            scale_n, rep_id),
          error = function(e) {
            data.frame(
              method = paste0("tabpfn_", variant),
              scale_n = scale_n,
              rep = rep_id,
              trait = target,
              split = "test",
              n_train = NA_integer_,
              n_val = NA_integer_,
              n_test = NA_integer_,
              n_features = NA_integer_,
              rmse = NA_real_,
              pearson_r = NA_real_,
              coverage_95 = NA_real_,
              width_95 = NA_real_,
              qhat = NA_real_,
              wall_sec = NA_real_,
              status = paste0("error: ", conditionMessage(e)),
              stringsAsFactors = FALSE
            )
          }
        )
        all_results[[length(all_results) + 1L]] <- row
      }
    }

    results <- do.call(rbind, all_results)
    saveRDS(list(results = results, metadata = metadata), out_rds)
    write_summary(results, metadata, out_md)
  }
}

results <- do.call(rbind, all_results)
saveRDS(list(results = results, metadata = metadata), out_rds)
write_summary(results, metadata, out_md)
log_line("Wrote ", out_rds)
log_line("Wrote ", out_md)
