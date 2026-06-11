#!/usr/bin/env Rscript
#
# script/bench_tabpfn_reconciliation.R
#
# Branch-local reconciliation benchmark for the Russell-style TabPFN result.
# The aim is to separate same-row cross-trait signal from phylogenetic signal.
#
# Regimes:
#   same_row                 other same-row traits only
#   phylo_only               Laplacian phylogenetic eigenvectors only
#   same_row_lappe           other same-row traits + Laplacian eigenvectors
#   shuffled_same_row_lappe  shuffled other-trait rows + Laplacian eigenvectors
#   same_row_lappe_nfa       same_row_lappe + nearest-target phylo aggregates
#
# Smoke:
#   PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_reconciliation_smoke \
#   PIGAUTO_TABPFN_SCALES=50 PIGAUTO_TABPFN_REPS=1 \
#   PIGAUTO_TABPFN_TARGETS=Mass \
#   PIGAUTO_TABPFN_REGIMES=same_row,phylo_only,shuffled_same_row_lappe \
#   PIGAUTO_TABPFN_RUN_PIGAUTO=false \
#     Rscript script/bench_tabpfn_reconciliation.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(".", quiet = TRUE)
})

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(
    sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
    ...,
    "\n",
    sep = ""
  )
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
out_stem <- env_chr("PIGAUTO_TABPFN_OUT_STEM", "bench_tabpfn_reconciliation")
if (!grepl("^[A-Za-z0-9_.-]+$", out_stem)) {
  stop("PIGAUTO_TABPFN_OUT_STEM must be a simple file stem, not a path.")
}
out_rds <- file.path(repo_root, "script", paste0(out_stem, ".rds"))
out_md <- file.path(repo_root, "script", paste0(out_stem, ".md"))
runner_py <- file.path(repo_root, "script", "run_tabpfn_phylo.py")

scales <- as.integer(env_vec("PIGAUTO_TABPFN_SCALES", c("50", "75", "300")))
regimes <- env_vec(
  "PIGAUTO_TABPFN_REGIMES",
  c(
    "same_row",
    "phylo_only",
    "same_row_lappe",
    "shuffled_same_row_lappe",
    "same_row_lappe_nfa"
  )
)
n_reps <- env_int("PIGAUTO_TABPFN_REPS", 1L)
seed0 <- env_int("PIGAUTO_TABPFN_SEED", 20260611L)
missing_frac <- env_num("PIGAUTO_TABPFN_MISSING_FRAC", 0.25)
val_frac <- env_num("PIGAUTO_TABPFN_VAL_FRAC", 0.5)
epochs <- env_int("PIGAUTO_TABPFN_EPOCHS", 500L)
nfa_k <- env_int("PIGAUTO_TABPFN_NFA_K", 15L)
dry_run <- env_lgl("PIGAUTO_TABPFN_DRY_RUN", FALSE)
run_pigauto <- env_lgl("PIGAUTO_TABPFN_RUN_PIGAUTO", !dry_run)
python <- env_chr("PIGAUTO_TABPFN_PYTHON", "python3")
device <- env_chr("PIGAUTO_TABPFN_DEVICE", "auto")
target_env <- env_vec("PIGAUTO_TABPFN_TARGETS", character())

valid_regimes <- c(
  "same_row",
  "phylo_only",
  "same_row_lappe",
  "shuffled_same_row_lappe",
  "same_row_lappe_nfa"
)
if (!all(regimes %in% valid_regimes)) {
  stop(
    "Unknown regimes: ",
    paste(setdiff(regimes, valid_regimes), collapse = ", ")
  )
}

metadata <- list(
  benchmark = "tabpfn_reconciliation",
  scales = scales,
  regimes = regimes,
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
  out_stem = out_stem,
  commit = tryCatch(
    system("git rev-parse HEAD", intern = TRUE),
    error = function(e) "unknown"
  )
)

rmse_vec <- function(truth, pred) {
  sqrt(mean((truth - pred)^2, na.rm = TRUE))
}

cor_vec <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (sum(ok) < 2L) {
    return(NA_real_)
  }
  stats::cor(truth[ok], pred[ok])
}

conformal_q <- function(scores, alpha = 0.05) {
  scores <- sort(scores[is.finite(scores)])
  if (!length(scores)) {
    return(NA_real_)
  }
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
    stop(
      "Requested n = ",
      n,
      " but AVONET source has only ",
      nrow(df),
      " rows."
    )
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
  if (length(train_rows) < 2L) {
    return(out)
  }

  for (i in seq_len(n)) {
    candidates <- setdiff(train_rows, i)
    if (!length(candidates)) {
      next
    }
    d <- D[i, candidates]
    ord <- order(d, decreasing = FALSE, na.last = NA)
    if (!length(ord)) {
      next
    }
    keep <- candidates[ord[seq_len(min(k, length(ord)))]]
    vals <- y[keep]
    out[i, "nfa_mean"] <- mean(vals, na.rm = TRUE)
    out[i, "nfa_sd"] <- if (length(vals) > 1L) {
      stats::sd(vals, na.rm = TRUE)
    } else {
      0
    }
    out[i, "nfa_min_dist"] <- min(D[i, keep], na.rm = TRUE)
    out[i, "nfa_n"] <- length(vals)
  }

  out
}

feature_matrix <- function(
  pd,
  graph,
  splits,
  target_col,
  y,
  train_rows,
  regime,
  nfa_k,
  seed
) {
  X <- pd$X_scaled
  X[c(splits$val_idx, splits$test_idx)] <- NA_real_
  X <- X[, setdiff(seq_len(ncol(X)), target_col), drop = FALSE]
  parts <- list()

  if (regime %in% c("same_row", "same_row_lappe", "same_row_lappe_nfa")) {
    parts$traits <- X
  }

  if (identical(regime, "shuffled_same_row_lappe")) {
    set.seed(seed)
    parts$traits_shuffled <- X[sample(seq_len(nrow(X))), , drop = FALSE]
  }

  if (
    regime %in%
      c(
        "phylo_only",
        "same_row_lappe",
        "same_row_lappe_nfa",
        "shuffled_same_row_lappe"
      )
  ) {
    coords <- graph$coords
    colnames(coords) <- paste0("lappe_", seq_len(ncol(coords)))
    parts$lappe <- coords
  }

  if (identical(regime, "same_row_lappe_nfa")) {
    parts$nfa <- neighbor_features(y, graph$D, train_rows, nfa_k)
  }

  if (!length(parts)) {
    stop("No predictors available for regime ", regime)
  }

  out <- do.call(cbind, parts)
  storage.mode(out) <- "double"
  out
}

run_tabpfn_predict <- function(train_x, train_y, pred_x, seed) {
  tmpdir <- tempfile("tabpfn_reconciliation_")
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
    "--train",
    shQuote(train_path),
    "--predict",
    shQuote(pred_path),
    "--target",
    ".target",
    "--out",
    shQuote(out_path),
    "--metadata",
    shQuote(meta_path),
    "--device",
    device,
    "--seed",
    as.character(seed)
  )
  cmd_out <- system2(python, cmd_args, stdout = TRUE, stderr = TRUE)
  status <- attr(cmd_out, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(cmd_out, collapse = "\n"))
  }

  utils::read.csv(out_path)$prediction
}

empty_row <- function(method, scale_n, rep_id, target, status) {
  data.frame(
    method = method,
    regime = if (startsWith(method, "tabpfn_")) {
      sub("^tabpfn_", "", method)
    } else {
      NA_character_
    },
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
    r_cal_bm = NA_real_,
    r_cal_gnn = NA_real_,
    r_cal_mean = NA_real_,
    wall_sec = NA_real_,
    status = status,
    stringsAsFactors = FALSE
  )
}

tabpfn_one_target <- function(
  pd,
  graph,
  splits,
  target,
  regime,
  seed,
  scale_n,
  rep_id
) {
  tm <- pd$trait_map[[target]]
  target_col <- tm$latent_cols[1L]
  n <- nrow(pd$X_scaled)
  y <- pd$X_scaled[, target_col]
  rows <- target_rows_from_split(splits, target_col, n)
  val_rows <- rows$val[is.finite(y[rows$val])]
  test_rows <- rows$test[is.finite(y[rows$test])]
  train_rows <- setdiff(which(is.finite(y)), c(val_rows, test_rows))

  base_row <- empty_row(
    paste0("tabpfn_", regime),
    scale_n,
    rep_id,
    target,
    if (dry_run) "dry_run" else "pending"
  )
  base_row$n_train <- length(train_rows)
  base_row$n_val <- length(val_rows)
  base_row$n_test <- length(test_rows)

  if (
    length(train_rows) < 5L ||
      length(val_rows) < 2L ||
      length(test_rows) < 2L
  ) {
    base_row$status <- "skipped_too_few_cells"
    return(base_row)
  }

  x <- feature_matrix(
    pd,
    graph,
    splits,
    target_col,
    y,
    train_rows,
    regime,
    nfa_k,
    seed
  )
  base_row$n_features <- ncol(x)
  if (dry_run) {
    return(base_row)
  }

  t0 <- proc.time()[["elapsed"]]
  fit_result <- tryCatch(
    {
      pred_rows <- c(val_rows, test_rows)
      pred_all <- run_tabpfn_predict(
        x[train_rows, , drop = FALSE],
        y[train_rows],
        x[pred_rows, , drop = FALSE],
        seed
      )
      val_pred <- pred_all[seq_along(val_rows)]
      test_pred <- pred_all[length(val_rows) + seq_along(test_rows)]
      qhat <- conformal_q(abs(y[val_rows] - val_pred), alpha = 0.05)
      list(ok = TRUE, test_pred = test_pred, qhat = qhat)
    },
    error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    }
  )
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

gate_map <- function(fit, pd, targets, field) {
  vals <- fit[[field]]
  if (is.null(vals)) {
    return(setNames(rep(NA_real_, length(targets)), targets))
  }
  out <- vapply(
    targets,
    function(target) {
      col <- pd$trait_map[[target]]$latent_cols[1L]
      if (length(vals) < col) {
        NA_real_
      } else {
        as.numeric(vals[[col]])
      }
    },
    numeric(1L)
  )
  out
}

standardize_eval <- function(
  ev,
  method,
  scale_n,
  rep_id,
  targets,
  wall_sec,
  gate_bm = NULL,
  gate_gnn = NULL,
  gate_mean = NULL
) {
  if (is.null(ev) || !nrow(ev)) {
    return(NULL)
  }
  keep <- ev$split == "test" & ev$trait %in% targets
  ev <- ev[keep, , drop = FALSE]
  if (!nrow(ev)) {
    return(NULL)
  }
  gates <- function(x) {
    if (is.null(x)) {
      rep(NA_real_, nrow(ev))
    } else {
      as.numeric(x[ev$trait])
    }
  }
  data.frame(
    method = method,
    regime = NA_character_,
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
    r_cal_bm = gates(gate_bm),
    r_cal_gnn = gates(gate_gnn),
    r_cal_mean = gates(gate_mean),
    wall_sec = wall_sec,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

pigauto_methods <- function(
  pd,
  tree,
  graph,
  splits,
  targets,
  seed,
  scale_n,
  rep_id
) {
  if (!run_pigauto) {
    return(NULL)
  }

  t0 <- proc.time()[["elapsed"]]
  baseline <- fit_baseline(pd, tree, splits = splits, graph = graph)
  ev_bl <- evaluate_imputation(
    baseline$mu,
    pd$X_scaled,
    splits,
    trait_map = pd$trait_map
  )
  baseline_wall <- proc.time()[["elapsed"]] - t0

  graph_fit <- graph
  graph_fit$D <- NULL
  invisible(gc(full = TRUE, verbose = FALSE))

  t1 <- proc.time()[["elapsed"]]
  fit <- fit_pigauto(
    pd,
    tree,
    splits = splits,
    baseline = baseline,
    graph = graph_fit,
    epochs = epochs,
    verbose = FALSE,
    seed = seed
  )
  pred <- stats::predict(fit, return_se = TRUE, n_imputations = 1L)
  ev_pg <- evaluate_imputation(pred, pd$X_scaled, splits)
  pigauto_wall <- proc.time()[["elapsed"]] - t1

  rbind(
    standardize_eval(
      ev_bl,
      "baseline",
      scale_n,
      rep_id,
      targets,
      baseline_wall
    ),
    standardize_eval(
      ev_pg,
      "pigauto",
      scale_n,
      rep_id,
      targets,
      pigauto_wall,
      gate_bm = gate_map(fit, pd, targets, "r_cal_bm"),
      gate_gnn = gate_map(fit, pd, targets, "r_cal_gnn"),
      gate_mean = gate_map(fit, pd, targets, "r_cal_mean")
    )
  )
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

write_summary <- function(results, metadata, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("# TabPFN Reconciliation Benchmark\n", con)
  writeLines(
    sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    con
  )
  writeLines(sprintf("- Git commit: `%s`", metadata$commit), con)
  writeLines(
    sprintf("- Scales: `%s`", paste(metadata$scales, collapse = ", ")),
    con
  )
  writeLines(
    sprintf("- Regimes: `%s`", paste(metadata$regimes, collapse = ", ")),
    con
  )
  writeLines(sprintf("- Replicates: `%d`", metadata$n_reps), con)
  writeLines(sprintf("- Dry run: `%s`", metadata$dry_run), con)
  writeLines("", con)

  writeLines("## Aim", con)
  writeLines("", con)
  writeLines(
    paste(
      "This benchmark tests whether Russell-style TabPFN gains are mostly",
      "explained by same-row cross-trait features, phylogenetic features,",
      "or their combination. The shuffled regime keeps the same feature",
      "distribution but breaks row-level cross-trait alignment."
    ),
    con
  )
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
  if (!nrow(scored)) {
    writeLines("## Test Summary", con)
    writeLines("", con)
    writeLines("No scored rows were produced. Inspect `Status Counts`.", con)
    return(invisible(NULL))
  }

  writeLines("## Mean RMSE By Method", con)
  writeLines("", con)
  groups <- split(
    scored,
    interaction(scored$method, scored$scale_n, drop = TRUE)
  )
  agg <- do.call(
    rbind,
    lapply(groups, function(x) {
      data.frame(
        method = x$method[1L],
        scale_n = x$scale_n[1L],
        rmse = mean_or_na(x$rmse),
        pearson_r = mean_or_na(x$pearson_r),
        coverage_95 = mean_or_na(x$coverage_95),
        wall_sec = mean_or_na(x$wall_sec),
        stringsAsFactors = FALSE
      )
    })
  )
  agg <- agg[order(agg$scale_n, agg$rmse, agg$method), , drop = FALSE]
  writeLines(capture.output(print(agg, row.names = FALSE, digits = 4)), con)
  writeLines("", con)

  writeLines("## Mean RMSE By Trait", con)
  writeLines("", con)
  groups <- split(
    scored,
    interaction(scored$method, scored$scale_n, scored$trait, drop = TRUE)
  )
  trait_agg <- do.call(
    rbind,
    lapply(groups, function(x) {
      data.frame(
        method = x$method[1L],
        scale_n = x$scale_n[1L],
        trait = x$trait[1L],
        rmse = mean_or_na(x$rmse),
        pearson_r = mean_or_na(x$pearson_r),
        coverage_95 = mean_or_na(x$coverage_95),
        stringsAsFactors = FALSE
      )
    })
  )
  trait_agg <- trait_agg[
    order(trait_agg$scale_n, trait_agg$trait, trait_agg$rmse, trait_agg$method),
    ,
    drop = FALSE
  ]
  writeLines(
    capture.output(print(trait_agg, row.names = FALSE, digits = 4)),
    con
  )
  writeLines("", con)

  gate_rows <- scored[scored$method == "pigauto", , drop = FALSE]
  if (nrow(gate_rows)) {
    writeLines("## Pigauto Gate Audit", con)
    writeLines("", con)
    gate_groups <- split(
      gate_rows,
      interaction(gate_rows$scale_n, gate_rows$trait, drop = TRUE)
    )
    gate_agg <- do.call(
      rbind,
      lapply(gate_groups, function(x) {
        data.frame(
          scale_n = x$scale_n[1L],
          trait = x$trait[1L],
          rmse = mean_or_na(x$rmse),
          r_cal_bm = mean_or_na(x$r_cal_bm),
          r_cal_gnn = mean_or_na(x$r_cal_gnn),
          r_cal_mean = mean_or_na(x$r_cal_mean),
          stringsAsFactors = FALSE
        )
      })
    )
    gate_agg <- gate_agg[
      order(gate_agg$scale_n, gate_agg$trait),
      ,
      drop = FALSE
    ]
    writeLines(
      capture.output(print(gate_agg, row.names = FALSE, digits = 4)),
      con
    )
    writeLines("", con)
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
    splits <- make_missing_splits(
      pd$X_scaled,
      missing_frac = missing_frac,
      val_frac = val_frac,
      seed = seed,
      trait_map = pd$trait_map
    )

    continuous_targets <- names(Filter(
      function(tm) {
        identical(tm$type, "continuous")
      },
      pd$trait_map
    ))
    if (length(target_env)) {
      continuous_targets <- intersect(continuous_targets, target_env)
    }
    if (!length(continuous_targets)) {
      stop("No continuous targets selected.")
    }

    log_line("Targets: ", paste(continuous_targets, collapse = ", "))

    all_results[[length(all_results) + 1L]] <- pigauto_methods(
      pd,
      dat$tree,
      graph,
      splits,
      continuous_targets,
      seed,
      scale_n,
      rep_id
    )

    for (regime in regimes) {
      for (target in continuous_targets) {
        log_line(
          "TabPFN ",
          regime,
          " target=",
          target,
          " n=",
          scale_n,
          " rep=",
          rep_id
        )
        row <- tryCatch(
          tabpfn_one_target(
            pd,
            graph,
            splits,
            target,
            regime,
            seed,
            scale_n,
            rep_id
          ),
          error = function(e) {
            empty_row(
              paste0("tabpfn_", regime),
              scale_n,
              rep_id,
              target,
              paste0("error: ", conditionMessage(e))
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
