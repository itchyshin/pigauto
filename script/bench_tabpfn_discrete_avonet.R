#!/usr/bin/env Rscript
#
# script/bench_tabpfn_discrete_avonet.R
#
# Experimental TabPFN classification benchmark on AVONET non-continuous traits.
#
# This is a branch-local benchmark, not package API. It asks a narrower
# question than pigauto's mixed-type workflow: can TabPFNClassifier predict
# AVONET categorical/ordinal trait classes from cross-trait and phylogenetic
# features on the same held-out cells used for pigauto evaluation?
#
# Smoke check:
#   PIGAUTO_TABPFN_SCALES=50 PIGAUTO_TABPFN_REPS=1 \
#     PIGAUTO_TABPFN_VARIANTS=lappe_nfa PIGAUTO_TABPFN_RUN_PIGAUTO=false \
#     Rscript script/bench_tabpfn_discrete_avonet.R
#
# Local overnight run:
#   PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_discrete_avonet_local \
#   PIGAUTO_TABPFN_SCALES=50,75,300 PIGAUTO_TABPFN_REPS=3 \
#   PIGAUTO_TABPFN_VARIANTS=plain,lappe,lappe_nfa,knn \
#   PIGAUTO_TABPFN_RUN_PIGAUTO=true \
#     Rscript script/bench_tabpfn_discrete_avonet.R

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
out_stem <- env_chr("PIGAUTO_TABPFN_OUT_STEM",
                    "bench_tabpfn_discrete_avonet")
if (!grepl("^[A-Za-z0-9_.-]+$", out_stem)) {
  stop("PIGAUTO_TABPFN_OUT_STEM must be a simple file stem, not a path.")
}
out_rds <- file.path(repo_root, "script", paste0(out_stem, ".rds"))
out_md <- file.path(repo_root, "script", paste0(out_stem, ".md"))
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
  benchmark = "tabpfn_discrete_avonet",
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
  out_stem = out_stem,
  commit = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                    error = function(e) "unknown")
)

mode_int <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_integer_)
  as.integer(names(sort(table(x), decreasing = TRUE))[1L])
}

accuracy_vec <- function(truth, pred) {
  ok <- !is.na(truth) & !is.na(pred)
  if (!any(ok)) return(NA_real_)
  mean(truth[ok] == pred[ok])
}

balanced_accuracy_vec <- function(truth, pred) {
  ok <- !is.na(truth) & !is.na(pred)
  if (!any(ok)) return(NA_real_)
  per_class <- tapply(truth[ok] == pred[ok], truth[ok], mean)
  mean(per_class, na.rm = TRUE)
}

class_labels_from_truth <- function(X, tm) {
  n <- nrow(X)
  labels <- rep(NA_integer_, n)
  if (tm$type == "categorical") {
    mat <- X[, tm$latent_cols, drop = FALSE]
    ok <- stats::complete.cases(mat)
    labels[ok] <- max.col(mat[ok, , drop = FALSE], ties.method = "first") - 1L
  } else if (tm$type == "ordinal") {
    z <- X[, tm$latent_cols[1L]]
    ok <- is.finite(z)
    cls <- round(z[ok] * tm$sd + tm$mean)
    cls <- pmax(pmin(cls, length(tm$levels) - 1L), 0L)
    labels[ok] <- as.integer(cls)
  } else if (tm$type == "binary") {
    z <- X[, tm$latent_cols[1L]]
    ok <- is.finite(z)
    labels[ok] <- as.integer(round(z[ok]))
  } else {
    stop("Unsupported discrete target type: ", tm$type)
  }
  labels
}

class_labels_from_pred <- function(pred_mat, tm) {
  n <- nrow(pred_mat)
  labels <- rep(NA_integer_, n)
  if (tm$type == "categorical") {
    mat <- pred_mat[, tm$latent_cols, drop = FALSE]
    ok <- stats::complete.cases(mat)
    labels[ok] <- max.col(mat[ok, , drop = FALSE], ties.method = "first") - 1L
  } else if (tm$type == "ordinal") {
    z <- pred_mat[, tm$latent_cols[1L]]
    ok <- is.finite(z)
    cls <- round(z[ok] * tm$sd + tm$mean)
    cls <- pmax(pmin(cls, length(tm$levels) - 1L), 0L)
    labels[ok] <- as.integer(cls)
  } else if (tm$type == "binary") {
    z <- pred_mat[, tm$latent_cols[1L]]
    ok <- is.finite(z)
    labels[ok] <- as.integer(1 / (1 + exp(-z[ok])) >= 0.5)
  } else {
    stop("Unsupported discrete target type: ", tm$type)
  }
  labels
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

neighbor_class_features <- function(labels, D, train_rows, k, n_classes) {
  n <- length(labels)
  props <- matrix(NA_real_, nrow = n, ncol = n_classes)
  colnames(props) <- paste0("nfa_p_class_", seq_len(n_classes) - 1L)
  extra <- matrix(NA_real_, nrow = n, ncol = 3L)
  colnames(extra) <- c("nfa_entropy", "nfa_min_dist", "nfa_n")

  train_rows <- train_rows[!is.na(labels[train_rows])]
  if (length(train_rows) < 2L) return(cbind(props, extra))

  for (i in seq_len(n)) {
    candidates <- setdiff(train_rows, i)
    if (!length(candidates)) next
    d <- D[i, candidates]
    ord <- order(d, decreasing = FALSE, na.last = NA)
    if (!length(ord)) next
    keep <- candidates[ord[seq_len(min(k, length(ord)))]]
    labs <- labels[keep]
    tab <- tabulate(labs + 1L, nbins = n_classes)
    p <- tab / sum(tab)
    props[i, ] <- p
    extra[i, "nfa_entropy"] <- -sum(ifelse(p > 0, p * log(p), 0))
    extra[i, "nfa_min_dist"] <- min(D[i, keep], na.rm = TRUE)
    extra[i, "nfa_n"] <- length(keep)
  }

  cbind(props, extra)
}

feature_matrix <- function(pd, graph, splits, target_cols, labels,
                           train_rows, variant, nfa_k, n_classes) {
  X <- pd$X_scaled
  X[c(splits$val_idx, splits$test_idx)] <- NA_real_
  X <- X[, setdiff(seq_len(ncol(X)), target_cols), drop = FALSE]

  parts <- list(traits = X)
  if (variant %in% c("lappe", "lappe_nfa")) {
    coords <- graph$coords
    colnames(coords) <- paste0("lappe_", seq_len(ncol(coords)))
    parts$lappe <- coords
  }
  if (variant %in% c("lappe_nfa", "knn")) {
    parts$nfa <- neighbor_class_features(labels, graph$D, train_rows,
                                         nfa_k, n_classes)
  }

  out <- do.call(cbind, parts)
  storage.mode(out) <- "double"
  out
}

run_tabpfn_classify <- function(train_x, train_y, pred_x, seed) {
  tmpdir <- tempfile("tabpfn_discrete_")
  dir.create(tmpdir, recursive = TRUE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  train_path <- file.path(tmpdir, "train.csv")
  pred_path <- file.path(tmpdir, "predict.csv")
  out_path <- file.path(tmpdir, "pred.csv")
  meta_path <- file.path(tmpdir, "metadata.json")

  train_frame <- data.frame(.target = as.integer(train_y),
                            train_x, check.names = FALSE)
  pred_frame <- data.frame(pred_x, check.names = FALSE)
  utils::write.csv(train_frame, train_path, row.names = FALSE, na = "")
  utils::write.csv(pred_frame, pred_path, row.names = FALSE, na = "")

  cmd_args <- c(
    shQuote(runner_py),
    "--task", "classification",
    "--train", shQuote(train_path),
    "--predict", shQuote(pred_path),
    "--target", ".target",
    "--out", shQuote(out_path),
    "--metadata", shQuote(meta_path),
    "--device", device,
    "--seed", as.character(seed)
  )

  cmd_out <- system2(python, cmd_args, stdout = TRUE, stderr = TRUE)
  status <- attr(cmd_out, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(cmd_out, collapse = "\n"))
  }

  as.integer(utils::read.csv(out_path)$prediction)
}

result_row <- function(method, scale_n, rep_id, target, type, n_train, n_val,
                       n_test, n_features, truth, pred, wall_sec, status) {
  data.frame(
    method = method,
    scale_n = scale_n,
    rep = rep_id,
    trait = target,
    type = type,
    split = "test",
    n_train = n_train,
    n_val = n_val,
    n_test = n_test,
    n_features = n_features,
    accuracy = accuracy_vec(truth, pred),
    balanced_accuracy = balanced_accuracy_vec(truth, pred),
    wall_sec = wall_sec,
    status = status,
    stringsAsFactors = FALSE
  )
}

tabpfn_one_target <- function(pd, graph, splits, target, variant, seed,
                              scale_n, rep_id) {
  tm <- pd$trait_map[[target]]
  n <- nrow(pd$X_scaled)
  labels <- class_labels_from_truth(pd$X_scaled, tm)
  rows <- target_rows_from_split(splits, tm$latent_cols[1L], n)
  val_rows <- rows$val[!is.na(labels[rows$val])]
  test_rows <- rows$test[!is.na(labels[rows$test])]
  train_rows <- setdiff(which(!is.na(labels)), c(val_rows, test_rows))
  n_classes <- length(tm$levels)

  x <- feature_matrix(pd, graph, splits, tm$latent_cols, labels, train_rows,
                      variant, nfa_k, n_classes)

  if (length(unique(labels[train_rows])) < 2L || length(test_rows) < 2L) {
    return(result_row(
      paste0("tabpfn_", variant), scale_n, rep_id, target, tm$type,
      length(train_rows), length(val_rows), length(test_rows), ncol(x),
      labels[test_rows], rep(NA_integer_, length(test_rows)), NA_real_,
      "skipped_too_few_classes"
    ))
  }
  if (dry_run) {
    return(result_row(
      paste0("tabpfn_", variant), scale_n, rep_id, target, tm$type,
      length(train_rows), length(val_rows), length(test_rows), ncol(x),
      labels[test_rows], rep(NA_integer_, length(test_rows)), NA_real_,
      "dry_run"
    ))
  }

  t0 <- proc.time()[["elapsed"]]
  fit_result <- tryCatch({
    pred_rows <- c(val_rows, test_rows)
    pred_all <- run_tabpfn_classify(x[train_rows, , drop = FALSE],
                                    labels[train_rows],
                                    x[pred_rows, , drop = FALSE],
                                    seed)
    list(ok = TRUE,
         test_pred = pred_all[length(val_rows) + seq_along(test_rows)])
  }, error = function(e) {
    list(ok = FALSE, message = conditionMessage(e))
  })
  wall <- proc.time()[["elapsed"]] - t0

  if (!isTRUE(fit_result$ok)) {
    return(result_row(
      paste0("tabpfn_", variant), scale_n, rep_id, target, tm$type,
      length(train_rows), length(val_rows), length(test_rows), ncol(x),
      labels[test_rows], rep(NA_integer_, length(test_rows)), wall,
      paste0("error: ", fit_result$message)
    ))
  }

  result_row(
    paste0("tabpfn_", variant), scale_n, rep_id, target, tm$type,
    length(train_rows), length(val_rows), length(test_rows), ncol(x),
    labels[test_rows], fit_result$test_pred, wall, "ok"
  )
}

standard_methods <- function(pd, tree, graph, splits, targets, seed,
                             scale_n, rep_id) {
  rows <- list()

  t0 <- proc.time()[["elapsed"]]
  baseline <- fit_baseline(pd, tree, splits = splits, graph = graph)
  baseline_wall <- proc.time()[["elapsed"]] - t0

  fit <- NULL
  pred_pg <- NULL
  pigauto_wall <- NA_real_
  if (run_pigauto) {
    graph_fit <- graph
    graph_fit$D <- NULL
    invisible(gc(full = TRUE, verbose = FALSE))
    t1 <- proc.time()[["elapsed"]]
    fit <- fit_pigauto(pd, tree, splits = splits, baseline = baseline,
                       graph = graph_fit, epochs = epochs, verbose = FALSE,
                       seed = seed)
    pred_pg <- stats::predict(fit, return_se = TRUE, n_imputations = 1L)
    pigauto_wall <- proc.time()[["elapsed"]] - t1
  }

  n <- nrow(pd$X_scaled)
  for (target in targets) {
    tm <- pd$trait_map[[target]]
    labels <- class_labels_from_truth(pd$X_scaled, tm)
    split_rows <- target_rows_from_split(splits, tm$latent_cols[1L], n)
    test_rows <- split_rows$test[!is.na(labels[split_rows$test])]
    val_rows <- split_rows$val[!is.na(labels[split_rows$val])]
    train_rows <- setdiff(which(!is.na(labels)), c(val_rows, test_rows))
    truth <- labels[test_rows]

    maj <- mode_int(labels[train_rows])
    rows[[length(rows) + 1L]] <- result_row(
      "majority", scale_n, rep_id, target, tm$type,
      length(train_rows), length(val_rows), length(test_rows), NA_integer_,
      truth, rep(maj, length(test_rows)), 0, "ok"
    )

    pred_bl <- class_labels_from_pred(baseline$mu, tm)
    rows[[length(rows) + 1L]] <- result_row(
      "baseline", scale_n, rep_id, target, tm$type,
      length(train_rows), length(val_rows), length(test_rows), NA_integer_,
      truth, pred_bl[test_rows], baseline_wall, "ok"
    )

    if (run_pigauto) {
      pred_pg_lab <- class_labels_from_pred(pred_pg$imputed_latent, tm)
      rows[[length(rows) + 1L]] <- result_row(
        "pigauto", scale_n, rep_id, target, tm$type,
        length(train_rows), length(val_rows), length(test_rows), NA_integer_,
        truth, pred_pg_lab[test_rows], pigauto_wall, "ok"
      )
    }
  }

  do.call(rbind, rows)
}

write_summary <- function(results, metadata, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("# TabPFN AVONET Discrete Benchmark\n", con)
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
    "AVONET categorical/ordinal targets are held out with pigauto's",
    "`make_missing_splits()` and scored as classification tasks on the same",
    "test rows used for baseline and pigauto. TabPFN feature sets use",
    "same-row latent traits, optional Laplacian phylogenetic eigenvectors,",
    "and optional nearest-training-target class-proportion features.",
    "This benchmark reports accuracy and balanced accuracy only; it does not",
    "attempt classification prediction sets."
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
  writeLines("## Test Summary", con)
  writeLines("", con)
  if (nrow(scored)) {
    groups <- split(scored, interaction(scored$method, scored$scale_n,
                                        scored$trait, drop = TRUE))
    agg <- do.call(rbind, lapply(groups, function(x) {
      data.frame(
        method = x$method[1L],
        scale_n = x$scale_n[1L],
        trait = x$trait[1L],
        type = x$type[1L],
        accuracy = mean(x$accuracy, na.rm = TRUE),
        balanced_accuracy = mean(x$balanced_accuracy, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    agg <- agg[order(agg$scale_n, agg$trait, -agg$accuracy, agg$method), ,
               drop = FALSE]
    writeLines(capture.output(print(agg, row.names = FALSE, digits = 4)), con)
  } else if (isTRUE(metadata$dry_run)) {
    writeLines("No scored rows yet. This is expected for dry runs.", con)
  } else {
    writeLines("No scored rows were produced. Inspect `Status Counts`.", con)
  }
}

all_results <- list()

for (scale_n in scales) {
  for (rep_id in seq_len(n_reps)) {
    seed <- seed0 + scale_n * 100L + rep_id
    log_line("Preparing AVONET discrete n=", scale_n, " rep=", rep_id)
    dat <- load_avonet_subset(scale_n, seed)
    pd <- preprocess_traits(dat$df, dat$tree, log_transform = TRUE)
    graph <- build_phylo_graph(dat$tree, k_eigen = "auto")
    splits <- make_missing_splits(pd$X_scaled, missing_frac = missing_frac,
                                  val_frac = val_frac, seed = seed,
                                  trait_map = pd$trait_map)

    targets <- names(Filter(function(tm) {
      tm$type %in% c("categorical", "ordinal", "binary")
    }, pd$trait_map))
    if (length(target_env)) {
      targets <- intersect(targets, target_env)
    }
    if (!length(targets)) {
      stop("No discrete targets selected.")
    }

    log_line("Targets: ", paste(targets, collapse = ", "))
    all_results[[length(all_results) + 1L]] <-
      standard_methods(pd, dat$tree, graph, splits, targets, seed,
                       scale_n, rep_id)

    for (variant in variants) {
      for (target in targets) {
        log_line("TabPFN ", variant, " target=", target,
                 " n=", scale_n, " rep=", rep_id)
        all_results[[length(all_results) + 1L]] <- tabpfn_one_target(
          pd, graph, splits, target, variant, seed, scale_n, rep_id
        )
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
