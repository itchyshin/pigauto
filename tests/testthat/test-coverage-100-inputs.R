test_that("held-out continuous, count and ordinal metrics use finite pairs", {
  truth <- cbind(x = c(1, 2, 3, NA), count = c(1, 2, 3, 4), ord = c(1, 2, 3, 4))
  pred <- cbind(x = c(2, 2, 5, 9), count = c(2, 2, 5, 4), ord = c(1, 3, 2, 4))
  map <- list(list(name = "x", type = "continuous", latent_cols = 1L),
              list(name = "count", type = "count", latent_cols = 2L),
              list(name = "ord", type = "ordinal", latent_cols = 3L))
  out <- eval_test_cells(pred, NULL, truth, rep(1:4, 3), rep(1:3, each = 4),
                         map, 4, "fixture", c(x = 1, count = 1, ord = 1))
  value <- function(trait, metric) out$value[out$trait == trait & out$metric == metric]
  expect_equal(value("x", "rmse"), sqrt(5/3))
  expect_equal(value("x", "mae"), 1)
  expect_equal(value("x", "conformal_coverage_95"), 2/3)
  expect_equal(value("count", "rmse"), sqrt(5/4))
  expect_equal(value("count", "mae"), 3/4)
  expect_equal(value("ord", "spearman_rho"), 0.8)
  expect_equal(unique(out$n_test[out$trait == "x"]), 3L)
  expect_equal(unique(out$n_test[out$trait == "count"]), 4L)
  expect_null(eval_test_cells(pred, NULL, truth, integer(), integer(), map, 4, "fixture"))
})

test_that("held-out discrete and compositional metrics match known values", {
  truth <- cbind(bin = c(0, 1, 1), a = c(1, 0, 0), b = c(0, 1, 0),
                 c = c(0, 0, 1), gate = c(0, 1, 1), mag = c(NA, log(3), log(5)),
                 clr1 = c(0, 0, 0), clr2 = c(0, 0, 0))
  pred <- truth
  pred[, 1] <- qlogis(c(.2, .8, .4))
  pred[, 2:4] <- rbind(c(2, 0, 0), c(1, 2, 0), c(0, 1, 2))
  pred[, 5] <- qlogis(c(.2, .8, .4))
  pred[, 6] <- log(3)
  pred[, 7:8] <- cbind(rep(log(2), 3), rep(-log(2), 3))
  map <- list(list(name = "bin", type = "binary", latent_cols = 1L),
              list(name = "cat", type = "categorical", latent_cols = 2:4, levels = c("a", "b", "c")),
              list(name = "zi", type = "zi_count", latent_cols = 5:6, mean = 0, sd = 1),
              list(name = "comp", type = "multi_proportion", latent_cols = 7:8,
                   n_latent = 2L, mean = c(0, 0), sd = c(1, 1)))
  out <- eval_test_cells(pred, NULL, truth, rep(1:3, 4), rep(c(1, 2, 5, 7), each = 3), map, 3, "fixture")
  value <- function(trait, metric) out$value[out$trait == trait & out$metric == metric]
  expect_equal(value("bin", "accuracy"), 2/3)
  expect_equal(value("bin", "brier"), (0.04 + 0.04 + 0.36)/3)
  expect_equal(value("cat", "accuracy"), 1)
  expect_equal(value("cat", "accuracy_b"), 1)
  expect_equal(value("zi", "zero_accuracy"), 2/3)
  expect_equal(value("zi", "rmse"), sqrt((.4^2 + .4^2 + 3.2^2)/3))
  expect_equal(value("comp", "aitchison"), sqrt(2) * log(2))
  expect_equal(value("comp", "rmse_clr"), log(2))
  expect_equal(value("comp", "simplex_mae"), .3)
  obj <- list(probabilities = list(bin = c(.9, .9, .9)))
  alt <- eval_test_cells(pred, obj, truth, 1:3, rep(1L, 3), map, 3, "fixture")
  expect_equal(alt$value[alt$metric == "accuracy"], 2/3)
  expect_equal(alt$value[alt$metric == "brier"], (.81 + .01 + .01)/3)
})

test_that("calibration and confusion matrices retain defined denominators", {
  cal <- calibration_df(c(0, 1, 1, 0), c(0, .5, 1, .75), 2)
  expect_equal(cal$n, c(2L, 2L))
  expect_equal(cal$obs_freq, c(.5, .5))
  expect_equal(cal$mean_pred, c(.25, .875))
  expect_warning(cal2 <- calibration_df(c(0, NA, 1), c(.1, .2, .9), 2), "1 incomplete")
  expect_equal(sum(cal2$n), 2L)
  expect_error(calibration_df(c(0, 2), c(.1, .9)), "truth")
  expect_error(calibration_df(0, 1.1), "prob")
  expect_error(calibration_df(0, .1, 0), "n_bins")
  cm <- confusion_matrix(c("a", "a", "b", "other"), c("a", "b", "b", "a"), levels = c("a", "b", "c"))
  expect_equal(cm$accuracy, 2/3)
  expect_equal(cm$per_class$precision[2], .5)
  expect_equal(cm$per_class$recall[1], .5)
  expect_true(is.na(cm$per_class$f1[3]))
})

test_that("trait files preserve species and declared ordinal order", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  write.csv(data.frame(species = c("sp1", "sp2"), grade = c("high", "low"), x = c(2, 4)), path, row.names = FALSE)
  out <- read_traits(path, ordinal_levels = list(grade = c("low", "high")))
  expect_equal(rownames(out), c("sp1", "sp2"))
  expect_equal(levels(out$grade), c("low", "high"))
  expect_true(is.ordered(out$grade))
  expect_error(read_traits(path, trait_cols = "missing"), "Trait columns not found")
  expect_error(read_traits(path, ordinal_levels = list(grade = c("low", "medium"))), "not included")
  expect_error(read_traits(data.frame(species = c("a", "a"), x = 1:2)), "Duplicate")
  expect_error(read_traits(data.frame(species = "a")), "No trait columns")
  expect_error(read_tree(paste0(path, ".missing")), "not found")
})

test_that("occurrence and climate summaries filter invalid observations", {
  records <- data.frame(decimalLatitude = c(0, 10, 20, 100, 30),
                        decimalLongitude = c(0, 20, 40, 10, 50),
                        hasGeospatialIssues = c(FALSE, FALSE, NA, FALSE, FALSE),
                        basisOfRecord = c("HUMAN_OBSERVATION", "HUMAN_OBSERVATION", "HUMAN_OBSERVATION", "HUMAN_OBSERVATION", "FOSSIL_SPECIMEN"))
  centroid <- .gbif_centroid_one(records)
  expect_equal(centroid, list(centroid_lat = 5, centroid_lon = 10, n_occurrences = 2L))
  expect_equal(.gbif_centroid_one(records[FALSE, ])$n_occurrences, 0L)
  expect_equal(.gbif_cache_key("A b/c"), "A_b_c")
  expect_equal(.wc_cache_key("A b/c"), "A_b_c")
  expect_error(.wc_cache_key(""), "non-empty")
  clim <- .wc_aggregate_one(data.frame(bio1 = c(1, 3, NA), bio2 = c(2, 6, 9)))
  expect_equal(clim$n_extracted, 2L)
  expect_equal(clim$bio_median, c(bio1 = 2, bio2 = 4))
  expect_equal(clim$bio_iqr, c(bio1 = 1, bio2 = 2))
  expect_true(all(is.na(.wc_aggregate_one(data.frame(bio1 = NA_real_))$bio_median)))
})

test_that("preprocessing rejects malformed groups, values and covariates", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  dat <- data.frame(x = c(1, 2, 3), y = c(.2, .5, .8), row.names = c("a", "b", "c"))
  expect_error(preprocess_traits(dat, tree, multi_proportion_groups = list(g = "x")), "needs >= 2")
  expect_error(preprocess_traits(dat, tree, multi_proportion_groups = list(g = c("x", "z"))), "not in")
  expect_error(preprocess_traits(dat, tree, multi_proportion_groups = list(x = c("x", "y"))), "collide")
  expect_error(preprocess_traits(dat, tree, multi_proportion_groups = list(g = c("x", "y"), h = c("x", "y"))), "more than one")
  bad <- dat; bad$x[1] <- -1
  expect_error(preprocess_traits(bad, tree, trait_types = c(x = "count")), "negative")
  bad <- dat; bad$y[1] <- 1.5
  expect_error(preprocess_traits(bad, tree, trait_types = c(y = "proportion")), "proportion")
  expect_error(preprocess_traits(dat, tree, trait_types = c(z = "continuous")), "not found")
  expect_error(preprocess_traits(dat, tree, covariates = data.frame(z = c(1, NA, 3), row.names = rownames(dat))), "fully observed")
  expect_error(preprocess_traits(dat, tree, covariates = data.frame(z = c(1, Inf, 3), row.names = rownames(dat))), "non-finite")
  expect_error(make_missing_splits(matrix(1, 3, 2), missing_frac = 1), "missing_frac")
  expect_error(make_missing_splits(matrix(1, 3, 2), val_frac = 0), "val_frac")
  expect_error(make_missing_splits(matrix(1, 3, 2), mechanism = "MAR_phylo"), "tree")
  expect_error(make_missing_splits(matrix(NA_real_, 3, 2)), "No observed")
})

test_that("wide evaluation scores each held-out trait and legacy matrix", {
  truth <- cbind(x = 1:3, cnt = 1:3, ord = 1:3, bin = c(0, 1, 1),
                 a = c(1, 0, 0), b = c(0, 1, 1), gate = c(0, 1, 1),
                 mag = c(NA, log(3), log(5)), c1 = c(0, 0, 0), c2 = c(0, 0, 0))
  pred <- truth
  pred[, 1:3] <- cbind(c(2, 2, 5), c(2, 2, 5), c(1, 3, 2))
  pred[, 4] <- qlogis(c(.2, .8, .4))
  pred[, 5:6] <- rbind(c(2, 0), c(0, 2), c(2, 0))
  pred[, 7] <- qlogis(c(.2, .8, .4)); pred[, 8] <- log(3)
  pred[, 9:10] <- cbind(rep(log(2), 3), rep(-log(2), 3))
  map <- list(list(name="x", type="continuous", latent_cols=1L),
              list(name="cnt", type="count", latent_cols=2L),
              list(name="ord", type="ordinal", latent_cols=3L),
              list(name="bin", type="binary", latent_cols=4L),
              list(name="cat", type="categorical", latent_cols=5:6),
              list(name="zi", type="zi_count", latent_cols=7:8, mean=0, sd=1),
              list(name="comp", type="multi_proportion", latent_cols=9:10,
                   n_latent=2L, mean=c(0,0), sd=c(1,1)))
  idx <- unlist(lapply(c(1,2,3,4,5,7,9), function(j) (j-1)*3 + 1:3))
  splits <- list(val_idx=idx, test_idx=idx)
  out <- evaluate_imputation(pred, truth, splits, pred_se=matrix(1, 3, 10), trait_map=map)
  val <- out[out$split == "val", ]
  row <- function(nm) val[val$trait == nm, ]
  expect_equal(nrow(out), 14L)
  expect_equal(row("x")$rmse, sqrt(5/3))
  expect_equal(row("x")$coverage_95, 2/3)
  expect_equal(row("cnt")$mae, 1)
  expect_equal(row("ord")$spearman_rho, .5)
  expect_equal(row("bin")$accuracy, 2/3)
  expect_equal(row("cat")$accuracy, 2/3)
  expect_equal(row("zi")$rmse, sqrt((.4^2 + .4^2 + 3.2^2)/3))
  expect_equal(row("comp")$simplex_mae, .3)
  legacy <- evaluate_imputation(pred[, 1:2], truth[, 1:2], list(val_idx=1:3, test_idx=4:6), matrix(1, 3, 2))
  expect_equal(legacy$trait, c("x", "cnt"))
  expect_equal(legacy$coverage_95, c(2/3, 2/3))
  expect_error(evaluate_imputation(pred, truth[-1, ], splits), "identical dimensions")
  expect_error(evaluate_imputation(as.data.frame(pred), truth, splits), "matrices")
  safe <- structure(list(imputed_latent=pred, trait_map=map, evaluation_safe_idx=1L), class="pigauto_pred")
  expect_error(evaluate_imputation(safe, truth, splits), "not held out")
  empty <- evaluate_imputation(pred, truth, list(val_idx=integer(), test_idx=integer()), trait_map=map)
  expect_equal(nrow(empty), 0L)
})

test_that("input checks report invalid identities and branch lengths", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  dat <- data.frame(x = c(1, 2, NA), row.names = c("a", "b", "c"))
  check <- function(d = dat, t = tree, ...) .check_pigauto_internal(d, t, probe_runtime = FALSE, ...)$check
  bad <- check(d = matrix(1), t = "tree")
  expect_equal(bad$status, "error")
  expect_setequal(bad$messages$code, c("traits_not_data_frame", "tree_not_phylo"))
  expect_true("ordinal_levels_invalid" %in% check(ordinal_levels = list(x = c("low", "high")))$messages$code)
  expect_true("species_column_missing" %in% check(species_col = "missing")$messages$code)
  no_len <- tree; no_len$edge.length <- NULL
  expect_true("branch_length_missing_or_wrong" %in% check(t = no_len)$messages$code)
  negative <- tree; negative$edge.length[1] <- -1
  expect_true("branch_length_invalid" %in% check(t = negative)$messages$code)
  zero <- tree; zero$edge.length[] <- 0
  expect_true("branch_length_all_zero" %in% check(t = zero)$messages$code)
  one_zero <- tree; one_zero$edge.length[1] <- 0
  expect_true("branch_length_zero" %in% check(t = one_zero)$messages$code)
  dup <- tree; dup$tip.label[1] <- dup$tip.label[2]
  expect_true("duplicate_tree_tip" %in% check(t = dup)$messages$code)
  expect_output(print(bad), "pigauto input check: error")
})

test_that("missingness weights reflect drivers and sparse fallbacks", {
  x <- cbind(c(-1, 0, 1), c(NA, NA, 2))
  mar <- generate_missing_probs(x, 3, 2, "MAR_trait", list(beta = 2), NULL, NULL)
  expect_true(mar[1] < mar[2] && mar[2] < mar[3])
  expect_equal(mar[1:3], mar[4:6])
  expect_null(generate_missing_probs(x, 3, 2, "MAR_trait", list(driver_col = 2), NULL, NULL))
  mnar <- generate_missing_probs(x, 3, 2, "MNAR", list(), NULL, NULL)
  expect_gt(mnar[1], mnar[2])
  expect_equal(mnar[4:6], rep(1, 3))
  map <- list(list(latent_cols = 1L), list(latent_cols = 2L))
  expect_equal(generate_missing_probs(x, 3, 2, "MNAR", list(), NULL, map), mnar)
  expect_error(generate_missing_probs(x, 3, 2, "unknown", list(), NULL, NULL), "Unknown mechanism")
  expect_error(generate_missing_probs(x, 3, 2, "MAR_phylo", list(), NULL, NULL), "tree")
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  set.seed(7)
  phy <- generate_missing_probs(matrix(1:8, 4, 2), 4, 2, "MAR_phylo", list(n_clades = 1), tree, NULL)
  expect_equal(sort(unique(phy)), c(.1, .7))
})

test_that("climate cache handles missing occurrences and complete raster cache", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  wc <- file.path(root, "wc2.1_10m"); dir.create(wc)
  file.create(file.path(wc, ".wc_complete"))
  file.create(file.path(wc, sprintf("wc2.1_10m_bio_%d.tif", 1:19)))
  expect_equal(.wc_download_rasters(root, verbose = FALSE, .download_fn = function(...) stop("network")), wc)
  expect_error(.wc_download_rasters(root, resolution = "bad", verbose = FALSE), "resolution")
  ext <- file.path(root, "extracts")
  out <- .wc_extract_one("a bird", root, ext)
  expect_equal(out$n_extracted, 0L)
  expect_true(all(is.na(out$bio_median)))
  expect_equal(.wc_extract_one("a bird", root, ext), out)
  saveRDS(list(centroid_lat = NA_real_, n_occurrences = 0L), file.path(root, "another.rds"))
  expect_equal(.wc_extract_one("another", root, ext)$n_extracted, 0L)
  saveRDS(list(centroid_lat = 1, centroid_lon = 2, n_occurrences = 1L), file.path(root, "valid.rds"))
  expect_error(.wc_extract_one("valid", root, ext), "rast_stack is NULL")
})

test_that("cross-validation folds and summaries preserve test units", {
  x <- cbind(c(1, 2, 3, NA), c(0, 1, 0, NA), c(1, 0, 1, NA))
  map <- list(list(name="x", type="continuous", latent_cols=1L),
              list(name="cat", type="categorical", latent_cols=2:3))
  folds <- make_cv_folds(x, 3, 42, map)
  expect_equal(sort(unlist(folds)), c(1:3, 5:7, 9:11))
  for (fold in folds) expect_setequal(fold[fold %in% 5:7] + 4L, fold[fold %in% 9:11])
  cell_folds <- make_cv_folds(x, 3, 42)
  expect_equal(sort(unlist(cell_folds)), which(!is.na(x)))
  expect_error(make_cv_folds(matrix(NA_real_, 2, 2), 2, 1), "No observed")
  wide <- data.frame(trait=c("x", "y"), type="continuous", rmse=c(1, 3), mae=c(2, NA))
  long <- pivot_eval_to_long(wide, 1, 2)
  expect_equal(nrow(long), 3L)
  expect_equal(long$value[long$trait == "x" & long$metric == "mae"], 2)
  expect_null(pivot_eval_to_long(data.frame(trait="x", type="continuous", rmse=NA_real_), 1, 1))
  agg <- aggregate_cv_results(rbind(long, transform(long, value = value + 2)))
  expect_equal(agg$mean[agg$trait == "x" & agg$metric == "rmse"], 2)
  expect_equal(agg$n_folds[agg$trait == "x" & agg$metric == "rmse"], 2L)
  expect_equal(nrow(aggregate_cv_results(long[FALSE, ])), 0L)
  expect_error(cross_validate(list(), NULL), "pigauto_data")
  expect_error(cross_validate(structure(list(), class="pigauto_data"), NULL), "phylo")
  tree <- ape::read.tree(text="((a:1,b:1):1,c:2);")
  expect_error(cross_validate(structure(list(), class="pigauto_data"), tree, k=1), ">= 2")
  expect_error(cross_validate(structure(list(), class="pigauto_data"), tree, k=2, holdout="species"), "k >= 3")
  expect_error(cross_validate(structure(list(), class="pigauto_data"), tree, seeds=NA), "finite numeric")
})

test_that("conformal fold coverage uses original trait scale", {
  x <- cbind(cont=c(0, 1), count=c(0, log(3)), prop=c(0, qlogis(.75)), ord=c(1, 2))
  map <- list(list(name="cont", type="continuous", latent_cols=1L, mean=1, sd=2, log_transform=FALSE),
              list(name="count", type="count", latent_cols=2L, mean=0, sd=1),
              list(name="prop", type="proportion", latent_cols=3L, mean=0, sd=1),
              list(name="ord", type="ordinal", latent_cols=4L, mean=0, sd=1))
  lo <- matrix(0, 2, 4, dimnames=list(NULL, colnames(x)))
  hi <- matrix(c(1, 2, 0, 2.1, .5, .6, 1, 3), 2, 4, dimnames=list(NULL, colnames(x)))
  out <- compute_conformal_coverage_fold(list(conformal_lower=lo, conformal_upper=hi), x, 1:8, map, 1, 1)
  expect_equal(out$coverage, c(.5, 1, .5, 1))
  expect_equal(out$n, rep(2L, 4))
  expect_null(compute_conformal_coverage_fold(list(), x, 1:8, map, 1, 1))
})

test_that("legacy scoring and persistence validate their inputs", {
  truth <- matrix(c(1,2,3,NA), 4, 1)
  pred <- matrix(c(2,2,5,0), 4, 1)
  out <- eval_test_cells(pred, NULL, truth, 1:4, rep(1L,4), NULL, 4, "legacy")
  expect_equal(out$trait, rep("trait1", 3))
  expect_equal(out$value[out$metric == "rmse"], sqrt(5/3))
  expect_equal(out$value[out$metric == "mae"], 1)
  expect_null(eval_test_cells(pred, NULL, truth, integer(), integer(), NULL, 4, "legacy"))
  expect_error(completed_data(list()), "pigauto_result")
  expect_error(save_pigauto(list(), tempfile()), "pigauto_fit")
  expect_error(load_pigauto(tempfile()), "not found")
  path <- tempfile(); on.exit(unlink(path))
  saveRDS(list(), path)
  expect_error(load_pigauto(path), "does not contain")
  fit <- structure(list(model_state=list(), marker=42), class="pigauto_fit")
  save_pigauto(fit, path, compress=FALSE)
  expect_equal(load_pigauto(path)$marker, 42)
})

test_that("benchmark summaries compare like metrics and tolerate sparse scenarios", {
  tbl <- data.frame(scenario=c("BM", "BM", "BM", "BM", "OU"),
                    trait=c("x", "x", "b", "b", "x"),
                    method=c("baseline", "pigauto", "baseline", "pigauto", "baseline"),
                    type=c("continuous", "continuous", "binary", "binary", "continuous"),
                    metric=c("rmse", "rmse", "accuracy", "accuracy", "rmse"),
                    mean=c(1, .8, .5, .75, 1.2), sd=rep(.1, 5))
  obj <- structure(list(scenarios=c("BM", "OU"), summary=tbl), class="pigauto_benchmark")
  expect_output(summary(obj), "20.0% RMSE", fixed=TRUE)
  expect_output(summary(obj), "75.0%", fixed=TRUE)
  path <- tempfile(fileext=".pdf"); grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add=TRUE)
  expect_null(plot(obj, metric="rmse"))
  expect_message(plot(obj, metric="brier"), "No results")
})

test_that("preprocessed data print describes grouped observations", {
  tree <- ape::read.tree(text="((a:1,b:1):1,c:2);")
  d <- data.frame(species=c("a", "a", "b", "c"), x=c(1, NA, 2, 3))
  pd <- preprocess_traits(d, tree, species_col="species")
  expect_output(print(pd), "Observations: 4", fixed=TRUE)
  expect_output(print(pd), "Missing values:", fixed=TRUE)
  legacy <- pd; legacy$trait_map <- NULL; legacy$log_transform <- FALSE
  expect_output(print(legacy), "Log-transformed: FALSE", fixed=TRUE)
})

test_that("small GNN-off benchmark returns scored held-out traits", {
  bench <- suppressWarnings(simulate_benchmark(n_species=12, n_traits=2, scenarios="BM",
                                                n_reps=1, seed=42, gnn=FALSE, verbose=FALSE))
  expect_s3_class(bench, "pigauto_benchmark")
  expect_setequal(unique(bench$results$scenario), "BM")
  expect_setequal(unique(bench$results$method), c("pigauto", "baseline"))
  expect_true(all(bench$results$n_test > 0L))
  expect_true(any(is.finite(bench$results$value[bench$results$metric == "rmse"])))
  expect_true(all(bench$summary$n_reps <= 1L))
})

test_that("small GNN-off cross-validation scores each fold", {
  tree <- ape::read.tree(text="(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
  d <- data.frame(x=as.numeric(1:8), y=c(1,3,2,5,4,8,7,6), row.names=tree$tip.label)
  pd <- preprocess_traits(d, tree)
  cv <- suppressWarnings(cross_validate(pd, tree, k=3, seeds=1, gnn=FALSE, verbose=FALSE, k_eigen=3))
  expect_s3_class(cv, "pigauto_cv")
  expect_setequal(unique(cv$results$fold), 1:3)
  expect_setequal(unique(cv$results$trait), c("x", "y"))
  expect_true(all(is.finite(cv$results$value[cv$results$metric == "rmse"])))
  expect_output(summary(cv), "cross-validation")
  expect_output(print(cv), "cross-validation")
  expect_error(make_grouped_cv_folds(pd, tree, k=9, seed=1, holdout="species"), "cannot exceed")
  bad <- pd; bad$species_names <- rev(bad$species_names)
  expect_error(make_grouped_cv_folds(bad, tree, k=3, seed=1, holdout="species"), "same order")
  expect_error(prepare_grouped_cv_data(list(), list()), "freshly preprocessed")
})

test_that("simulation variants preserve tip alignment and requested traits", {
  tree <- ape::read.tree(text="((a:1,b:1):1,c:2);")
  set.seed(91)
  shifted <- simulate_regime_shift(tree, n_traits=2, sigma=1, shift_magnitude=2)
  expect_equal(rownames(shifted), tree$tip.label)
  expect_equal(ncol(shifted), 2L)
  set.seed(91)
  nonlinear <- simulate_nonlinear(tree, n_traits=5, sigma=1)
  expect_equal(rownames(nonlinear), tree$tip.label)
  expect_equal(names(nonlinear), paste0("trait", 1:5))
  expect_true(all(is.finite(as.matrix(nonlinear))))
})

test_that("preprocessing encodes covariates and reports malformed inputs", {
  tree <- ape::read.tree(text="((a:1,b:1):1,c:2);")
  dat <- data.frame(x=c(1,2,3), flag=c(TRUE,FALSE,TRUE), row.names=tree$tip.label)
  cov <- data.frame(group=c("a","b","a"), constant=c(5,5,5), row.names=tree$tip.label)
  pd <- preprocess_traits(dat, tree, covariates=cov, log_cols="x")
  expect_true(pd$trait_map$x$log_transform)
  expect_equal(pd$trait_map$flag$type, "binary")
  expect_true(all(pd$covariates[, "constant"] == 0))
  expect_equal(rowSums(pd$covariates[, c("group.a", "group.b")]), rep(1,3))
  expect_error(preprocess_traits(dat, tree, multi_proportion_groups=list(c("x", "flag"))), "named list")
  expect_error(preprocess_traits(dat, tree, multi_proportion_groups=list(g=c("x", "flag"))), "non-numeric")
  expect_error(preprocess_traits(dat, tree, species_col="absent"), "not found")
  expect_error(preprocess_traits(dat, tree, covariates=1), "data.frame or matrix")
  expect_error(preprocess_traits(dat, tree, trait_types="continuous"), "named")
  bad <- cov; bad$group <- as.Date("2020-01-01") + 0:2
  expect_error(preprocess_traits(dat, tree, covariates=bad), "unsupported class")
  expect_error(detect_trait_types(data.frame(x=I(list(1,2,3)))), "Cannot determine")
})

test_that("trait and tree readers reject absent or incompatible input", {
  expect_error(read_traits(tempfile()), "Trait file not found")
  expect_error(read_traits(matrix(1, 2, 2)), "data.frame or a path")
  expect_error(read_traits(data.frame(x = 1)), "Column 'species' not found")
  converted <- read_traits(data.frame(species = c("a", "b", "c"),
                                     group = c("z", "a", "z"), rank = c("low", "high", "low")),
                           factor_cols = "group", ordered_cols = "rank")
  expect_identical(rownames(converted), c("a", "b", "c"))
  expect_identical(levels(converted$group), c("a", "z"))
  expect_true(is.ordered(converted$rank))
  expect_identical(as.character(converted$rank), c("low", "high", "low"))
})

test_that("preprocessing validates basic objects and encodes constant scales", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  d <- data.frame(x = c(2, 2, 2), row.names = tree$tip.label)
  expect_error(preprocess_traits(as.matrix(d), tree), "must be a data.frame")
  expect_error(preprocess_traits(d, list()), "must be a phylo object")
  for (type in c("count", "proportion")) {
    values <- if (type == "count") c(2, 2, 2) else c(.4, .4, .4)
    p <- preprocess_traits(data.frame(x = values, row.names = tree$tip.label), tree,
                           trait_types = c(x = type))
    expect_identical(p$trait_map[[1]]$sd, 1)
    expect_equal(as.numeric(p$X_scaled), rep(0, 3))
  }
  ord <- preprocess_traits(data.frame(x = ordered(rep("low", 3), levels = c("low", "medium", "high")),
                                      row.names = tree$tip.label), tree)
  expect_identical(ord$trait_map[[1]]$sd, 1)
  expect_equal(as.numeric(ord$X_scaled), rep(0, 3))
  covd <- preprocess_traits(d, tree, trait_types = c(x = "continuous"),
                            covariates = data.frame(flag = c(TRUE, FALSE, TRUE), row.names = tree$tip.label))
  expect_true(all(is.finite(covd$covariates)))
  expect_equal(nrow(covd$covariates), 3L)
})

test_that("scalar metric helpers preserve empty-pair and trait-map contracts", {
  expect_equal(accuracy_vec(c(1, NA, 0), c(1, 1, 1)), .5)
  expect_identical(accuracy_vec(NA_real_, 1), NA_real_)
  expect_identical(brier_vec(NA_real_, .5), NA_real_)
  expect_identical(mae_vec(c(NA, Inf), c(1, 2)), NA_real_)
  expect_equal(mae_vec(c(1, 3, Inf), c(2, 5, 0)), 1.5)
  map <- list(list(name = "x", type = "continuous", latent_cols = 1L),
              list(name = "class", type = "categorical", latent_cols = 2:4))
  expect_identical(get_trait_for_latent_col(3, map), map[[2]])
  expect_null(get_trait_for_latent_col(5, map))
  expect_identical(get_latent_cols_by_type(map, "categorical"), 2:4)
  expect_identical(get_latent_cols_by_type(map, "binary"), integer())
  expect_error(calibration_df(factor(c(0, 1)), c(.1, .9)), "0/1 or logical")
  expect_error(calibration_df(c(0, 1), c("low", "high")), "numeric vector")
  expect_error(calibration_df(c(0, 1), .5), "same length")
  expect_warning(expect_error(calibration_df(c(NA, NA), c(.1, .9)),
                             "no complete"), "2 incomplete pairs")
})
test_that("empty metric inputs and latent map lookups retain their contracts", {
  expect_true(is.na(accuracy_vec(c(NA, 1), c(0, NA))))
  expect_true(is.na(brier_vec(c(NA, 1), c(0, NA))))
  expect_true(is.na(mae_vec(c(Inf, NA), c(0, 1))))
  map <- list(list(type = "continuous", latent_cols = 1L),
              list(type = "categorical", latent_cols = 2:4),
              list(type = "continuous", latent_cols = 5L))
  expect_identical(get_trait_for_latent_col(3L, map), map[[2]])
  expect_null(get_trait_for_latent_col(6L, map))
  expect_identical(get_latent_cols_by_type(map, "continuous"), c(1L, 5L))
  expect_identical(get_latent_cols_by_type(map, "binary"), integer())
  x <- matrix(c(NA_real_, 2), ncol = 1)
  result <- evaluate_continuous_legacy(x, x, list(val_idx = 1L, test_idx = integer()), NULL)
  expect_null(result)
  cv <- structure(list(k = 2L, n_reps = 1L, summary = data.frame()), class = "pigauto_cv")
  expect_output(expect_identical(print(cv), cv), "No results")
})

test_that("mixed evaluation skips traits with no selected or usable truth", {
  types <- c("continuous", "binary", "categorical", "zi_count", "multi_proportion")
  cols <- list(1L, 2L, 3:4, 5:6, 7:8)
  map <- Map(function(type, lc) list(name = type, type = type, latent_cols = lc,
                                    mean = 0, sd = 1, n_latent = length(lc)), types, cols)
  truth <- matrix(NA_real_, 2, 8)
  pred <- matrix(0, 2, 8)
  for (j in c(1L, 2L, 3L, 5L, 7L)) {
    idx <- 1L + (j - 1L) * 2L
    out <- evaluate_imputation(pred, truth, list(val_idx = idx, test_idx = integer()), trait_map = map)
    expect_equal(nrow(out), 0L)
  }
  truth[2, 1] <- 4
  pred[2, 1] <- 5
  out <- evaluate_imputation(pred, truth, list(val_idx = 2L, test_idx = integer()), trait_map = map)
  expect_identical(out$trait, "continuous")
  expect_equal(out$rmse, 1)
  expect_true(is.na(out$pearson_r))
})

test_that("sparse folds and conformal scoring omit unavailable cells", {
  x <- cbind(empty = c(NA, NA, NA), one = c(2, NA, NA))
  map <- list(list(name = "empty", type = "continuous", latent_cols = 1L, mean = 0, sd = 1),
              list(name = "one", type = "continuous", latent_cols = 2L, mean = 0, sd = 1))
  folds <- make_cv_folds(x, k = 3L, seed = 42, trait_map = map)
  expect_equal(sort(unlist(folds)), 4L)
  expect_equal(sum(lengths(folds) == 0L), 2L)
  bounds <- list(conformal_lower = matrix(0, 3, 2, dimnames = list(NULL, colnames(x))),
                 conformal_upper = matrix(3, 3, 2, dimnames = list(NULL, colnames(x))))
  expect_null(compute_conformal_coverage_fold(bounds, x, 1L, map, 1L, 1L))
  result <- compute_conformal_coverage_fold(bounds, x, c(1L, 4L), map, 2L, 3L)
  expect_identical(result$trait, "one")
  expect_equal(result$coverage, 1)
  expect_equal(result$n, 1L)
  expect_equal(result$fold, 2L)
  expect_equal(result$rep, 3L)
  cv <- structure(list(k = 3L, n_reps = 1L, summary = data.frame(
    trait = c("grade", "flag"), type = c("ordinal", "binary"),
    metric = c("spearman_rho", "accuracy"), mean = c(.7, .8), sd = c(.1, .2), n_folds = 3L)),
    class = "pigauto_cv")
  txt <- paste(capture.output(summary(cv)), collapse = " ")
  expect_match(txt, "0.700.*rho")
  expect_match(txt, "0.800")
  empty <- cv; empty$summary <- cv$summary[FALSE, ]
  expect_output(expect_identical(summary(empty), empty$summary), "No results")
})

test_that("evaluation entry points reject invalid inputs before fitting", {
  expect_error(evaluate(list()), "pigauto_fit")
  expect_error(evaluate(structure(list(), class = "pigauto_fit")), "No test split")
  expect_error(compare_methods(list(), NULL), "pigauto_data")
  data <- structure(list(), class = "pigauto_data")
  expect_error(compare_methods(data, NULL), "phylo")
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  expect_error(compare_methods(data, tree, seeds = c(1, Inf)), "finite numeric")
})

test_that("empty occurrence and raster summaries retain their named shape", {
  rejected <- data.frame(decimalLatitude = c(100, NA_real_),
                         decimalLongitude = c(10, 20))
  expect_identical(.gbif_centroid_one(rejected),
                   list(centroid_lat = NA_real_, centroid_lon = NA_real_, n_occurrences = 0L))
  values <- data.frame(bio1 = numeric(), bio2 = numeric())
  result <- .wc_aggregate_one(values)
  expect_identical(result$bio_median, c(bio1 = NA_real_, bio2 = NA_real_))
  expect_identical(result$bio_iqr, c(bio1 = NA_real_, bio2 = NA_real_))
  expect_identical(result$n_extracted, 0L)
})

test_that("method comparison scores a supplied holdout with a null seed", {
  tree <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
  d <- data.frame(x = as.numeric(1:8), y = c(1,3,2,5,4,8,7,6), row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  splits <- make_missing_splits(pd$X_scaled, seed = 41, trait_map = pd$trait_map)
  messages <- character()
  expect_warning(out <- withCallingHandlers(
    compare_methods(pd, tree, splits = splits, seeds = NULL, gnn = FALSE,
                    verbose = TRUE, k_eigen = 3),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }), "Small validation set")
  expect_true(any(grepl("Evaluating", messages)))
  expect_setequal(unique(out$method), c("pigauto", "baseline"))
  expect_identical(unique(out$rep), 1L)
  expect_true(all(out$n_test > 0L))
  expect_true(any(is.finite(out$value[out$metric == "rmse"])))
})

test_that("missing external packages fail without network access or cache damage", {
  root <- tempfile("external-boundary-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  sp <- "Example species"
  path <- file.path(root, paste0(.gbif_cache_key(sp), ".rds"))
  original <- list(species = sp, centroid_lat = 2, centroid_lon = 3, n_occurrences = 1L)
  saveRDS(original, path)
  real_require <- base::requireNamespace
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (package %in% c("rgbif", "terra")) return(FALSE)
      real_require(package, ...)
    }, .package = "base")
  expect_error(.gbif_fetch_one(sp, root, refresh_cache = TRUE), "requires the 'rgbif' package")
  expect_identical(readRDS(path), original)
  expect_error(pull_worldclim_per_species(sp, root, file.path(root, "wc")),
               "requires the 'terra' package")
})

test_that("failed WorldClim download leaves existing cache intact", {
  root <- tempfile("wc-failure-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  wc <- file.path(root, "wc2.1_10m"); dir.create(wc)
  existing <- file.path(wc, "wc2.1_10m_bio_1.tif")
  writeLines("previous bytes", existing)
  calls <- 0L
  expect_message(expect_error(.wc_download_rasters(root, verbose = TRUE,
    .download_fn = function(url, destfile, mode) {
      calls <<- calls + 1L
      stop("offline fixture")
    }), "offline fixture"), "downloading")
  expect_identical(calls, 1L)
  expect_identical(readLines(existing), "previous bytes")
  expect_false(file.exists(file.path(wc, ".wc_complete")))
})

test_that("WorldClim raster failure falls back to cached extraction", {
  if (requireNamespace("terra", quietly = TRUE)) {
    root <- tempfile("wc-raster-"); dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    wc <- file.path(root, "wc2.1_10m"); dir.create(wc)
    paths <- file.path(wc, sprintf("wc2.1_10m_bio_%d.tif", 1:19))
    invisible(vapply(paths, function(path) { writeLines("invalid raster", path); TRUE }, logical(1)))
    file.create(file.path(wc, ".wc_complete"))
    extracts <- file.path(root, "extracts"); dir.create(extracts)
    med <- setNames(as.numeric(1:19), paste0("bio", 1:19))
    cached <- list(species = "bird", bio_median = med, bio_iqr = med / 2, n_extracted = 3L)
    saveRDS(cached, file.path(extracts, "bird.rds"))
    testthat::local_mocked_bindings(rast = function(...) stop("invalid raster fixture"),
                                    .package = "terra")
    expect_warning(out <- pull_worldclim_per_species("bird", root, root, verbose = FALSE),
                   "terra::rast failed.*cache-only")
    expect_equal(unname(unlist(out[1, paste0("bio", 1:19, "_median")])), unname(med))
    expect_equal(out$n_extracted, 3L)
    expect_identical(readRDS(file.path(extracts, "bird.rds")), cached)
  }
})

test_that("GBIF occurrence failure preserves a previous cache entry", {
  if (requireNamespace("rgbif", quietly = TRUE)) {
    root <- tempfile("gbif-offline-"); dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    sp <- "Example species"
    path <- file.path(root, paste0(.gbif_cache_key(sp), ".rds"))
    old <- list(species = sp, centroid_lat = 1, centroid_lon = 2, n_occurrences = 1L)
    saveRDS(old, path)
    testthat::local_mocked_bindings(name_backbone = function(...) list(usageKey = 1L),
                                    .package = "rgbif")
    testthat::local_mocked_bindings(.gbif_fetch_page = function(...) NULL,
                                    .package = "pigauto")
    expect_error(.gbif_fetch_one(sp, root, refresh_cache = TRUE, sleep_ms = 0L),
                 "occurrence request failed")
    expect_identical(readRDS(path), old)
  }
})

test_that("requesting GBIF points refreshes a centroid-only cache", {
  if (requireNamespace("rgbif", quietly = TRUE)) {
    root <- tempfile("gbif-points-"); dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    sp <- "Example species"
    path <- file.path(root, paste0(.gbif_cache_key(sp), ".rds"))
    saveRDS(list(species = sp, centroid_lat = 1, centroid_lon = 2, n_occurrences = 1L), path)
    testthat::local_mocked_bindings(name_backbone = function(...) list(matchType = "NONE"),
                                    .package = "rgbif")
    expect_message(out <- .gbif_fetch_one(sp, root, store_points = TRUE, sleep_ms = 0L),
                   "per-occurrence points requested but not found")
    expect_equal(out$n_occurrences, 0L)
    expect_identical(readRDS(path)$match_type, "NO_MATCH")
  }
})

test_that("runtime check reports an absent torch package", {
  real_require <- base::requireNamespace
  testthat::local_mocked_bindings(requireNamespace = function(package, ...) {
    if (package == "torch") return(FALSE)
    real_require(package, ...)
  }, .package = "base")
  result <- .check_pigauto_runtime()
  expect_false(result$torch_package)
  expect_true(result$unavailable)
  expect_identical(result$device, "unavailable")
  expect_match(result$detail, "package is not installed")
})

test_that("runtime check reports a missing torch runtime", {
  real_require <- base::requireNamespace
  testthat::local_mocked_bindings(requireNamespace = function(package, ...) {
    if (package == "torch") return(TRUE)
    real_require(package, ...)
  }, .package = "base")
  testthat::local_mocked_bindings(.check_pigauto_torch_is_installed = function() FALSE,
                                  .package = "pigauto")
  result <- .check_pigauto_runtime()
  expect_true(result$torch_package)
  expect_false(result$torch_is_installed)
  expect_true(result$unavailable)
  expect_identical(result$detail, "torch runtime is not installed")
})

test_that("tree check identifies an unrooted polytomy", {
  tree <- ape::read.tree(text = "(a:1,b:1,c:1,d:1);")
  d <- data.frame(x = c(1, 2, 3, 4), row.names = tree$tip.label)
  result <- .check_pigauto_internal(d, tree, probe_runtime = FALSE)$check
  expect_true("tree_unrooted" %in% result$messages$code)
  expect_true("tree_polytomy" %in% result$messages$code)
})

test_that("zero requested CV repetitions returns a typed empty summary", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  d <- data.frame(x = c(1, 2, 4, 3), row.names = tree$tip.label)
  pd <- preprocess_traits(d, tree, log_transform = FALSE)
  cv <- cross_validate(pd, tree, k = 2L, seeds = numeric(), verbose = FALSE,
                       gnn = FALSE, k_eigen = 2L)
  expect_s3_class(cv, "pigauto_cv")
  expect_equal(nrow(cv$results), 0L)
  expect_identical(names(cv$results), c("fold", "rep", "trait", "type", "metric", "value"))
  expect_type(cv$results$value, "double")
  expect_equal(nrow(cv$summary), 0L)
})
