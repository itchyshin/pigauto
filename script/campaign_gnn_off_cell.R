# script/campaign_gnn_off_cell.R
#
# One cell of the with/without-GNN campaign (plan: ~/.claude/plans/valiant-forging-gray.md,
# section ARC B). Runs every arm on ONE (dgp, n, seed) cell with the SAME user-level 30% MCAR
# mask, scores on the masked cells, and writes one rds. A driver fans cells out over cores.
#
# Usage:
#   Rscript script/campaign_gnn_off_cell.R --dgp bm_mixed --n 100 --seed 1 --out results/ \
#           [--arms gnn_on,gnn_off,gnn_off_pure,rphylopars,bace,floor] \
#           [--epochs 2000] [--bace_nitt 50000 --bace_burnin 10000 --bace_thin 25] [--smoke]
#
# DGPs: bm_mixed (4 continuous BM + 1 binary + 1 categorical(3) on ape::rcoal(n)),
#       ou_mixed (simulate_non_bm OU continuous + same discrete), avonet (bundled avonet300/tree300,
#       n ignored). Arms and metrics as locked in the plan. `--smoke` shrinks epochs/chains so the
#       whole cell runs in seconds; it is for checking the invocation, never for results.

suppressPackageStartupMessages({
  library(pigauto)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) return(TRUE)
  args[i + 1L]
}
dgp    <- get_arg("--dgp", "bm_mixed")
n      <- as.integer(get_arg("--n", 100L))
seed   <- as.integer(get_arg("--seed", 1L))
out    <- get_arg("--out", "results")
arms   <- strsplit(get_arg("--arms", "gnn_on,gnn_off,gnn_off_pure,rphylopars,bace,floor"), ",")[[1]]
# gnn_on_full (GNN-on predicting from the tax-free baseline) is derived from the gnn_on fit, not refit.
epochs <- as.integer(get_arg("--epochs", 2000L))
smoke  <- isTRUE(get_arg("--smoke", FALSE))
bace_nitt   <- as.integer(get_arg("--bace_nitt", 50000L))
bace_burnin <- as.integer(get_arg("--bace_burnin", 10000L))
bace_thin   <- as.integer(get_arg("--bace_thin", 25L))
miss_frac <- 0.30
if (smoke) { epochs <- 20L; bace_nitt <- 600L; bace_burnin <- 100L; bace_thin <- 5L }
dir.create(out, showWarnings = FALSE, recursive = TRUE)
tag <- sprintf("%s_n%d_s%d", dgp, n, seed)
log_line <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

# ---- data ---------------------------------------------------------------------------------
make_dgp <- function(dgp, n, seed) {
  set.seed(seed)
  if (dgp == "avonet") {
    e <- new.env(); utils::data("avonet300", package = "pigauto", envir = e)
    utils::data("tree300", package = "pigauto", envir = e)
    df <- e$avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
    return(list(df = df, tree = e$tree300))
  }
  tree <- ape::rcoal(n)   # ultrametric: BM-appropriate, and BACE/MCMCglmm requires it
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  sp <- tree$tip.label
  cont <- if (dgp == "ou_mixed") {
    pigauto::simulate_non_bm(tree, n_traits = 4L, scenario = "OU", seed = seed)
  } else {
    data.frame(row.names = sp,
               t1 = ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 0),
               t2 = ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 1),
               t3 = ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 2),
               t4 = ape::rTraitCont(tree, model = "BM", sigma = 1, root.value = 3))
  }
  cont <- as.data.frame(cont)[sp, , drop = FALSE]
  names(cont) <- paste0("c", seq_len(ncol(cont)))
  lat_b <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  bin <- factor(ifelse(lat_b > stats::median(lat_b), "yes", "no"))
  lat_k <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  cat3 <- factor(cut(lat_k, breaks = stats::quantile(lat_k, c(0, 1/3, 2/3, 1)),
                     labels = c("A", "B", "C"), include.lowest = TRUE))
  df <- cbind(cont, data.frame(bin = bin, cat3 = cat3, row.names = sp))
  list(df = df, tree = tree)
}

d <- make_dgp(dgp, n, seed)
truth <- d$df; tree <- d$tree
set.seed(seed + 1000L)
mask <- matrix(FALSE, nrow(truth), ncol(truth), dimnames = dimnames(truth))
for (v in names(truth)) {
  obs <- which(!is.na(truth[[v]])); hide <- sample(obs, ceiling(miss_frac * length(obs)))
  mask[hide, v] <- TRUE
}
df_miss <- truth; for (v in names(truth)) df_miss[mask[, v], v] <- NA
cont_traits <- names(truth)[vapply(truth, is.numeric, logical(1))]
log_line("cell %s: n=%d traits=%d (continuous %d) masked=%d", tag, nrow(truth), ncol(truth),
         length(cont_traits), sum(mask))

# ---- scoring ------------------------------------------------------------------------------
score_arm <- function(arm, completed, lower = NULL, upper = NULL) {
  rows <- list()
  for (v in names(truth)) {
    idx <- which(mask[, v]); if (!length(idx)) next
    if (is.numeric(truth[[v]])) {
      tr <- truth[[v]][idx]; pr <- completed[[v]][idx]
      train <- truth[[v]][!mask[, v] & !is.na(truth[[v]])]
      z <- sqrt(mean(((tr - pr) / stats::sd(train))^2))
      cov <- NA_real_
      if (!is.null(lower) && v %in% colnames(lower)) {
        rn <- rownames(truth)[idx]; cov <- mean(tr >= lower[rn, v] & tr <= upper[rn, v])
      }
      rows[[v]] <- data.frame(arm = arm, trait = v, metric = "zRMSE", value = z, coverage = cov)
    } else {
      acc <- mean(as.character(truth[[v]][idx]) == as.character(completed[[v]][idx]))
      rows[[v]] <- data.frame(arm = arm, trait = v, metric = "accuracy", value = acc, coverage = NA_real_)
    }
  }
  do.call(rbind, rows)
}

# ---- arms ---------------------------------------------------------------------------------
run_pigauto <- function(arm) {
  extra <- switch(arm,
    gnn_on       = list(gnn = TRUE,  epochs = epochs),
    gnn_off      = list(gnn = FALSE),
    gnn_off_pure = list(gnn = FALSE, safety_floor = FALSE, phylo_signal_gate = FALSE))
  res <- do.call(pigauto::impute, c(list(traits = df_miss, tree = tree, verbose = FALSE, seed = seed), extra))
  path <- res$fit$baseline$path
  pred <- res$prediction
  comp <- res$completed[rownames(truth), names(truth)]
  out <- list(completed = comp, lower = pred$conformal_lower, upper = pred$conformal_upper, path = path)
  if (arm == "gnn_on") {
    # Plan arm 2: the SAME GNN-on fit predicting from the tax-free baseline (baseline_override),
    # so the GNN effect and the held-out-cell tax can be separated. No refit.
    bf <- pigauto::fit_baseline(res$data, tree, splits = NULL,
                                lambda_mode = res$fit$model_config$lambda_mode %||% "fixed_1",
                                joint_solver = res$fit$model_config$joint_solver %||% "inhouse")
    pred_f <- stats::predict(res$fit, return_se = TRUE, baseline_override = bf)
    comp_f <- comp
    for (v in names(truth)) comp_f[mask[, v], v] <- pred_f$imputed[rownames(truth)[mask[, v]], v]
    out$derived <- list(gnn_on_full = list(completed = comp_f, lower = pred_f$conformal_lower,
                                           upper = pred_f$conformal_upper, path = path))
  }
  out
}

run_rphylopars <- function() {
  df4 <- df_miss[, cont_traits, drop = FALSE]
  df_in <- data.frame(species = rownames(df4), df4, stringsAsFactors = FALSE)
  fit <- Rphylopars::phylopars(df_in, tree = tree, model = "BM", phylo_correlated = TRUE,
                               pheno_correlated = TRUE, REML = TRUE)
  comp <- truth; comp[] <- NA
  rec <- fit$anc_recon[rownames(df4), cont_traits, drop = FALSE]
  for (v in cont_traits) { comp[[v]] <- df_miss[[v]]; comp[mask[, v], v] <- rec[mask[, v], v] }
  list(completed = comp)
}

run_bace <- function() {
  tree_b <- tree; if (any(tree_b$edge.length == 0)) tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
  df_b <- df_miss; df_b$Species <- rownames(df_miss)
  all_traits <- setdiff(names(df_b), "Species")
  fixformula <- lapply(all_traits, function(v) paste0(v, " ~ ", paste(setdiff(all_traits, v), collapse = " + ")))
  outb <- BACE::bace(fixformula = fixformula, ran_phylo_form = "~ 1 |Species", phylo = tree_b,
                     data = df_b, nitt = bace_nitt, burnin = bace_burnin, thin = bace_thin,
                     runs = 2L, n_final = 5L, verbose = FALSE, skip_conv = TRUE, ovr_categorical = TRUE)
  sets <- if ("imputed_datasets" %in% names(outb)) outb$imputed_datasets else
          if ("imputed_data" %in% names(outb)) outb$imputed_data else list(outb$data)
  comp <- df_miss
  for (v in names(comp)) {
    idx <- which(mask[, v]); if (!length(idx)) next
    draws <- sapply(sets, function(s) s[[v]][idx])
    if (!is.matrix(draws)) draws <- matrix(draws, ncol = length(sets))
    if (is.numeric(comp[[v]])) comp[idx, v] <- apply(draws, 1, stats::median)
    else comp[idx, v] <- apply(draws, 1, function(x) names(which.max(table(as.character(x)))))
  }
  list(completed = comp)
}

run_floor <- function() {
  comp <- df_miss
  for (v in names(comp)) {
    obs <- df_miss[[v]][!is.na(df_miss[[v]])]
    fill <- if (is.numeric(comp[[v]])) mean(obs) else names(which.max(table(as.character(obs))))
    comp[mask[, v], v] <- fill
  }
  list(completed = comp)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
results <- list(); walls <- list(); errors <- list(); paths <- list()
for (arm in arms) {
  t0 <- Sys.time()
  r <- tryCatch({
    if (arm %in% c("gnn_on", "gnn_off", "gnn_off_pure")) run_pigauto(arm)
    else if (arm == "rphylopars") run_rphylopars()
    else if (arm == "bace") run_bace()
    else if (arm == "floor") run_floor()
    else stop("unknown arm ", arm)
  }, error = function(e) e)
  walls[[arm]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error")) { errors[[arm]] <- conditionMessage(r); log_line("%s ERROR %s", arm, errors[[arm]]); next }
  results[[arm]] <- score_arm(arm, r$completed, r$lower, r$upper)
  paths[[arm]] <- r$path
  for (dn in names(r$derived)) {
    dr <- r$derived[[dn]]
    results[[dn]] <- score_arm(dn, dr$completed, dr$lower, dr$upper); paths[[dn]] <- dr$path
    walls[[dn]] <- 0
  }
  log_line("%s done in %.1f s", arm, walls[[arm]])
}
tab <- do.call(rbind, results); rownames(tab) <- NULL
tab$dgp <- dgp; tab$n <- nrow(truth); tab$seed <- seed
cell <- list(tag = tag, dgp = dgp, n = nrow(truth), seed = seed, arms = arms, smoke = smoke,
             epochs = epochs, bace = c(nitt = bace_nitt, burnin = bace_burnin, thin = bace_thin),
             results = tab, walls = unlist(walls), errors = errors, paths = paths,
             pigauto_version = as.character(utils::packageVersion("pigauto")),
             host = Sys.info()[["nodename"]], time = Sys.time())
saveRDS(cell, file.path(out, paste0(tag, if (smoke) "_smoke" else "", ".rds")))
print(tab, digits = 3); print(round(unlist(walls), 1))
if (length(errors)) { cat("ERRORS:\n"); print(errors) }
