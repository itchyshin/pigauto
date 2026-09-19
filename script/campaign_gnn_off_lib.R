# script/campaign_gnn_off_lib.R
# Shared pieces of the campaign runners: DGPs, the seeded user-level mask, and scoring.
# Sourced by script/campaign_gnn_off_cell.R and script/campaign_solver_cell.R.
# score_arm() reads `truth` and `mask` from the calling environment (set by the runner).

# ---- data ---------------------------------------------------------------------------------
make_dgp <- function(dgp, n, seed) {
  set.seed(seed)
  if (dgp == "avonet") {
    e <- new.env(); utils::data("avonet300", package = "pigauto", envir = e)
    utils::data("tree300", package = "pigauto", envir = e)
    df <- e$avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
    return(list(df = df, tree = e$tree300))
  }
  if (dgp == "bace_dgp") {
    # BACE's own simulator (as in script/sim_bace_dgp.R): gaussian response y with
    # phylogenetic signal 0.4 on three predictors (gaussian, gaussian, binary; signal 0.3),
    # birth-death tree. All four columns are treated as traits to impute; the campaign's
    # uniform 30% MCAR mask replaces BACE's response-only missingness.
    sim <- BACE::sim_bace(response_type = "gaussian",
                          predictor_types = c("gaussian", "gaussian", "binary"),
                          beta_resp = 0.5 * c(0.7, 0.4, 0.5), phylo_signal = c(0.4, 0.3, 0.3, 0.3),
                          n_cases = as.integer(n), n_species = as.integer(n),
                          missingness = c(0, 0, 0, 0), birth = 0.8, death = 0.4)
    df <- sim$complete_data; tree <- sim$tree
    rownames(df) <- as.character(df$Species); df$Species <- NULL
    df$x3 <- factor(as.character(df$x3), ordered = FALSE)   # binary, not ordered
    df <- df[tree$tip.label, , drop = FALSE]
    return(list(df = df, tree = tree))
  }
  if (dgp == "types_mixed") {
    # Every pigauto trait type on one tree, each from its own BM latent (independent liabilities):
    # 2 continuous, 1 count (Poisson, log link), 1 proportion (logit-normal), 1 binary, 1 ordinal
    # (4 levels), 1 categorical (3 levels). Used for the per-type small tests.
    tree <- ape::rcoal(n); sp <- tree$tip.label
    lat <- function() ape::rTraitCont(tree, model = "BM", sigma = 1)
    l_cnt <- lat(); l_prp <- lat(); l_bin <- lat(); l_ord <- lat(); l_cat <- lat()
    ord_breaks <- stats::quantile(l_ord, c(0, .25, .5, .75, 1))
    df <- data.frame(row.names = sp,
      c1 = lat(), c2 = lat(),
      cnt = as.integer(stats::rpois(n, exp(1.5 + 0.8 * l_cnt))),
      prp = stats::plogis(l_prp + stats::rnorm(n, 0, 0.3)),
      bin = factor(ifelse(l_bin > stats::median(l_bin), "yes", "no")),
      ord = factor(cut(l_ord, ord_breaks, labels = c("L1", "L2", "L3", "L4"), include.lowest = TRUE),
                   levels = c("L1", "L2", "L3", "L4"), ordered = TRUE),
      cat3 = factor(cut(l_cat, stats::quantile(l_cat, c(0, 1/3, 2/3, 1)), labels = c("A", "B", "C"), include.lowest = TRUE)))
    df$prp <- pmin(pmax(df$prp, 1e-4), 1 - 1e-4)
    return(list(df = df, tree = tree, trait_types = c(prp = "proportion")))
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

# Build the cell: truth, tree, the user-level MCAR mask (seeded by seed + 1000), df_miss.
make_cell <- function(dgp, n, seed, miss_frac = 0.30) {
  d <- make_dgp(dgp, n, seed)
  truth <- d$df; tree <- d$tree
  set.seed(seed + 1000L)
  mask <- matrix(FALSE, nrow(truth), ncol(truth), dimnames = dimnames(truth))
  for (v in names(truth)) {
    obs <- which(!is.na(truth[[v]])); hide <- sample(obs, ceiling(miss_frac * length(obs)))
    mask[hide, v] <- TRUE
  }
  df_miss <- truth; for (v in names(truth)) df_miss[mask[, v], v] <- NA
  list(truth = truth, tree = tree, mask = mask, df_miss = df_miss,
       cont_traits = names(truth)[vapply(truth, is.numeric, logical(1))],
       trait_types = d$trait_types)
}

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


# ---- frequentist stack: Rphylopars on the continuous-family columns (count on log1p, proportion on
# logit, back-transformed), castor Mk (ML hidden-state prediction) on each discrete trait separately.
run_freq <- function(df_miss, truth, mask, tree, cont_traits, trait_types = NULL) {
  comp <- truth; comp[] <- NA
  # continuous family
  is_prop <- names(truth) %in% names(trait_types)[trait_types == "proportion"]
  is_cnt  <- vapply(truth, is.integer, logical(1))
  tf <- function(v, x) if (is_prop[match(v, names(truth))]) stats::qlogis(x) else if (is_cnt[match(v, names(truth))]) log1p(x) else x
  itf <- function(v, x) if (is_prop[match(v, names(truth))]) stats::plogis(x) else if (is_cnt[match(v, names(truth))]) pmax(expm1(x), 0) else x
  df4 <- df_miss[, cont_traits, drop = FALSE]
  for (v in cont_traits) df4[[v]] <- tf(v, as.numeric(df4[[v]]))
  df_in <- data.frame(species = rownames(df4), df4, stringsAsFactors = FALSE)
  fit <- Rphylopars::phylopars(df_in, tree = tree, model = "BM", phylo_correlated = TRUE, pheno_correlated = TRUE, REML = TRUE)
  rec <- fit$anc_recon[rownames(df4), cont_traits, drop = FALSE]
  for (v in cont_traits) { comp[[v]] <- df_miss[[v]]; comp[mask[, v], v] <- itf(v, rec[mask[, v], v]) }
  # discrete traits: castor Mk, equal rates for binary/categorical, stepwise (SUEDE) for ordinal
  for (v in setdiff(names(truth), cont_traits)) {
    f <- df_miss[[v]]; lev <- levels(f); tip <- as.integer(f)  # NA where masked
    tip_full <- tip[match(tree$tip.label, rownames(df_miss))]
    rm <- if (is.ordered(f)) "SUEDE" else "ER"
    h <- castor::hsp_mk_model(tree, tip_states = tip_full, Nstates = length(lev), rate_model = rm, Ntrials = 3, Nthreads = 1)
    pred <- max.col(h$likelihoods[seq_along(tree$tip.label), , drop = FALSE])
    pred <- pred[match(rownames(df_miss), tree$tip.label)]
    out <- as.character(f); out[mask[, v]] <- lev[pred[mask[, v]]]
    comp[[v]] <- factor(out, levels = lev, ordered = is.ordered(f))
  }
  list(completed = comp)
}

# ---- phylogeny + machine-learning hybrid (the approach of Gendre, Hauffe, Pimiento and Silvestro 2024, MEE,
# TDIP's missForest_phylo): missForest on the trait table plus phylogenetic eigenvectors. Eigenvectors are the
# principal coordinates of the cophenetic distance matrix (PVR-style), keeping axes that explain `variance_fraction`
# of the variance (TDIP's default idea), capped at `max_axes`.
phylo_eigenvectors <- function(tree, variance_fraction = 0.9, max_axes = 20L) {
  D <- ape::cophenetic.phylo(tree); n <- nrow(D)
  J <- diag(n) - matrix(1 / n, n, n); B <- -0.5 * J %*% (D^2) %*% J
  e <- eigen(B, symmetric = TRUE); keep <- e$values > 1e-8
  vals <- e$values[keep]; vecs <- e$vectors[, keep, drop = FALSE]
  k <- min(max_axes, which(cumsum(vals) / sum(vals) >= variance_fraction)[1])
  ev <- vecs[, seq_len(k), drop = FALSE] %*% diag(sqrt(vals[seq_len(k)]), k, k)
  rownames(ev) <- rownames(D); colnames(ev) <- paste0("pev", seq_len(k)); ev
}
run_mf_phylo <- function(df_miss, truth, mask, tree, variance_fraction = 0.9) {
  ev <- phylo_eigenvectors(tree, variance_fraction)[rownames(df_miss), , drop = FALSE]
  X <- cbind(df_miss, as.data.frame(ev))
  for (v in names(df_miss)) if (is.ordered(X[[v]])) X[[v]] <- factor(as.character(X[[v]]), levels = levels(df_miss[[v]]))
  imp <- missForest::missForest(X, maxiter = 10, ntree = 100, verbose = FALSE)$ximp
  comp <- truth; comp[] <- NA
  for (v in names(df_miss)) {
    x <- imp[[v]]
    if (is.factor(df_miss[[v]])) comp[[v]] <- factor(as.character(x), levels = levels(df_miss[[v]]), ordered = is.ordered(df_miss[[v]]))
    else if (is.integer(df_miss[[v]])) comp[[v]] <- as.integer(round(x)) else comp[[v]] <- x
  }
  list(completed = comp)
}
