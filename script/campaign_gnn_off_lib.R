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
       cont_traits = names(truth)[vapply(truth, is.numeric, logical(1))])
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

