# Option-C tests for the discrete-lambda / safety-gate decision (Fisher review, 2026-10-05).
# Scores with the four-arm harness's own run_arms()/score_arm() (sourced from the patched Totoro copy).
#   --kind alldisc : types_mixed cell with only bin, ord, cat3 kept (no continuous trait -> label-propagation path)
#   --kind bien    : BIEN plants (4,745 species, 5 continuous traits), MCAR 30% mask from --seed
# Arms: gnn_off (default), gnn_off_nogate, gnn_off_pure, floor. Discrete lambda comes from R_PROFILE_USER.
args <- commandArgs(TRUE); ga <- function(k, d = NULL) { i <- match(k, args); if (is.na(i)) d else args[i + 1] }
kind <- ga("--kind"); seed <- as.integer(ga("--seed")); n <- as.integer(ga("--n", "300"))
lambda <- as.numeric(ga("--lambda", "0.3")); rho <- as.numeric(ga("--rho", "0")); out <- ga("--out"); arms <- strsplit(ga("--arms"), ",")[[1]]
RNGkind("L'Ecuyer-CMRG")
here <- "~/pigauto_sim/dlam/script"
suppressMessages(source(file.path(here, "campaign_gnn_off_lib.R")))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (kind == "alldisc") {
  cell <- make_cell("types_mixed", n, seed, miss_frac = 0.3, miss = "mcar", lambda = lambda, rho = rho,
                    thresholds = "fixed", driver = TRUE)
  keep <- c("bin", "ord", "cat3")
  cell$truth <- cell$truth[, keep, drop = FALSE]; cell$mask <- cell$mask[, keep, drop = FALSE]
  cell$df_miss <- cell$df_miss[, keep, drop = FALSE]; cell$cont_traits <- character(0)
  if (!is.null(cell$trait_types)) cell$trait_types <- cell$trait_types[intersect(names(cell$trait_types), keep)]
  tag <- sprintf("alldisc_l%s_r%s_n%d_s%d", format(lambda), format(rho), n, seed)
} else if (kind == "bien") {
  # same construction as script/bench_bien.R (wide table from the per-trait mean tables; tree tips with spaces)
  tm <- readRDS("~/pigauto_sim/dlam/bien/bien_trait_means.rds"); tree <- readRDS("~/pigauto_sim/dlam/bien/bien_tree.rds"); if (!inherits(tree, "phylo")) tree <- tree$scenario.3
  sp <- Reduce(union, lapply(tm, function(d) if (is.null(d)) character(0) else d$species))
  wide <- data.frame(species = sp)
  for (nm in names(tm)) { d <- tm[[nm]]; vc <- intersect(c("mean_value", "trait_value", "value"), names(d))[1]
    wide[[nm]] <- if (is.null(d) || is.na(vc)) NA_real_ else suppressWarnings(as.numeric(d[[vc]][match(wide$species, d$species)])) }
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  keep <- intersect(wide$species, tree$tip.label)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  wide <- wide[match(tree$tip.label, wide$species), , drop = FALSE]
  df <- wide[, setdiff(names(wide), "species"), drop = FALSE]; rownames(df) <- wide$species
  df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]; tree <- ape::keep.tip(tree, rownames(df))
  set.seed(seed + 1000L)
  mask <- matrix(FALSE, nrow(df), ncol(df), dimnames = dimnames(df))
  for (v in names(df)) { obs <- which(!is.na(df[[v]])); mask[sample(obs, ceiling(0.3 * length(obs))), v] <- TRUE }
  dm <- df; for (v in names(df)) dm[mask[, v], v] <- NA
  cell <- list(truth = df, tree = tree, mask = mask, df_miss = dm, cont_traits = names(df), trait_types = NULL, covs = NULL)
  tag <- sprintf("bien_s%d", seed)
} else stop("unknown --kind")
t0 <- proc.time()[["elapsed"]]
res <- run_arms(cell, arms, list(seed = seed))
saveRDS(list(results = res$results %||% res, kind = kind, seed = seed, n = n, lambda = lambda, rho = rho,
             wall = proc.time()[["elapsed"]] - t0, option = getOption("pigauto.discrete_lambda"),
             pigauto = as.character(utils::packageVersion("pigauto"))), file.path(out, paste0(tag, ".rds")))
cat("done", tag, round(proc.time()[["elapsed"]] - t0), "s\n")
