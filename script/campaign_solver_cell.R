# script/campaign_solver_cell.R
#
# Arc C (2026-09-19): which baseline solver closes the AVONET300 continuous-trait gap to raw
# Rphylopars, and does it hold on the simulated DGPs? GNN-off arms only, so a cell runs in seconds.
# Same DGPs, seeds and masks as the with/without-GNN campaign (script/campaign_gnn_off_lib.R).
#
# Usage: Rscript script/campaign_solver_cell.R --dgp avonet --n 300 --seed 1 --out results_solver/
#
# Arms (all pigauto with gnn = FALSE unless stated):
#   inhouse_pure    joint_solver = "inhouse", sigma_method = "single_pass" (the current default), pure baseline
#   inhouse_def     same, default safety machinery
#   fisher_pure     joint_solver = "inhouse", sigma_method = "fisher_ml" (opt-in on main, reached here
#                   through a namespace wrapper around fit_joint_solver(); diagnostic only), pure baseline
#   fisher_def      same, default safety machinery
#   rphylo_pure     joint_solver = "rphylopars", pure baseline
#   rphylo_def      joint_solver = "rphylopars", default safety machinery
#   bayes_pure      joint_solver = "inhouse", lambda_mode = "bayes", pure baseline
#   raw_rphylopars  Rphylopars::phylopars(model = "BM", pheno_correlated = TRUE, REML = TRUE) on the
#                   continuous columns only
#   raw_rphylo_nopheno  same with pheno_correlated = FALSE (is the phenotypic-variance term the difference?)
#   cont_only_pure  pigauto gnn = FALSE pure on the continuous columns only (mixed-type path removed)
#   floor           mean / mode

suppressPackageStartupMessages({ library(pigauto); library(ape) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) return(TRUE); args[i + 1L]
}
dgp <- get_arg("--dgp", "avonet"); n <- as.integer(get_arg("--n", 300L)); seed <- as.integer(get_arg("--seed", 1L))
out <- get_arg("--out", "results_solver"); miss_frac <- 0.30
arms <- strsplit(get_arg("--arms", "inhouse_pure,inhouse_def,fisher_pure,fisher_def,rphylo_pure,rphylo_def,bayes_pure,raw_rphylopars,raw_rphylo_nopheno,cont_only_pure,floor"), ",")[[1]]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
tag <- sprintf("%s_n%d_s%d", dgp, n, seed)
log_line <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "campaign_gnn_off_lib.R"))
cell_data <- make_cell(dgp, n, seed, miss_frac)
truth <- cell_data$truth; tree <- cell_data$tree; mask <- cell_data$mask
df_miss <- cell_data$df_miss; cont_traits <- cell_data$cont_traits
log_line("cell %s: n=%d traits=%d (continuous %d) masked=%d", tag, nrow(truth), ncol(truth), length(cont_traits), sum(mask))

# Diagnostic-only wrapper: force sigma_method on every fit_joint_solver() call inside pigauto.
with_sigma_method <- function(method, expr) {
  ns <- asNamespace("pigauto"); orig <- get("fit_joint_solver", ns)
  unlockBinding("fit_joint_solver", ns)
  assign("fit_joint_solver", function(L, tree, joint_solver = "inhouse", predict_method = "per_column",
                                      sigma_method = "single_pass", joint_refine_iter = 0L) {
    orig(L = L, tree = tree, joint_solver = joint_solver, predict_method = predict_method,
         sigma_method = method, joint_refine_iter = joint_refine_iter)
  }, envir = ns)
  on.exit({ assign("fit_joint_solver", orig, envir = ns); lockBinding("fit_joint_solver", ns) }, add = TRUE)
  force(expr)
}

pure <- list(safety_floor = FALSE, phylo_signal_gate = FALSE)
run_pig <- function(extra, sigma = NULL, traits = df_miss) {
  call <- function() do.call(pigauto::impute, c(list(traits = traits, tree = tree, gnn = FALSE, verbose = FALSE, seed = seed), extra))
  res <- if (is.null(sigma)) call() else with_sigma_method(sigma, call())
  comp <- truth; comp[] <- NA
  for (v in names(truth)) if (v %in% names(res$completed)) comp[[v]] <- res$completed[rownames(truth), v]
  list(completed = comp, lower = res$prediction$conformal_lower, upper = res$prediction$conformal_upper,
       path = res$fit$baseline$path, r_bm = res$fit$r_cal_bm, r_mean = res$fit$r_cal_mean)
}
run_raw <- function(pheno) {
  df4 <- df_miss[, cont_traits, drop = FALSE]
  df_in <- data.frame(species = rownames(df4), df4, stringsAsFactors = FALSE)
  fit <- Rphylopars::phylopars(df_in, tree = tree, model = "BM", phylo_correlated = TRUE,
                               pheno_correlated = pheno, REML = TRUE)
  comp <- truth; comp[] <- NA
  rec <- fit$anc_recon[rownames(df4), cont_traits, drop = FALSE]
  for (v in cont_traits) { comp[[v]] <- df_miss[[v]]; comp[mask[, v], v] <- rec[mask[, v], v] }
  list(completed = comp)
}
run_floor <- function() {
  comp <- df_miss
  for (v in names(comp)) { obs <- df_miss[[v]][!is.na(df_miss[[v]])]
    comp[mask[, v], v] <- if (is.numeric(comp[[v]])) mean(obs) else names(which.max(table(as.character(obs)))) }
  list(completed = comp)
}

results <- list(); walls <- list(); errors <- list(); paths <- list(); weights <- list()
for (arm in arms) {
  t0 <- Sys.time()
  r <- tryCatch(switch(arm,
    inhouse_pure = run_pig(c(pure, list(joint_solver = "inhouse"))),
    inhouse_def  = run_pig(list(joint_solver = "inhouse")),
    fisher_pure  = run_pig(c(pure, list(joint_solver = "inhouse")), sigma = "fisher_ml"),
    fisher_def   = run_pig(list(joint_solver = "inhouse"), sigma = "fisher_ml"),
    rphylo_pure  = run_pig(c(pure, list(joint_solver = "rphylopars"))),
    rphylo_def   = run_pig(list(joint_solver = "rphylopars")),
    bayes_pure   = run_pig(c(pure, list(joint_solver = "inhouse", lambda_mode = "bayes"))),
    raw_rphylopars = run_raw(TRUE),
    raw_rphylo_nopheno = run_raw(FALSE),
    cont_only_pure = run_pig(c(pure, list(joint_solver = "inhouse")), traits = df_miss[, cont_traits, drop = FALSE]),
    floor = run_floor(),
    stop("unknown arm ", arm)), error = function(e) e)
  walls[[arm]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error")) { errors[[arm]] <- conditionMessage(r); log_line("%s ERROR %s", arm, errors[[arm]]); next }
  results[[arm]] <- score_arm(arm, r$completed, r$lower, r$upper)
  paths[[arm]] <- r$path; weights[[arm]] <- list(r_bm = r$r_bm, r_mean = r$r_mean)
  log_line("%s done in %.1f s", arm, walls[[arm]])
}
tab <- do.call(rbind, results); rownames(tab) <- NULL
tab$dgp <- dgp; tab$n <- nrow(truth); tab$seed <- seed
saveRDS(list(tag = tag, dgp = dgp, n = nrow(truth), seed = seed, arms = arms, results = tab,
             walls = unlist(walls), errors = errors, paths = paths, weights = weights,
             pigauto_version = as.character(utils::packageVersion("pigauto")), host = Sys.info()[["nodename"]],
             time = Sys.time()), file.path(out, paste0(tag, ".rds")))
print(tab[tab$metric == "zRMSE", c("arm", "trait", "value", "coverage")], digits = 3)
if (length(errors)) { cat("ERRORS:\n"); print(errors) }
