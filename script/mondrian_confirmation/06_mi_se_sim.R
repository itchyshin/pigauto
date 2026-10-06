#!/usr/bin/env Rscript
# Usage: Rscript 06_mi_se_sim.R <rep> <outdir>
#
# One paired rep of the Mondrian-vs-split multi_impute() Rubin-SE
# comparison pre-registered in useful/mondrian-mi-se-justification.md.
# Writes results/<outdir>/rep_<rep>.rds containing BOTH methods (split and
# mondrian) for that one (tree, mask, truth) draw, plus a complete-data
# reference arm. Do not run more than one rep locally -- see that memo's
# "Simulation design" section for the wall-time estimate and the DRAC
# job-array recommendation for the full 500-rep campaign.
#
# DGP: two correlated (rho = 0.7) BM traits on an n-tip tree, reusing the
# structure (cophenetic-scaled V, MAR_phylo clade mechanism) of
# ~/pigauto_regime_map/mech_cell.R, simplified from that script's four
# traits to the two needed here: x (subject to MAR_phylo missingness) and
# y (always fully observed). Downstream model:
# gls(y ~ x, correlation = corBrownian(tree)), pooled via pool_mi().
#
# n / epochs default to the pre-registered n = 1000, epochs = 500 and are
# overridable via env vars MI_SE_SIM_N / MI_SE_SIM_EPOCHS for a cheap local
# smoke run (e.g. n = 200, epochs = 20) without changing the two-argument
# CLI contract above.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("expected: rep outdir", call. = FALSE)
rep_i   <- as.integer(args[[1L]])
outdir  <- args[[2L]]
if (!is.finite(rep_i)) stop("rep must be an integer", call. = FALSE)

n      <- as.integer(Sys.getenv("MI_SE_SIM_N", "1000"))
epochs <- as.integer(Sys.getenv("MI_SE_SIM_EPOCHS", "500"))
m      <- as.integer(Sys.getenv("MI_SE_SIM_M", "20"))
m_miss <- 0.3
true_rho <- 0.7
true_beta <- true_rho  # y, x both unit-variance BM traits -> slope == correlation

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressMessages({ library(pigauto); library(ape); library(nlme) })
if (requireNamespace("torch", quietly = TRUE)) {
  try(torch::torch_set_num_threads(1L), silent = TRUE)
  try(torch::torch_set_num_interop_threads(1L), silent = TRUE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(outdir, sprintf("rep_%d.rds", rep_i))
if (file.exists(out_f)) { cat("SKIP", out_f, "\n"); quit(save = "no") }

seed <- 20260923L + rep_i
set.seed(seed)

# ---- DGP: 2 correlated BM traits, x/y (mech_cell.R structure, p = 2) ------
tree <- ape::rtree(n)
V <- ape::vcv(tree); V <- V / max(V)
L <- chol(V + 1e-8 * diag(n))
Sig <- matrix(true_rho, 2, 2); diag(Sig) <- 1
Z <- t(L) %*% matrix(rnorm(n * 2), n, 2) %*% chol(Sig)
truth <- data.frame(row.names = tree$tip.label, x = Z[, 1], y = Z[, 2])

# ---- MAR_phylo mechanism on x only (mech_cell.R's clade construction) -----
nodes <- (n + 2L):(n + tree$Nnode)
sizes <- vapply(nodes, function(nd) length(
  ape::extract.clade(tree, nd)$tip.label), integer(1))
cand <- nodes[sizes >= floor(0.15 * n) & sizes <= ceiling(0.35 * n)]
picked <- if (length(cand) >= 2L) sample(cand, 2L) else
          nodes[order(abs(sizes - 0.25 * n))[1:2]]
in_clade <- rep(FALSE, n); names(in_clade) <- tree$tip.label
for (nd in picked) in_clade[ape::extract.clade(tree, nd)$tip.label] <- TRUE

pvec <- ifelse(in_clade[tree$tip.label], 7, 1)
pvec <- pvec * (m_miss * n) / sum(pvec)
pvec <- pmin(pvec, 0.95)
mask_x <- runif(n) < pvec
if (sum(!mask_x) < 20L) {
  keep <- sample(which(mask_x), sum(mask_x) - (n - 20L))
  mask_x[keep] <- FALSE
}

df <- truth
df$x[mask_x] <- NA
species <- tree$tip.label

# ---- reference arm: complete-data gls on the (unobserved-in-practice) --
# truth, establishing the target sampling SD with no imputation
# uncertainty at all.
ref_dat <- truth
ref_dat$species <- species
ref_fit <- tryCatch(
  nlme::gls(y ~ x, correlation = ape::corBrownian(phy = tree, form = ~species),
            data = ref_dat, method = "ML"),
  error = function(e) e
)
ref_out <- if (inherits(ref_fit, "error")) {
  list(estimate = NA_real_, std.error = NA_real_, failed = TRUE,
       error = conditionMessage(ref_fit))
} else {
  co <- summary(ref_fit)$tTable
  list(estimate = unname(co["x", "Value"]),
       std.error = unname(co["x", "Std.Error"]), failed = FALSE)
}

# ---- helper: one method's multi_impute() -> per-draw gls -> pool_mi() ----
run_arm <- function(conformal_method) {
  t0 <- proc.time()[["elapsed"]]
  mi <- tryCatch(
    suppressMessages(pigauto::multi_impute(
      df, tree, m = m, draws_method = "conformal",
      conformal_method = conformal_method,
      epochs = epochs, verbose = FALSE, seed = seed, gnn = TRUE
    )),
    error = function(e) e
  )
  wall <- proc.time()[["elapsed"]] - t0
  if (inherits(mi, "error")) {
    return(list(failed = TRUE, error = conditionMessage(mi), wall_s = wall))
  }

  fits <- lapply(mi$datasets, function(dat) {
    dat$species <- rownames(dat)
    tryCatch(
      nlme::gls(y ~ x, correlation = ape::corBrownian(phy = tree, form = ~species),
                data = dat, method = "ML"),
      error = function(e) NULL
    )
  })
  ok <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2L) {
    return(list(failed = TRUE, error = "fewer than 2 successful per-draw fits",
                wall_s = wall))
  }
  pooled <- tryCatch(pigauto::pool_mi(fits[ok]), error = function(e) e)
  if (inherits(pooled, "error")) {
    return(list(failed = TRUE, error = conditionMessage(pooled), wall_s = wall))
  }
  row_x <- pooled[pooled$term == "x", , drop = FALSE]

  # ---- per-missing-x-cell draw PIT / coverage, by Mondrian stratum ------
  miss_sp <- species[mask_x]
  draws_x <- vapply(mi$datasets, function(dat) as.numeric(dat[miss_sp, "x"]),
                    numeric(length(miss_sp)))
  if (is.null(dim(draws_x))) draws_x <- matrix(draws_x, nrow = 1L)
  rownames(draws_x) <- miss_sp
  truth_x <- truth[miss_sp, "x"]

  pit <- vapply(seq_along(miss_sp), function(i)
    mean(draws_x[i, ] <= truth_x[i]), numeric(1))
  lo <- apply(draws_x, 1L, stats::quantile, probs = 0.025, na.rm = TRUE)
  hi <- apply(draws_x, 1L, stats::quantile, probs = 0.975, na.rm = TRUE)
  covered <- truth_x >= lo & truth_x <= hi

  stratum <- rep(NA_character_, length(miss_sp))
  if (identical(conformal_method, "mondrian") &&
      identical(mi$fit$conformal_method, "mondrian") &&
      !is.null(mi$fit$conformal_mondrian[["x"]]) &&
      !isTRUE(mi$fit$conformal_mondrian[["x"]]$fallback)) {
    D_sq <- mi$fit$graph$D_sq
    tm_x <- mi$fit$trait_map[[which(vapply(mi$fit$trait_map, `[[`,
                                           character(1), "name") == "x")]]
    n_sp <- length(mi$fit$species_names)
    cs <- pigauto:::mondrian_cell_scores(mi$fit, D_sq, mi$fit$trait_map, n_sp)
    if (!is.null(cs)) {
      obs_mask_train <- !is.na(mi$fit$X_scaled)
      hold <- c(mi$fit$splits$val_idx, mi$fit$splits$test_idx)
      if (length(hold)) obs_mask_train[hold] <- FALSE
      obs_idx <- which(obs_mask_train[, tm_x$latent_cols[1L]])
      target_idx <- match(miss_sp, mi$fit$species_names)
      if (length(obs_idx) > 0L) {
        locality <- pigauto:::mondrian_locality(D_sq, obs_idx, target_idx, k = 5L)
        mo <- mi$fit$conformal_mondrian[["x"]]
        stratum <- ifelse(is.finite(locality),
                          ifelse(locality > mo$threshold, "far", "near"),
                          NA_character_)
      }
    }
  }

  list(
    failed = FALSE, wall_s = wall,
    estimate = row_x$estimate, std.error = row_x$std.error,
    df = row_x$df, fmi = row_x$fmi, riv = row_x$riv,
    cell = data.frame(species = miss_sp, pit = pit, covered = covered,
                      stratum = stratum, stringsAsFactors = FALSE)
  )
}

split_out    <- run_arm("split")
mondrian_out <- run_arm("mondrian")

result <- list(
  rep = rep_i, seed = seed, n = n, epochs = epochs, m = m,
  true_beta = true_beta, m_miss_target = m_miss,
  n_missing_x = sum(mask_x),
  reference = ref_out,
  split = split_out,
  mondrian = mondrian_out
)
saveRDS(result, out_f)

fail_note <- function(x) if (isTRUE(x$failed)) paste0(" FAILED(", x$error, ")") else ""
cat(sprintf(
  "OK rep=%d split_est=%.4f%s mondrian_est=%.4f%s wall_split=%.1fs wall_mondrian=%.1fs\n",
  rep_i,
  if (isTRUE(split_out$failed)) NA_real_ else split_out$estimate, fail_note(split_out),
  if (isTRUE(mondrian_out$failed)) NA_real_ else mondrian_out$estimate, fail_note(mondrian_out),
  if (isTRUE(split_out$failed)) NA_real_ else split_out$wall_s,
  if (isTRUE(mondrian_out$failed)) NA_real_ else mondrian_out$wall_s
))
