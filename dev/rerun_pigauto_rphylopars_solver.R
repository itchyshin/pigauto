# Confirmation for the joint_solver switch: re-run ONLY the pigauto arm of
# script/bench_external_comparators.R (same masks, seeds 2027-2031; same z-RMSE
# metric) with joint_solver = "rphylopars", side by side with the committed
# 5-seed results for pigauto-inhouse / raw rphylopars / phylolm.
# Acceptance (docs/dev-log/2026-08-16-continuous-gap-diagnosis.md): recovery
# toward the arm-C level (0.37-0.87 z-RMSE), i.e. most of the solver gap closes.

suppressPackageStartupMessages({ library(pigauto); library(ape) })
try(torch::torch_set_num_threads(1L), silent = TRUE)
t0 <- proc.time()[["elapsed"]]
say <- function(...) { cat(sprintf("[%7.1fs] ", proc.time()[["elapsed"]] - t0),
                           ..., "\n", sep = ""); flush.console() }

MISS_FRAC <- 0.30; N_REPS <- 5L; MASK_SEED_BASE <- 2026L
CONT <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

e <- new.env(parent = emptyenv())
utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300",   package = "pigauto", envir = e)
df <- e$avonet300; tree <- e$tree300
rownames(df) <- df$Species_Key; df$Species_Key <- NULL
stopifnot(all(rownames(df) == tree$tip.label),
          "joint_solver" %in% names(formals(pigauto::impute)))

rows <- list()
for (rep_id in seq_len(N_REPS)) {
  seed <- MASK_SEED_BASE + rep_id
  say(sprintf("rep %d/%d (mask seed %d)", rep_id, N_REPS, seed))
  set.seed(seed)
  mask <- matrix(FALSE, nrow(df), ncol(df), dimnames = list(rownames(df), names(df)))
  for (v in names(df)) {
    obs_idx <- which(!is.na(df[[v]]))
    mask[sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC)), v] <- TRUE
  }
  df_miss <- df
  for (v in names(df)) df_miss[[v]][mask[, v]] <- NA

  res <- tryCatch(
    pigauto::impute(df_miss, tree, seed = seed, verbose = FALSE,
                    n_imputations = 1L, joint_solver = "rphylopars"),
    error = function(er) er)
  if (inherits(res, "error")) {
    say("FAIL: ", conditionMessage(res)); next
  }
  for (v in CONT) {
    idx <- which(mask[, v])
    obs <- df[[v]][!mask[, v]]
    s <- stats::sd(obs, na.rm = TRUE)
    rz <- sqrt(mean(((as.numeric(res$completed[[v]][idx]) - df[[v]][idx]) / s)^2,
                    na.rm = TRUE))
    rows[[length(rows) + 1L]] <- data.frame(rep = rep_id, trait = v, rmse_z = rz)
  }
}
d <- do.call(rbind, rows)
agg <- do.call(rbind, lapply(split(d, d$trait), function(g)
  data.frame(trait = g$trait[1], n = nrow(g), rmse_z = mean(g$rmse_z),
             mcse = stats::sd(g$rmse_z) / sqrt(nrow(g)))))
rownames(agg) <- NULL
say("pigauto with joint_solver='rphylopars' (mean +/- MCSE over 5 masks):")
print(agg, digits = 3)
saveRDS(list(raw = d, agg = agg), "dev/rerun_pigauto_rphylopars_solver.rds")
say("wrote dev/rerun_pigauto_rphylopars_solver.rds")
say("Reference (committed bench, same masks): pigauto-inhouse Mass 1.594 Beak 0.912")
say("Tarsus 1.220 Wing 0.688 ; raw rphylopars Mass 1.360 Beak 0.445 Tarsus 0.639 Wing 0.409")
