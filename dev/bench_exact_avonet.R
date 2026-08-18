# Slice 1b real-data gate: does cross-trait refinement (predict_method)
# move pigauto's continuous traits toward the rphylopars numbers WITHOUT a
# dependency? Same 5 masks / metric as script/bench_external_comparators.md.
#
# Reference on these masks (mean z-RMSE over 5 seeds):
#   pigauto inhouse, no refinement : Mass 1.594  Beak 0.912  Tarsus 1.220  Wing 0.688
#   pigauto joint_solver=rphylopars: Mass 1.295  Beak 0.602  Tarsus 0.873  Wing 0.449
#   raw rphylopars                 : Mass 1.360  Beak 0.445  Tarsus 0.639  Wing 0.409
#
# PASS = refinement lands at or near the rphylopars-solver row, dependency-free.

suppressPackageStartupMessages({ library(pigauto); library(ape) })
try(torch::torch_set_num_threads(1L), silent = TRUE)
t0 <- proc.time()[["elapsed"]]
say <- function(...) { cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
                           ..., "\n", sep = ""); flush.console() }

stopifnot("predict_method" %in% names(formals(pigauto::impute)))
MISS <- 0.30; CONT <- c("Mass","Beak.Length_Culmen","Tarsus.Length","Wing.Length")
ITERS <- c("per_column", "exact")

e <- new.env(); utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300", package = "pigauto", envir = e)
df <- e$avonet300; tree <- e$tree300
rownames(df) <- df$Species_Key; df$Species_Key <- NULL

rows <- list()
for (rep_id in 1:5) {
  seed <- 2026L + rep_id
  set.seed(seed)
  mask <- matrix(FALSE, nrow(df), ncol(df), dimnames=list(rownames(df), names(df)))
  for (v in names(df)) { obs <- which(!is.na(df[[v]]))
    mask[sample(obs, ceiling(length(obs) * MISS)), v] <- TRUE }
  dm <- df; for (v in names(df)) dm[[v]][mask[, v]] <- NA
  for (it in ITERS) {
    say(sprintf("rep %d seed %d predict_method=%s", rep_id, seed, it))
    res <- tryCatch(pigauto::impute(dm, tree, seed = seed, verbose = FALSE,
                                    n_imputations = 1L, predict_method = it),
                    error = function(er) er)
    if (inherits(res, "error")) { say("FAIL: ", conditionMessage(res)); next }
    for (v in CONT) {
      idx <- which(mask[, v]); s <- stats::sd(df[[v]][!mask[, v]], na.rm = TRUE)
      rows[[length(rows)+1L]] <- data.frame(rep = rep_id, iter = it, trait = v,
        rmse_z = sqrt(mean(((as.numeric(res$completed[[v]][idx]) - df[[v]][idx])/s)^2, na.rm=TRUE)))
    }
  }
}
d <- do.call(rbind, rows)
agg <- do.call(rbind, lapply(split(d, list(d$trait, d$iter), drop=TRUE), function(g)
  data.frame(trait=g$trait[1], iter=g$iter[1], n=nrow(g), rmse_z=mean(g$rmse_z),
             mcse=stats::sd(g$rmse_z)/sqrt(nrow(g)))))
rownames(agg) <- NULL
agg <- agg[order(agg$trait, agg$iter), ]
say("results:"); print(agg, digits = 3, row.names = FALSE)
saveRDS(list(raw = d, agg = agg), "dev/bench_refine_avonet.rds")
say("wrote dev/bench_refine_avonet.rds")
