# script/rubin_pig_crash_repro.R
#
# Reproduce one pig_post campaign failure ("leading principal minor ... is not positive", CHOLMOD) and record the
# call stack at the error, so the fix can target the right site in R/mi_posterior.R. Same dataset, seed offset (909)
# and call as the campaign (rubin_cell.R).
#
#   Rscript script/rubin_pig_crash_repro.R --n 100 --lambda 0.3 --rho 0 --seed 5 --out <dir>

RNGkind("L'Ecuyer-CMRG")
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1L] }
n <- as.integer(get_arg("--n")); lambda <- as.numeric(get_arg("--lambda")); rho <- as.numeric(get_arg("--rho"))
seed <- as.integer(get_arg("--seed")); out <- get_arg("--out")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
f_out <- file.path(out, sprintf("crash_n%d_l%s_r%s_s%d.rds", n, format(lambda), format(rho), seed))

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
suppressMessages({ for (f in c("campaign_gnn_off_lib.R", "rubin_lib.R", "rubin_freq.R", "rubin_pigauto.R")) source(file.path(here, f)) })
cell <- make_cell("types_mixed", n, seed, miss_frac = 0.30, miss = "mcar", lambda = lambda, rho = rho,
                  thresholds = "fixed", driver = TRUE)
calls <- NULL; frames_info <- NULL
set.seed(seed + 909L)
res <- tryCatch(withCallingHandlers(mi_pig_post(cell, 20, seed = seed + 909L), error = function(e) {
  sc <- sys.calls()
  calls <<- vapply(sc, function(cl) paste(deparse(cl, width.cutoff = 200L)[1], collapse = ""), "")
  # Sigma_E / Sigma_P at the failing frame, when visible
  for (i in rev(seq_along(sc))) {
    fr <- sys.frame(i)
    if (exists("Sigma_E", envir = fr, inherits = FALSE)) {
      frames_info <<- list(frame = calls[i], Sigma_E = get("Sigma_E", envir = fr),
                           SE_new = if (exists("SE_new", envir = fr, inherits = FALSE)) get("SE_new", envir = fr) else NULL)
      break
    }
  }
}), error = function(e) e)
saveRDS(list(n = n, lambda = lambda, rho = rho, seed = seed, failed = inherits(res, "error"),
             message = if (inherits(res, "error")) conditionMessage(res) else NA, calls = calls, frame = frames_info), f_out)
cat(basename(f_out), if (inherits(res, "error")) "FAILED" else "ok", "\n")
