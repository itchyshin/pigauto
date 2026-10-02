# script/rubin_pig_gate.R
#
# Gate checks for the pig_post arm (.unlazy/pigauto-rubin-arm/GATES.md). Each mode prints its success
# marker only after every assertion passes, otherwise stops.
#
#   Rscript script/rubin_pig_gate.R arms     <smoke dir>
#   Rscript script/rubin_pig_gate.R identity <smoke dir> [pool dir]
#   Rscript script/rubin_pig_gate.R repro    <smoke dir> [pool dir]
#   Rscript script/rubin_pig_gate.R prov     <smoke dir> <expected pigauto sha>
#   Rscript script/rubin_pig_gate.R prerun   <pre-run dir>

`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; dir <- args[2]
pool <- if (length(args) >= 3 && mode %in% c("identity", "repro")) args[3] else path.expand("~/pigauto_rubin_pool/freq")
files <- list.files(dir, "^rubin_.*\\.rds$", full.names = TRUE, recursive = TRUE)
if (length(files) < 1L) stop("no rubin rds under ", dir)

stored_twin <- function(x) {
  tag <- sub("_smoke$", "", x$tag)
  f <- Sys.glob(file.path(pool, "*", paste0(tag, ".rds")))
  if (length(f) != 1L) stop(length(f), " stored campaign rds for ", tag, " under ", pool, " (need exactly 1)")
  y <- readRDS(f); cat("stored twin:", f, "git_hash", y$git_hash, "\n"); y
}

one_file <- function() { if (length(files) != 1L) stop(length(files), " rds under ", dir, " (need exactly 1)"); readRDS(files) }

if (identical(mode, "arms")) {
  x <- one_file(); est <- x$estimands
  want <- c("freqA", "freqB", "bace", "bace_chain", "bace_resid", "pig_post")
  got <- intersect(want, unique(est$arm))
  ok <- vapply(got, function(a) {
    e <- est[est$arm == a, ]
    nrow(e) == 2L && setequal(e$estimand, c("slope", "cor")) && all(is.finite(e$estimate)) && all(is.finite(e$se))
  }, logical(1))
  if (length(x$errors)) stop("errors recorded: ", paste(names(x$errors), collapse = ", "))
  if (!all(want %in% got) || !all(ok)) stop("arms missing or non-finite: ", paste(setdiff(want, got[ok]), collapse = ", "))
  cat("ARMS_OK", sum(ok), "\n")
} else if (identical(mode, "identity")) {
  x <- one_file(); y <- stored_twin(x)
  # Cross-platform (Mac vs the cluster that wrote the stored file): numeric traits agree to floating-point
  # rounding (measured 1.7e-15), so they are compared at 1e-12; the mask and discrete traits must be identical.
  if (!identical(x$mask, y$mask)) stop("mask differs from the stored campaign dataset")
  if (!identical(names(x$truth), names(y$truth)) || !identical(rownames(x$truth), rownames(y$truth))) stop("truth layout differs")
  for (v in names(x$truth)) {
    a <- x$truth[[v]]; b <- y$truth[[v]]
    same <- if (is.numeric(a)) isTRUE(max(abs(a - b)) <= 1e-12) else identical(a, b)
    if (!same) stop("truth column ", v, " differs from the stored campaign dataset")
  }
  cat("IDENTITY_OK\n")
} else if (identical(mode, "repro")) {
  x <- one_file(); y <- stored_twin(x)
  key <- c("estimate", "se", "lower", "upper")
  for (a in c("complete", "freqA", "freqB")) {
    ex <- x$estimands[x$estimands$arm == a, ]; ey <- y$estimands[y$estimands$arm == a, ]
    ex <- ex[order(ex$estimand), key]; ey <- ey[order(ey$estimand), key]
    # relative 1e-6: cross-platform rounding measured up to 2.4e-9 (freqA); a disturbed arm moves estimates ~1e-2
    if (nrow(ex) != 2L || !isTRUE(all.equal(unname(as.matrix(ex)), unname(as.matrix(ey)), tolerance = 1e-6)))
      stop(a, " does not reproduce the stored campaign row")
  }
  cat("REPRO_OK\n")
} else if (identical(mode, "prov")) {
  want_sha <- args[3]
  if (is.na(want_sha) || !grepl("^[0-9a-f]{40}$", want_sha)) stop("prov needs the expected 40-hex pigauto sha")
  for (f in files) {
    d <- readRDS(f)$diag$pig_post
    if (is.null(d)) stop("no pig_post diag in ", basename(f))
    if (is.na(d$pigauto_sha) || !identical(d$pigauto_sha, want_sha)) stop("pigauto sha ", d$pigauto_sha, " != ", want_sha, " in ", basename(f))
    if (!identical(d$mi_workflow, "pigauto_posterior_mi_v1")) stop("not proper posterior MI in ", basename(f))
    if (!isTRUE(d$converged)) stop("not converged (R-hat ", d$rhat_max, ", ESS ", d$ess_min, ") in ", basename(f))
    if (length(d$logged)) stop("log-transformed traits in ", basename(f))
  }
  cat("PROV_OK\n")
} else if (identical(mode, "prerun")) {
  n_ok <- 0L; fails <- character(0); rows <- list()
  for (f in files) {
    x <- readRDS(f); e <- x$estimands[x$estimands$arm == "pig_post", ]; d <- x$diag$pig_post
    ok <- nrow(e) == 2L && all(is.finite(e$estimate))
    if (ok) n_ok <- n_ok + 1L else fails <- c(fails, sprintf("%s [%s]", basename(f), paste(unlist(x$errors), collapse = "; ")))
    rows[[f]] <- data.frame(n = x$n, lambda = x$lambda, rho = x$rho, seed = x$seed, ok = ok,
                            converged = isTRUE(d$converged), n_ext = d$n_extensions %||% NA, wall_s = x$walls[["pig_post"]] %||% NA)
  }
  tab <- do.call(rbind, rows)
  print(aggregate(cbind(converged, n_ext, wall_s) ~ n, tab, mean))
  cat("files:", length(files), " pig_post failures:", length(fails), "\n"); if (length(fails)) cat(fails, sep = "\n")
  cat("PRERUN", n_ok, "\n")
} else stop("unknown mode ", mode)
