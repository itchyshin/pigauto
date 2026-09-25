#!/usr/bin/env Rscript
# Usage: Rscript 09_gate_smoke.R
# Fast gate: builds a 30-tip AVONET input from pigauto's bundled
# avonet300/tree300, runs 01_run_masked_confirmation.R (epochs = 5) and
# 02_summarise_masked_confirmation.R into a tempdir for both arms, and
# prints SMOKE_OK only when all four artefacts (mask_receipt.rds,
# split.rds, mondrian.rds, summary.rds) exist and the summary metrics table
# has rows. Intended as a <1 min sanity check, not evidence for the
# pre-registered study.
# NOTE: this parent process only orchestrates child `Rscript` calls (01/02)
# via system2() -- it must NOT initialize torch itself. Loading libtorch here
# starts its internal thread pool, and forking via system2() afterwards can
# corrupt the child process's startup on macOS (observed as spurious R
# errors unrelated to any code in this file). Thread caps are still applied
# in-process by 01_run_masked_confirmation.R itself.
suppressPackageStartupMessages(library(ape))

`%||%` <- function(a, b) if (is.null(a)) b else a
# Deliberately relative to the current working directory, matching every
# other script in this directory (run from the repo root). An absolute path
# here is avoided on purpose -- passing one to the child `Rscript` calls
# below was observed to break child startup in this environment.
script_dir <- "script/mondrian_confirmation"
if (!file.exists(file.path(script_dir, "01_run_masked_confirmation.R"))) {
  stop("run this script from the pigauto repo root", call. = FALSE)
}

e <- new.env(parent = emptyenv())
utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300",   package = "pigauto", envir = e)
df <- e$avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL
tree <- e$tree300

set.seed(20260818L)
keep <- sort(sample(tree$tip.label, 30L))
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
df <- df[tree$tip.label, , drop = FALSE]
stopifnot(identical(rownames(df), tree$tip.label))

input_file <- tempfile(fileext = ".rds")
saveRDS(list(data = df, tree = tree, dataset = "avonet300_smoke"), input_file)

run_arm <- function(arm) {
  run_dir <- tempfile(paste0("gate_smoke_", arm, "_"))
  summary_file <- tempfile(fileext = ".rds")
  status <- system2("Rscript", shQuote(c(
    file.path(script_dir, "01_run_masked_confirmation.R"),
    input_file, run_dir, arm, "20260818", "5", "split,mondrian"
  )))
  if (!identical(status, 0L)) stop("01_run_masked_confirmation.R failed for arm ", arm, call. = FALSE)
  status <- system2("Rscript", shQuote(c(
    file.path(script_dir, "02_summarise_masked_confirmation.R"),
    run_dir, summary_file
  )))
  if (!identical(status, 0L)) stop("02_summarise_masked_confirmation.R failed for arm ", arm, call. = FALSE)

  artefacts <- c(
    file.path(run_dir, "mask_receipt.rds"),
    file.path(run_dir, "split.rds"),
    file.path(run_dir, "mondrian.rds"),
    summary_file
  )
  ok <- all(file.exists(artefacts))
  if (ok) {
    summ <- readRDS(summary_file)
    ok <- ok && !is.null(summ$metrics) && nrow(summ$metrics) > 0L
  }
  ok
}

ok_mcar <- run_arm("mcar")
ok_structured <- run_arm("structured")

if (isTRUE(ok_mcar) && isTRUE(ok_structured)) {
  cat("SMOKE_OK (mcar)\n")
  cat("SMOKE_OK (structured)\n")
} else {
  stop(sprintf("smoke failed: mcar=%s structured=%s", ok_mcar, ok_structured), call. = FALSE)
}
