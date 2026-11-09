args <- commandArgs(TRUE)
stopifnot(length(args) == 6L)
root <- normalizePath(args[1])
identity_file <- file.path(root, "SOURCE_SHA")
source_sha <- if (file.exists(identity_file)) readLines(identity_file, warn=FALSE)[1] else
  system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout=TRUE)[1]
stopifnot(identical(source_sha,"033649b93115c6d8191e791f06b1ddf52eb82d19"))
runner_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
runner_md5 <- unname(tools::md5sum(runner_file))
prepared_md5 <- unname(tools::md5sum(args[2]))
prepared <- readRDS(args[2])
out <- args[3]
fold <- as.integer(args[4])
seed <- as.integer(args[5])
epochs <- as.integer(args[6])
stopifnot(fold %in% 1:5, seed %in% c(0L, prepared$fit_seeds), epochs > 0)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
devtools::load_all(root, quiet=TRUE)
gnn <- seed != 0L
if (gnn) {
  torch::torch_set_num_threads(2L)
  torch::torch_set_num_interop_threads(1L)
}
d <- prepared$data
hidden <- prepared$fold == fold
traits <- data.frame(species=d$species, CTmax=d$CTmax)
traits$CTmax[hidden] <- NA_real_
stopifnot(all(is.na(traits$CTmax[hidden])), all(is.finite(d$CTmax[hidden])))
job <- sprintf("fold%d-seed%d-epochs%d", fold, seed, epochs)
path <- file.path(out, paste0(job, ".rds"))
if (file.exists(path)) stop("Result already exists: ", path)
warnings <- character()
started <- Sys.time()
record <- tryCatch(withCallingHandlers({
  fit <- pigauto::impute(traits, prepared$tree, species_col="species",
    log_transform=FALSE, covariates=prepared$covariates, gnn=gnn,
    seed=if (gnn) seed else 51000L, epochs=epochs, verbose=TRUE)
  # Completions retain input order; prediction bounds use internal row order.
  internal <- match(which(hidden), fit$data$input_row_order)
  stopifnot(!anyNA(internal))
  pred <- data.frame(row_n=d$row_n[hidden], species=d$species[hidden],
    observed=d$CTmax[hidden], predicted=fit$completed$CTmax[hidden],
    lower=fit$prediction$conformal_lower[internal, "CTmax"],
    upper=fit$prediction$conformal_upper[internal, "CTmax"])
  stopifnot(nrow(pred)==sum(hidden), all(is.finite(pred$predicted)),
            all(is.finite(pred$lower)), all(is.finite(pred$upper)),
            all(pred$lower <= pred$upper))
  list(status="ok", predictions=pred, gates=fit$fit$calibrated_gates,
       history=fit$fit$history, config=fit$fit$model_config)
}, warning=function(w) {
  warnings <<- c(warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
}), error=function(e) list(status="failed", error=conditionMessage(e)))
record$job <- job
record$fold <- fold
record$seed <- seed
record$epochs_requested <- epochs
record$n_expected <- sum(hidden)
record$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units="secs"))
record$warnings <- warnings
record$source_sha <- source_sha
record$runner_md5 <- runner_md5
record$prepared_md5 <- prepared_md5
record$split_seed <- prepared$split_seed
record$session <- sessionInfo()
saveRDS(record, path)
cat("RESULT", job, record$status, record$elapsed_seconds, "seconds\n")
if (record$status != "ok") quit(status=1)
