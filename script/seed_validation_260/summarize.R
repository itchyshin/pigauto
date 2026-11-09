args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(3L,4L))
partial <- length(args)==4L && identical(args[4],"--partial")
stopifnot(length(args)==3L || partial)
script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
expected_runner_md5 <- unname(tools::md5sum(file.path(dirname(script_file),"run.R")))
prepared <- readRDS(args[1])
prepared_md5 <- unname(tools::md5sum(args[1]))
results <- args[2]
out <- args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
jobs <- expand.grid(fold=1:5, seed=c(0L,prepared$fit_seeds))
records <- vector("list",nrow(jobs))
attempts <- jobs
attempts$status <- "missing"
attempts$error <- ""
attempts$elapsed_seconds <- NA_real_
attempts$n_expected <- vapply(jobs$fold,function(f)sum(prepared$fold==f),integer(1))
for(i in seq_len(nrow(jobs))) {
  p <- file.path(results,sprintf("fold%d-seed%d-epochs2000.rds",jobs$fold[i],jobs$seed[i]))
  if(!file.exists(p)) next
  r <- tryCatch({
    r <- readRDS(p)
    stopifnot(identical(as.integer(r$fold),as.integer(jobs$fold[i])),
      identical(as.integer(r$seed),as.integer(jobs$seed[i])),
      identical(as.integer(r$n_expected),as.integer(attempts$n_expected[i])),
      identical(r$epochs_requested,2000L),identical(r$split_seed,prepared$split_seed),
      identical(r$source_sha,"033649b93115c6d8191e791f06b1ddf52eb82d19"),
      identical(r$prepared_md5,prepared_md5),length(r$runner_md5)==1L,
      identical(r$runner_md5,expected_runner_md5),
      length(r$status)==1L,r$status %in% c("ok","failed"))
    if(r$status=="ok") {
    expected <- prepared$data[prepared$fold==r$fold,c("row_n","species","CTmax")]
    pred <- r$predictions
    stopifnot(nrow(pred)==nrow(expected),!anyDuplicated(pred$row_n),
      identical(pred$row_n,expected$row_n),identical(pred$species,expected$species),
      identical(pred$observed,expected$CTmax),all(is.finite(pred$predicted)))
    }
    r
  },error=function(e)e)
  if(inherits(r,"error")) {
    attempts$status[i] <- "invalid"
    attempts$error[i] <- conditionMessage(r)
  } else {
    attempts$status[i] <- r$status
    attempts$elapsed_seconds[i] <- r$elapsed_seconds
    records[[i]] <- r
  }
}
digests <- unique(vapply(Filter(Negate(is.null),records),function(r)r$runner_md5,character(1)))
if(length(digests)>1L) {
  attempts$status[attempts$status=="ok"] <- "invalid"
  attempts$error[attempts$error==""] <- "Mixed runner digests"
}
write.csv(attempts,file.path(out,"attempts.csv"),row.names=FALSE)
if(partial) {
  if(any(attempts$status %in% c("invalid","failed"))) stop("Prior attempts need investigation")
  quit(status=0)
}
if(any(attempts$status!="ok")) stop("Incomplete or failed attempts; see attempts.csv")
by_seed <- lapply(c(0L,prepared$fit_seeds),function(s) {
  x <- do.call(rbind,lapply(records[jobs$seed==s],`[[`,"predictions"))
  x[order(x$row_n),]
})
names(by_seed) <- as.character(c(0L,prepared$fit_seeds))
truth <- by_seed[[1]]
for(x in by_seed) stopifnot(identical(x$row_n,truth$row_n),identical(x$observed,truth$observed))
P <- sapply(by_seed[-1],`[[`,"predicted")
spread <- apply(P,1,function(x)diff(range(x)))
write.csv(data.frame(n=nrow(P),species=length(unique(truth$species)),
  median_seed_range=median(spread),p90_seed_range=unname(quantile(spread,.9)),
  max_seed_range=max(spread),n_range_over_3=sum(spread>3)),
  file.path(out,"seed_spread.csv"),row.names=FALSE)
methods <- lapply(by_seed,`[[`,"predicted")
names(methods) <- paste0("seed_",names(methods))
methods$trio_51001_51003 <- rowMeans(P[,1:3,drop=FALSE])
methods$trio_51004_51006 <- rowMeans(P[,4:6,drop=FALSE])
metric <- list()
for(name in names(methods)) for(without_newt in c(FALSE,TRUE)) {
  keep <- if(without_newt) truth$species!="Notophthalmus viridescens" else rep(TRUE,nrow(truth))
  se <- (methods[[name]][keep]-truth$observed[keep])^2
  metric[[length(metric)+1L]] <- data.frame(method=name,without_newt=without_newt,
    n=sum(keep),n_species=length(unique(truth$species[keep])),
    cell_rmse=sqrt(mean(se)),species_balanced_rmse=sqrt(mean(tapply(se,truth$species[keep],mean))))
}
write.csv(do.call(rbind,metric),file.path(out,"point_metrics.csv"),row.names=FALSE)
interval <- lapply(names(by_seed),function(s) {
  x <- by_seed[[s]]
  ok <- is.finite(x$lower)&is.finite(x$upper)&x$lower<=x$upper
  data.frame(seed=s,n_expected=nrow(x),n_interval=sum(ok),
    covered=sum(x$observed[ok]>=x$lower[ok]&x$observed[ok]<=x$upper[ok]),
    mean_width=if(any(ok))mean(x$upper[ok]-x$lower[ok]) else NA_real_)
})
write.csv(do.call(rbind,interval),file.path(out,"individual_intervals.csv"),row.names=FALSE)
saveRDS(list(source_sha=records[[1]]$source_sha,split_seed=prepared$split_seed,
  input_md5=prepared$input_md5,runner_md5=lapply(records,`[[`,"runner_md5"),
  warnings=lapply(records,`[[`,"warnings"),session=sessionInfo()),file.path(out,"provenance.rds"))
cat("ACCOUNTED",nrow(jobs),"fits; all output finite\n")
