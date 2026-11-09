# Post-hoc diagnostic after the predeclared summaries; no fitting or tuning.
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
p <- readRDS(args[1])
seen <- logical(nrow(p$data))
for(f in 1:5) seen[p$fold==f] <- p$data$species[p$fold==f] %in% p$data$species[p$fold!=f]
by_seed <- lapply(c(0L,p$fit_seeds),function(s) {
  x <- do.call(rbind,lapply(1:5,function(f) {
    r <- readRDS(file.path(args[2],sprintf("fold%d-seed%d-epochs2000.rds",f,s)))
    stopifnot(r$status=="ok",identical(r$prepared_md5,unname(tools::md5sum(args[1]))))
    r$predictions
  }))
  x[order(x$row_n),]
})
truth <- by_seed[[1]]
stopifnot(!anyDuplicated(truth$row_n),nrow(truth)==nrow(p$data))
for(x in by_seed) stopifnot(identical(x$row_n,truth$row_n),identical(x$observed,truth$observed))
seen <- seen[match(truth$row_n,p$data$row_n)]
P <- sapply(by_seed[-1],`[[`,"predicted")
methods <- lapply(by_seed,`[[`,"predicted")
names(methods) <- paste0("seed_",c(0L,p$fit_seeds))
methods$trio_51001_51003 <- rowMeans(P[,1:3,drop=FALSE])
methods$trio_51004_51006 <- rowMeans(P[,4:6,drop=FALSE])
out <- list()
for(method in names(methods)) for(group in c(FALSE,TRUE)) {
  keep <- seen==group
  se <- (methods[[method]][keep]-truth$observed[keep])^2
  out[[length(out)+1L]] <- data.frame(method=method,species_has_training_data=group,
    n=sum(keep),n_species=length(unique(truth$species[keep])),
    cell_rmse=sqrt(mean(se)),species_balanced_rmse=sqrt(mean(tapply(se,truth$species[keep],mean))))
}
write.csv(do.call(rbind,out),args[3],row.names=FALSE)
cat("POSTHOC_STRATA_ACCOUNTED",nrow(truth),"rows\n")
