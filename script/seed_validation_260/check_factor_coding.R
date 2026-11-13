args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
devtools::load_all(args[1],quiet=TRUE)
p <- readRDS(file.path(args[2],"prepared.rds"))
path <- file.path(args[2],"original_cv_folds.rds")
if(!file.exists(path)) download.file(paste0("https://raw.githubusercontent.com/daniel1noble/BACE/",
  "b4b3cde02c4b171947042100e26b13b9f7ac425d/vignettes/amphibian_ctmax/data/cv_folds.rds"),path,mode="wb")
folds <- readRDS(path)
d <- p$data
numeric_covs <- as.data.frame(model.matrix(~ .,p$covariates)[,-1,drop=FALSE])
warnings <- list()
rows <- lapply(1:5,function(k) {
  traits <- data.frame(species=d$species,CTmax=d$CTmax)
  traits$CTmax[d$row_n %in% folds[[k]]$rows_masked] <- NA_real_
  fit <- function(covs,label) withCallingHandlers(pigauto::impute(traits,p$tree,
    species_col="species",covariates=covs,gnn=FALSE,log_transform=FALSE,
    seed=2026L+k,verbose=FALSE),warning=function(w) {
      warnings[[paste(k,label)]] <<- c(warnings[[paste(k,label)]],conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  tryCatch({
    a <- fit(p$covariates,"factor")
    b <- fit(numeric_covs,"numeric")
    data.frame(fold=k,status="ok",error="",lambda_factor=unname(a$baseline$lambda_per_trait[1]),
      lambda_numeric=unname(b$baseline$lambda_per_trait[1]),
      max_completed_difference=max(abs(a$completed$CTmax-b$completed$CTmax)),
      max_se_difference=max(abs(a$prediction$se-b$prediction$se)))
  },error=function(e)data.frame(fold=k,status="failed",error=conditionMessage(e),
     lambda_factor=NA_real_,lambda_numeric=NA_real_,max_completed_difference=NA_real_,
     max_se_difference=NA_real_))
})
out <- do.call(rbind,rows)
dir.create(args[3],recursive=TRUE,showWarnings=FALSE)
write.csv(out,file.path(args[3],"factor_coding.csv"),row.names=FALSE)
saveRDS(list(warnings=warnings,input_md5=tools::md5sum(path),
  source_sha=readLines(file.path(args[1],"SOURCE_SHA")),session=sessionInfo()),
  file.path(args[3],"factor_coding_provenance.rds"))
print(out)
stopifnot(all(out$status=="ok"),all(is.finite(out$lambda_factor)),
          all(abs(out$lambda_factor-out$lambda_numeric)<1e-5),
          all(out$max_completed_difference<1e-5),all(out$max_se_difference<1e-5))
cat("FACTOR CODING APPLICATION CHECK PASS\n")
