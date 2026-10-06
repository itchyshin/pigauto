args <- commandArgs(TRUE); indir <- args[[1L]]; out <- args[[2L]]
rows <- list()
backends <- if (length(args)>=3L) args[[3L]] else c("drm","gllvm")
stopifnot(all(backends %in% c("drm","gllvm")))
for (backend in backends) {
  loadNamespace(if(backend=="drm") "drmTMB" else "gllvmTMB")
  seeds <- 2026100601:2026100610
  files <- file.path(indir,paste0(backend,"-",seeds,".rds"))
  if (!all(file.exists(files))) stop("Missing registered results: ",paste(files[!file.exists(files)],collapse=", "))
  records <- lapply(files,readRDS)
  stopifnot(identical(vapply(records,function(x)x$seed,integer(1)),seeds),
            all(vapply(records,function(x)length(x$fits)==20L && x$m==20L && x$n==300L,logical(1))))
  terms <- names(records[[1L]]$truth)
  for (term in terms) {
    q <- vapply(records,function(x)x$oracle$estimate[x$oracle$term==term],numeric(1))
    se <- vapply(records,function(x)sqrt(x$oracle$totalvar[x$oracle$term==term]),numeric(1))
    truth <- records[[1L]]$truth[[term]];bias <- q-truth
    covered <- vapply(records,function(x) {
      r <- x$oracle[x$oracle$term==term,];r$conf.low<=truth && r$conf.high>=truth
    },logical(1))
    full <- vapply(records,function(x) {
      if (backend=="drm") {
        bb <- coef(x$fullfit); vv <- unlist(bb,use.names=FALSE)
        names(vv) <- unlist(lapply(names(bb),function(k)paste0(k,":",names(bb[[k]]))))
        vv[[term]]
      } else {
        td <- getS3method("tidy","gllvmTMB_multi",envir=asNamespace("gllvmTMB"))(x$fullfit,effects="fixed")
        td$estimate[td$term==term]
      }
    },numeric(1))
    mcse <- sd(bias)/sqrt(10);lo <- mean(bias)-qt(.975,9)*mcse;hi <- mean(bias)+qt(.975,9)*mcse
    ci <- binom.test(sum(covered),10)$conf.int
    rows[[length(rows)+1L]] <- data.frame(backend=backend,term=term,seeds=10L,m=20L,n=300L,truth=truth,
      mean_estimate=mean(q),bias=mean(bias),bias_mcse=mcse,bias_low=lo,bias_high=hi,
      empirical_sd=sd(q),mean_pooled_se=mean(se),coverage_count=sum(covered),coverage_n=10L,
      coverage_low=ci[1],coverage_high=ci[2],full_data_bias=mean(full-truth),
      bias_margin_pass=lo>-.15 && hi<.15)
  }
}
res <- do.call(rbind,rows);write.csv(res,out,row.names=FALSE);print(res,row.names=FALSE)
if (!all(res$bias_margin_pass)) stop("BOUND_RECOVERY_FAILED_OR_INCONCLUSIVE: bias interval outside +/- .15; no silent expansion")
cat("BOUNDED_RECOVERY_OK; ten seeds are not a coverage certificate\n")
