validation_job_hash <- function(job) digest::digest(as.list(job),algo="sha256")
validation_check_resume <- function(record,job,source,runtime) {
 if(!identical(record$source_identity,source)||!identical(record$runtime_identity,runtime)||!identical(record$job_spec_hash,validation_job_hash(job)))stop("Stale or mismatched job result")
 invisible(TRUE)
}
validation_check_approval <- function(receipt,source,runtime,jobs,workers,budget,timeout,context) {
 want <- list(approved=TRUE,source_identity=source,runtime_identity=runtime,jobs_hash=digest::digest(jobs,algo="sha256"),workers=as.integer(workers),budget_sec=as.numeric(budget),job_timeout_sec=as.numeric(timeout),execution_context=context)
 for(nm in names(want))if(!identical(receipt[[nm]],want[[nm]]))stop("Approval does not bind ",nm)
 invisible(TRUE)
}
validation_check_resources <- function(context,host,workers,slurm_id,cpus) {
 if(context=="totoro") {
   if(!grepl("totoro",host,ignore.case=TRUE))stop("Totoro context requires the Totoro host")
 } else if(context=="drac") {
   if(!nzchar(slurm_id)||!is.finite(suppressWarnings(as.numeric(cpus)))||as.numeric(cpus)<workers)stop("DRAC needs a sufficient scheduler allocation")
 } else stop("Explicit totoro or drac execution context required")
 invisible(TRUE)
}
