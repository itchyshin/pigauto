#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
arg <- function(name, default=NULL) { i <- match(name,args); if(is.na(i)) default else { if(i==length(args)) stop("Missing value: ",name); args[i+1L] } }
root <- normalizePath(arg("--root", "."), mustWork=TRUE)
setwd(root)
base <- file.path("script","validation_20261018")
source(file.path(base,"orchestration.R")); source(file.path(base,"manifest.R")); source(file.path(base,"freeze.R"))
mode <- arg("--mode","manifest")
if (!mode %in% c("manifest","pilot","production","job","approval-check")) stop("Invalid mode")
manifest <- validation_manifest(root); validate_manifest(manifest)
jobs <- validation_jobs(manifest)
out <- arg("--out",file.path(base,"results",mode))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
identity <- validation_source_identity(root)
identity$runtime <- validation_runtime_identity()
if(mode=="manifest") {
 write.csv(manifest,file.path(out,"scenarios.csv"),row.names=FALSE)
 write.csv(jobs,file.path(out,"jobs.csv"),row.names=FALSE)
 saveRDS(jobs,file.path(out,"jobs.rds"))
 saveRDS(identity,file.path(out,"identity.rds")); cat("MANIFEST_FROZEN scenarios=",nrow(manifest)," jobs=",nrow(jobs)," hash=",identity$hash,"\n",sep=""); quit(status=0)
}
if(mode=="job") {
 job <- readRDS(arg("--job-file")); source(file.path(base,"runners.R"))
 started <- proc.time()[["elapsed"]]
 if(requireNamespace("torch",quietly=TRUE))try({torch::torch_set_num_threads(1L);torch::torch_set_num_interop_threads(1L)},silent=TRUE)
 result <- tryCatch(validation_run_job(job,root=root,pilot=FALSE), error=function(e) .vr_envelope(job,"failed",.vr_row(job,job$adapter,"attempt","failed",failure_reason=conditionMessage(e)),error=conditionMessage(e)))
 for(nm in c("job_id","scenario_id","rep","seed","slice","adapter")) result[[nm]] <- job[[nm]][1]
 result$job_spec_hash <- validation_job_hash(job)
 result$source_identity <- identity$hash; result$source_manifest <- identity
 result$runtime_identity <- identity$runtime$hash
 result$settings$mode <- arg("--record-mode","pilot")
 result$wall_sec <- proc.time()[["elapsed"]]-started
 result$session_info <- utils::sessionInfo(); result$host <- Sys.info()[["nodename"]]
 dest <- arg("--result-file"); tmp <- paste0(dest,".tmp-",Sys.getpid()); saveRDS(result,tmp); if(!file.rename(tmp,dest)) stop("Atomic result rename failed")
 cat("JOB_RECORDED ",job$job_id," status=",result$status,"\n",sep=""); quit(status=0)
}
if(mode %in% c("production","approval-check")) {
 approval <- arg("--approval")
 if(is.null(approval)||!file.exists(approval)) stop("Production requires the measured-pilot approval receipt")
 receipt <- readRDS(approval)
 workers <- as.integer(arg("--workers","1")); worker_index <- as.integer(arg("--worker-index","1"))
 budget <- as.numeric(arg("--budget-sec","1800")); job_timeout <- as.numeric(arg("--job-timeout","600"))
 context <- Sys.getenv("VALIDATION_CONTEXT","")
 if(anyNA(c(workers,worker_index,budget,job_timeout))||workers<1L||workers>48L||worker_index<1L||worker_index>workers||budget<=0||job_timeout<=0)stop("Invalid bounded worker assignment")
 validation_check_resources(context,Sys.info()[["nodename"]],workers,Sys.getenv("SLURM_JOB_ID"),Sys.getenv("SLURM_CPUS_PER_TASK"))
 validation_check_approval(receipt,identity$hash,identity$runtime$hash,jobs,workers,budget,job_timeout,context)
 if(mode=="approval-check"){cat("APPROVAL_BOUND\n");quit(status=0)}
 jobs <- jobs[((seq_len(nrow(jobs))-1L) %% workers)+1L==worker_index,,drop=FALSE]
 out <- file.path(out,sprintf("worker-%03d",worker_index));dir.create(out,recursive=TRUE,showWarnings=FALSE)
}
if(mode=="pilot") {
 required <- c("prediction","zi_count","multi_proportion","multiobs","multiobs_mixed","covariate","tree","posterior","analysis","real")
 if(!all(required %in% jobs$adapter))stop("Pilot manifest lacks a required adapter")
 selected <- !duplicated(paste(jobs$adapter,ifelse(is.na(jobs$model),"",jobs$model)))
 selected <- selected | (jobs$adapter=="real" & jobs$rep==1L)
 for(adapter in c("prediction","posterior")) {
   idx <- which(jobs$adapter==adapter & jobs$n==1000)
   if(!length(idx))stop("Pilot lacks large-n ",adapter)
   selected[idx[1]] <- TRUE
 }
 idx <- which(jobs$slice=="prediction_low_signal")
 if(!length(idx))stop("Pilot lacks low-signal sentinel")
 selected[idx[1]] <- TRUE
 if ("--supplement" %in% args) {
   selected <- jobs$adapter == "analysis" & jobs$rep == 1L
   selected <- selected | (jobs$adapter == "real" & jobs$scenario %in% c("avonet","bien") & jobs$rep == 1L)
   for (adapter in c("tree","posterior")) selected[which(jobs$adapter == adapter)[1L]] <- TRUE
   selected[which(jobs$adapter == "prediction" & jobs$mechanism == "mar")[1L]] <- TRUE
 }
 jobs <- jobs[selected,,drop=FALSE]
 if(!all(c("lm","glm","lmer") %in% jobs$model[jobs$adapter=="analysis"]))stop("Pilot lacks an analysis engine")
}

write.csv(jobs,file.path(out,"jobs.csv"),row.names=FALSE)
saveRDS(jobs,file.path(out,"jobs.rds")); saveRDS(identity,file.path(out,"identity.rds"))
limit <- as.numeric(arg("--job-timeout","600")); budget <- as.numeric(arg("--budget-sec","1800"))
if(!is.finite(limit)||limit<=0||!is.finite(budget)||budget<=0)stop("Invalid time budget")
started <- proc.time()[["elapsed"]]
for(i in seq_len(nrow(jobs))) {
 job <- jobs[i,,drop=FALSE]; dest <- file.path(out,paste0(job$job_id,".rds"))
 if(file.exists(dest)) {
  old <- readRDS(dest)
  validation_check_resume(old,job,identity$hash,identity$runtime$hash)
  next
 }
 if(proc.time()[["elapsed"]]-started>=budget) { cat("PILOT_BUDGET_REACHED\n"); break }
 jf <- tempfile("validation-job-",tmpdir=out,fileext=".rds"); saveRDS(job,jf)
 log <- file.path(out,paste0(job$job_id,".log"))
 remaining <- budget-(proc.time()[["elapsed"]]-started)
 cmdargs <- c(file.path(base,"run.R"),"--root",root,"--mode","job","--job-file",jf,"--result-file",dest,"--record-mode",mode)
 child_start <- proc.time()[["elapsed"]]
 timeout_bin <- Sys.which("timeout")
 if(!nzchar(timeout_bin))stop("GNU timeout is required for bounded jobs")
 timeout_args <- c("--foreground","--signal=KILL",paste0(max(1,min(limit,remaining)),"s"),file.path(R.home("bin"),"Rscript"),cmdargs)
 status <- system2(timeout_bin,shQuote(timeout_args),stdout=log,stderr=log,timeout=0)
 unlink(jf)
 if(!file.exists(dest)) {
  source(file.path(base,"runners.R"))
  reason <- if(status %in% c(124,137))"Worker timed out or was killed" else paste("Worker exited",status)
  result <- .vr_envelope(job,"failed",.vr_row(job,job$adapter,"attempt","failed",failure_reason=reason),error=reason)
  result$job_spec_hash <- validation_job_hash(job);result$source_identity <- identity$hash;result$runtime_identity <- identity$runtime$hash;result$wall_sec <- proc.time()[["elapsed"]]-child_start;saveRDS(result,dest)
 }
 record <- readRDS(dest);record$process_wall_sec <- proc.time()[["elapsed"]]-child_start
 tmp <- paste0(dest,".parent-tmp");saveRDS(record,tmp);if(!file.rename(tmp,dest))stop("Parent result rename failed")
 cat("PILOT_PROGRESS ",i,"/",nrow(jobs)," ",job$job_id,"\n",sep=""); flush.console()
}
cat("RUN_ACCOUNTED mode=",mode," elapsed_sec=",proc.time()[["elapsed"]]-started,"\n",sep="")
