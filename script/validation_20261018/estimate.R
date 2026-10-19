#!/usr/bin/env Rscript
source("script/validation_20261018/manifest.R")
args <- commandArgs(TRUE)
inputs <- if(length(args)) args else "script/validation_20261018/results/pilot"
out <- inputs[1L]
m <- validation_manifest(".")
pilot <- do.call(rbind,lapply(inputs,function(input) {
 jobs <- readRDS(file.path(input,"jobs.rds"))
 records <- lapply(file.path(input,paste0(jobs$job_id,".rds")),readRDS)
 identity <- readRDS(file.path(input,"identity.rds"))
 if(any(!vapply(records,function(x)identical(x$source_identity,identity$hash),logical(1))))stop("Mixed raw source identity in timing input")
 jobs$wall_sec <- vapply(records,function(x)if(!is.null(x$process_wall_sec))x$process_wall_sec else x$wall_sec,numeric(1))
 jobs$status <- vapply(records,`[[`,character(1),"status")
 jobs$timing_input <- input
 jobs$raw_source_identity <- identity$hash
 jobs
}))
write.csv(pilot,file.path(out,"pilot_timings.csv"),row.names=FALSE)
rows <- lapply(seq_len(nrow(m)),function(i){
 s <- m[i,,drop=FALSE];ix <- which(pilot$adapter==s$adapter)
 if(!is.na(s$model))ix <- ix[!is.na(pilot$model[ix]) & pilot$model[ix]==s$model]
 if(s$adapter=="real")ix <- ix[pilot$scenario[ix]==s$scenario]
 if(!length(ix))stop("Missing pilot route for ",s$scenario_id)
 if(!is.na(s$n)) {
   exact <- ix[!is.na(pilot$n[ix]) & pilot$n[ix]==s$n]
   if(length(exact))ix <- exact
 }
 mechanism_pooled <- FALSE
 if(!is.na(s$mechanism)) {
   matched <- ix[!is.na(pilot$mechanism[ix]) & tolower(pilot$mechanism[ix]) == tolower(s$mechanism)]
   if(length(matched))ix <- matched else mechanism_pooled <- TRUE
 }
 pilot_had_failures <- any(pilot$status[ix] != "ok")
 successful <- ix[pilot$status[ix] == "ok"]
 pilot_failed <- !length(successful)
 if(length(successful))ix <- successful
 t <- max(pilot$wall_sec[ix]);scale <- 1
 if(!is.na(s$n) && all(!is.na(pilot$n[ix])))scale <- max(1,(s$n/max(pilot$n[ix]))^3)
 data.frame(scenario_id=s$scenario_id,adapter=s$adapter,reps=s$reps,
    seconds_per_job=t*scale,estimated_serial_sec=s$reps*t*scale,
    pilot_failed=pilot_failed,pilot_had_failures=pilot_had_failures,
    mechanism_pooled=mechanism_pooled,extrapolated=scale>1)
})
tab <- do.call(rbind,rows);write.csv(tab,file.path(out,"production_estimate.csv"),row.names=FALSE)
serial <- sum(tab$estimated_serial_sec);workers <- 48L
lower <- serial/workers/3600;upper <- 2*serial/(workers*.75)/3600
text <- c("# Measured pilot: provisional production estimate","",
 sprintf("Planned jobs: %d. Nominal workers: %d. Measured serial-equivalent work: %.1f CPU-hours.",sum(m$reps),workers,serial/3600),
 sprintf("Provisional elapsed range: %.1f to %.1f hours at 48 single-thread workers, before queue time.",lower,upper),
 "Timing-only inputs retain separate raw source identities. Successful measured routes supersede failed adapter startup timings; all timing attempts remain in pilot_timings.csv.",
 "The upper planning estimate applies a 2x timing margin and 75% utilization. This is a pilot extrapolation, not a runtime guarantee. Scenarios without exact size matches use conservative cubic scaling. Failed/unavailable pilot routes remain flagged; their timings do not establish successful production feasibility.",
 sprintf("Scenarios without a successful measured pilot: %d. Scenarios retaining a failed/unavailable timing attempt: %d. Size-extrapolated scenarios: %d. Mechanism-pooled timings: %d.",sum(tab$pilot_failed),sum(tab$pilot_had_failures),sum(tab$extrapolated),sum(tab$mechanism_pooled)),
 "No production approval receipt has been created. Recheck host load, aggregate user CPU use and measured memory before scheduling.")
writeLines(text,file.path(out,"estimate.md"));cat(paste(text,collapse="\n"),"\n")
