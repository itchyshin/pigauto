# Synthetic harness checks only. No imputation model is fitted.
base <- tempfile("seed260-resume-")
dir.create(base)
runner <- file.path(base,"pigauto_seed260","runner")
dir.create(runner,recursive=TRUE)
src <- file.path(base,"pigauto_seed260","source")
dir.create(src)
sha <- "033649b93115c6d8191e791f06b1ddf52eb82d19"
writeLines(sha,file.path(src,"SOURCE_SHA"))
file.copy(c("script/seed_validation_260/run.R","script/seed_validation_260/summarize.R"),runner)
a <- file.path(base,"sim-archive/pigauto/seed260-20261108")
dir.create(file.path(a,"input"),recursive=TRUE)
res <- file.path(a,"production")
dir.create(res)
p <- list(data=data.frame(row_n=1:10,species=rep(letters[1:5],each=2),CTmax=seq_len(10)+.1),
          fold=rep(1:5,each=2),fit_seeds=51001:51006,split_seed=2601109L,input_md5="synthetic-only")
pp <- file.path(a,"input/prepared.rds")
saveRDS(p,pp)
hash <- unname(tools::md5sum(pp))
rh <- unname(tools::md5sum(file.path(runner,"run.R")))
for(s in c(0L,p$fit_seeds)) for(f in 1:5) {
  d <- p$data[p$fold==f,]
  r <- list(fold=f,seed=s,n_expected=nrow(d),epochs_requested=2000L,split_seed=p$split_seed,
    source_sha=sha,prepared_md5=hash,runner_md5=rh,status="ok",elapsed_seconds=1,
    predictions=data.frame(row_n=d$row_n,species=d$species,observed=d$CTmax,
                           predicted=d$CTmax+.1,lower=d$CTmax-1,upper=d$CTmax+1))
  saveRDS(r,file.path(res,sprintf("fold%d-seed%d-epochs2000.rds",f,s)))
}
run_summary <- function(label) system2("Rscript",c(shQuote(file.path(runner,"summarize.R")),
  shQuote(pp),shQuote(res),shQuote(file.path(a,label))),stdout=TRUE,stderr=TRUE)
stopifnot(is.null(attr(run_summary("good"),"status")))
m <- read.csv(file.path(a,"good/point_metrics.csv"))
stopifnot(nrow(m)==18,all(abs(m$cell_rmse-.1)<1e-8))
path <- file.path(res,"fold1-seed0-epochs2000.rds")
orig <- readRDS(path)
for(kind in c("prepared","runner","corrupt","missing")) {
  r <- orig
  if(kind=="prepared") {r$prepared_md5<-"wrong";saveRDS(r,path)}
  if(kind=="runner") {r$runner_md5<-paste(rep("b",32),collapse="");saveRDS(r,path)}
  if(kind=="corrupt") writeLines("broken",path)
  if(kind=="missing") unlink(path)
  x <- suppressWarnings(run_summary(kind))
  stopifnot(attr(x,"status")!=0)
  ledger <- read.csv(file.path(a,kind,"attempts.csv"))
  stopifnot(nrow(ledger)==35,any(ledger$status!="ok"))
  saveRDS(orig,path)
}
launcher <- normalizePath("script/seed_validation_260/launch.py")
python <- paste0("from pathlib import Path\nfrom unittest.mock import patch\nimport runpy,sys\n",
  "sys.argv=['launch.py','--plan-only']\nwith patch.object(Path,'home',return_value=Path(",
  encodeString(base,quote='"'),")):\n runpy.run_path(",encodeString(launcher,quote='"'),",run_name='__main__')")
plan <- function() system2("python3",c("-c",shQuote(python)),stdout=TRUE,stderr=TRUE)
x <- plan()
stopifnot(is.null(attr(x,"status")),any(grepl('"reused": 35',x,fixed=TRUE)))
unlink(path)
x <- plan()
stopifnot(is.null(attr(x,"status")),any(grepl('"queued": 1',x,fixed=TRUE)))
writeLines("wrong",file.path(src,"SOURCE_SHA"))
x <- suppressWarnings(plan())
stopifnot(attr(x,"status")!=0)
saveRDS(orig,path)
run_strata <- function(label) system2("Rscript",c(
  shQuote("script/seed_validation_260/stratify.R"),shQuote(pp),shQuote(res),
  shQuote(file.path(a,paste0(label,".csv")))),stdout=TRUE,stderr=TRUE)
stopifnot(is.null(attr(run_strata("strata-good"),"status")))
for(kind in c("source","runner","truth","interval")) {
  r <- orig
  if(kind=="source") r$source_sha <- "wrong"
  if(kind=="runner") r$runner_md5 <- "wrong"
  if(kind=="truth") r$predictions$observed[1] <- -99
  if(kind=="interval") r$predictions$lower[1] <- Inf
  saveRDS(r,path)
  x <- suppressWarnings(run_strata(paste0("strata-",kind)))
  stopifnot(!is.null(attr(x,"status")),attr(x,"status")!=0)
}
saveRDS(orig,path)
cat("HARNESS SELFTEST PASS: complete summary, four invalid/missing controls, resume, source gate\n")
