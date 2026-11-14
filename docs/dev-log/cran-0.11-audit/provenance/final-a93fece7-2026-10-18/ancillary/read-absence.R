args <- commandArgs(TRUE)
root <- normalizePath(args[1],mustWork=TRUE)
library(testthat)
objects <- readRDS(file.path(root,"absence-results.rds"))
d <- do.call(rbind,lapply(objects$results,as.data.frame))
flat <- d[,c("file","test","failed","skipped","error","warning","passed")]
stopifnot(!any(vapply(flat,is.list,logical(1))))
write.csv(flat,file.path(root,"absence-results-flat.csv"),row.names=FALSE)
summary <- jsonlite::read_json(file.path(root,"absence-summary.json"),simplifyVector=TRUE)
stopifnot(sum(flat$passed)==summary$pass,sum(flat$failed)==summary$fail,
          sum(flat$warning)==summary$warn,sum(flat$skipped)==summary$skip,
          length(objects$problems)==summary$fail,length(objects$warnings)==summary$warn,
          length(objects$skips)==summary$skip,!any(flat$error))
jsonlite::write_json(list(pass=sum(flat$passed),fail=sum(flat$failed),warnings=sum(flat$warning),
                         skips=sum(flat$skipped),tests=nrow(flat),
                         scope="Read retained actual result objects only; no tests rerun"),
                    file.path(root,"absence-rds-readback.json"),pretty=TRUE,auto_unbox=TRUE)
cat("ACTUAL_ABSENCE_RDS_AND_FLAT_COUNTS_AGREE\n")
