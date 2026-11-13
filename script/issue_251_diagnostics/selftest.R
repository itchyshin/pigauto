source("script/issue_251_diagnostics/subset_intervals.R")
x <- list(tree1=list(completed_draws=setNames(lapply(1:6,function(i)
  seq(-.2,.2,length.out=100)+i/100),paste0("mi",1:6)),
  reference_draws=seq(-.3,.3,length.out=600)))
a <- subset_intervals(x,c(2,6),3)
stopifnot(nrow(a)==4,all(is.finite(a$width_ratio)),sum(a$m==6)==1)
bad <- x
bad[[1]]$completed_draws[[1]] <- 2:3
stopifnot(inherits(try(subset_intervals(bad),silent=TRUE),"try-error"))
bad <- x
bad[[1]]$reference_draws <- rep(0,20)
stopifnot(inherits(try(subset_intervals(bad),silent=TRUE),"try-error"))
for(grid in list(c(2,20),Inf,NA_real_))
  stopifnot(inherits(try(subset_intervals(x,grid),silent=TRUE),"try-error"))
cat("INTERVAL SUBSET SELFTEST PASS\n")
