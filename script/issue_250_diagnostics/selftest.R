source("script/issue_250_diagnostics/diagnose_ties.R")
tree <- ape::read.tree(text="((a:1,b:1):1,(c:1,d:1):1);")
d <- data.frame(y=c(1,2,NA,NA),row.names=letters[1:4])
r <- list(tree=tree,completed=data.frame(y=c(1,2,1.5,1.5),row.names=letters[1:4]),
  data=list(multi_obs=FALSE),fit=list(baseline_full=list(
    predict_method_by_trait=c(y="per_column"),lambda_per_trait=c(y=1))))
b <- list(traits=d,tree=tree,result=r,trait="y")
x <- diagnose_ties(b)
stopifnot(x$summary$n_exact_tied==2,x$pairs$max_observed_correlation_difference==0)
b$result$completed$y[4] <- 1.5+1e-10
x <- diagnose_ties(b)
stopifnot(x$summary$n_exact_tied==0,x$summary$n_rounded_tied==2,nrow(x$pairs)==0)
b$result$completed$y[4] <- NA
stopifnot(inherits(try(diagnose_ties(b),silent=TRUE),"try-error"))
b$result <- r
b$tree$edge.length[1] <- 2
stopifnot(inherits(try(diagnose_ties(b),silent=TRUE),"try-error"))
b$tree <- tree
b$traits$y[1] <- 99
stopifnot(inherits(try(diagnose_ties(b),silent=TRUE),"try-error"))
cat("TIE DIAGNOSTIC SELFTEST PASS\n")
