sp <- commandArgs(TRUE)[1]
source("script/mi_gls/regimes.R"); rg <- regimes[, c("regime_id","lambda","n","mechanism","missing")]
a <- read.csv(file.path(sp,"sim_summary_iw.csv")); b <- read.csv(file.path(sp,"sim_summary.csv"))
k <- c("regime_id","method","downstream"); v <- c("paired_bias","paired_bias_mcse","coverage","se_ratio","complete_coverage","n_converged")
m <- merge(a[a$method=="posterior_full",c(k,v)], b[b$method=="posterior_full",c(k,v)], by=k, suffixes=c(".iw",".sep")); m <- merge(m, rg, by="regime_id")
m$g <- m$regime_id >= 17
for (g in c(TRUE,FALSE)) { x <- m[m$g==g,]; cat(if (g) "\nGATED 17-40" else "\nSTRESS 1-16", nrow(x), "rows\n")
 cat(" |bias| mean iw", mean(abs(x$paired_bias.iw)), "sep", mean(abs(x$paired_bias.sep)), "; max iw", max(abs(x$paired_bias.iw)), "sep", max(abs(x$paired_bias.sep)), "\n")
 cat(" coverage range iw", range(x$coverage.iw), "sep", range(x$coverage.sep), "\n")
 cat(" coverage - complete min iw", min(x$coverage.iw-x$complete_coverage.iw), "sep", min(x$coverage.sep-x$complete_coverage.sep), "\n")
 cat(" converged iw", sum(x$n_converged.iw), "sep", sum(x$n_converged.sep), "of", 200*nrow(x), "\n") }
x <- m[m$g & m$lambda==1,]; x$z.iw <- x$paired_bias.iw/x$paired_bias_mcse.iw; x$z.sep <- x$paired_bias.sep/x$paired_bias_mcse.sep
cat("\nlambda = 1 gated rows\n"); print(x[order(x$regime_id),c("regime_id","downstream","n","missing","paired_bias.iw","z.iw","paired_bias.sep","z.sep","coverage.iw","coverage.sep")], digits=3, row.names=FALSE)
x <- m[m$g & m$lambda!=1,]; cat("\nlambda != 1 gated: bias iw", range(x$paired_bias.iw), " sep", range(x$paired_bias.sep), "; cov iw", range(x$coverage.iw), "sep", range(x$coverage.sep), "\n")
c1 <- read.csv(file.path(sp,"cell_coverage_iw.csv")); c2 <- read.csv(file.path(sp,"cell_coverage.csv"))
f <- function(d) { d <- d[d$method=="posterior_full" & d$mechanism=="MCAR" & d$regime_id>=17,]; d$cov <- d$covered_sum/d$n; d }
cat("\nG7 gated MCAR rows: iw", range(f(c1)$cov), " sep", range(f(c2)$cov), " n rows", nrow(f(c1)), nrow(f(c2)), "\n")
