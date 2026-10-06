# Build-excluded, pre-registered bounded Gaussian recovery.
args <- commandArgs(TRUE)
backend <- args[[1L]]
seed <- as.integer(args[[2L]])
outdir <- args[[3L]]
stopifnot(backend %in% c("drm", "gllvm"), seed %in% 2026100601:2026100610)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
library(pigauto)
close_enough <- function(actual, expected, tol, label) {
  err <- max(abs(actual - expected) / pmax(1, abs(expected)))
  if (!is.finite(err) || err > tol) stop(label, ": relative error ", err)
  err
}
# Explicit, independent deterministic streams; data seed is the registered seed.
streams <- c(data = seed, tree = seed - 100000L,
             mask = seed - 200000L, imputation = seed - 300000L)
n <- 300L; m <- 20L
set.seed(streams[["tree"]]); tree <- ape::rtree(n)
set.seed(streams[["data"]]); x <- rnorm(n)
if (backend == "drm") {
  full <- data.frame(x = x, y = 1 + .7*x + .8*rnorm(n), row.names = tree$tip.label)
  truth <- c("mu:(Intercept)" = 1, "mu:x" = .7, "sigma:(Intercept)" = log(.8))
} else {
  alpha <- c(.5,-.5,1,0); beta <- c(.7,-.4,.3,.5); loading <- c(.6,.4,-.5,.7)
  z <- rnorm(n); eps <- matrix(rnorm(n*4L), n)
  y <- matrix(alpha, n, 4L, byrow = TRUE) + outer(x,beta) + outer(z,loading) + .8*eps
  full <- data.frame(x=x,y, row.names=tree$tip.label)
  names(full) <- c("x",paste0("y",1:4))
  truth <- setNames(c(alpha,beta), c(paste0("traity",1:4),paste0("traity",1:4,":x")))
}
masked <- full
set.seed(streams[["mask"]])
for (nm in setdiff(names(full), "x")) masked[sample.int(n,60L),nm] <- NA_real_
stopifnot(all(colSums(is.na(masked))[-1L] == 60L), !anyNA(masked$x))
fit_one <- if (backend == "drm") function(d) {
  drmTMB::drmTMB(drmTMB::bf(y ~ x, sigma ~ 1), data=d, family=gaussian(), REML=FALSE, engine="tmb")
} else function(d) {
  tr <- paste0("y",1:4)
  long <- data.frame(unit=factor(rep(rownames(d),times=4L)),
                     trait=factor(rep(tr,each=n),levels=tr),
                     x=rep(d$x,times=4L),value=unlist(d[tr],use.names=FALSE))
  stopifnot(nrow(long)==4L*n, all(table(long$unit)==4L))
  gllvmTMB::gllvmTMB(value ~ 0 + trait + trait:x + latent(0 + trait | unit,d=1,unique=FALSE),
                    data=long,trait="trait",unit="unit",family=gaussian(),REML=FALSE,engine="tmb",silent=TRUE)
}
extract_one <- function(fit) {
  if (backend == "drm") {
    b <- coef(fit); terms <- unlist(lapply(names(b), function(k) paste0(k,":",names(b[[k]]))))
    q <- setNames(unlist(b,use.names=FALSE),terms); V <- vcov(fit)
    list(q=q,se=sqrt(diag(V))[match(terms,rownames(V))])
  } else {
    td <- getS3method("tidy","gllvmTMB_multi",envir=asNamespace("gllvmTMB"))(fit,effects="fixed")
    list(q=setNames(td$estimate,td$term),se=setNames(td$std.error,td$term))
  }
}
validate_one <- function(fit,d) {
  stopifnot(!inherits(fit,"pigauto_mi_error"), identical(as.integer(fit$opt$convergence),0L))
  sdrep <- if (backend=="drm") fit$sdr else fit$sd_report
  stopifnot(isTRUE(sdrep$pdHess))
  ex <- extract_one(fit); stopifnot(setequal(names(ex$q),names(truth)),all(is.finite(ex$se)),all(ex$se>0))
  X <- cbind(1,d$x); inv <- solve(crossprod(X))
  if (backend=="drm") {
    ols <- solve(crossprod(X), crossprod(X,d$y))
    resid <- d$y - X%*%ols; sigma <- sqrt(sum(resid^2)/n)
    oracle_q <- setNames(c(ols,log(sigma)),names(truth))
    oracle_se <- c(sqrt(diag(inv)*sigma^2),sqrt(1/(2*n)))
  } else {
    Y <- as.matrix(d[paste0("y",1:4)]); ols <- solve(crossprod(X),crossprod(X,Y))
    # Marginal trait covariance uses fitted loading and Gaussian residual SD;
    # matrix calculation is independent of the package coefficient/SE extractor.
    Lambda <- fit$report$Lambda_B
    if (is.null(Lambda)) stop("No fitted Lambda_B available for covariance oracle")
    sigma <- fit$report$sigma_eps
    if (length(sigma)!=1L) stop("Expected scalar Gaussian observation SD")
    S <- Lambda%*%t(Lambda)+diag(sigma^2,4L)
    oracle_q <- setNames(c(ols[1,],ols[2,]),names(truth))
    oracle_se <- c(sqrt(diag(S)*inv[1,1]),sqrt(diag(S)*inv[2,2]))
  }
  c(coef=close_enough(ex$q[names(truth)],oracle_q,1e-5,"fit coefficients vs OLS"),
    se=close_enough(ex$se,oracle_se,1e-5,"fit SE vs Gaussian information"))
}
t0 <- proc.time()[[3L]]
fullfit <- fit_one(full); full_oracle <- validate_one(fullfit,full)
cat("FULL_DATA_OK",backend,"elapsed_s",proc.time()[[3L]]-t0,"\n"); flush.console()
cat("POSTERIOR_START",backend,"at",format(Sys.time(),tz="America/Edmonton"),"\n"); flush.console()
mi <- multi_impute(masked,tree,m=m,draws_method="posterior",log_transform=FALSE,
                   posterior_control=list(param_uncertainty="full",residual_prior="sep"),
                   seed=streams[["imputation"]],verbose=TRUE)
stopifnot(isTRUE(mi$posterior$converged),identical(mi$mi_workflow,"pigauto_posterior_mi_v1"))
cat("POSTERIOR_OK",backend,"wall",proc.time()[[3L]]-t0,"\n");flush.console()
fits <- with_imputations(mi,fit_one,.progress=FALSE,.on_error="continue")
stopifnot(length(fits)==m,!any(vapply(fits,inherits,logical(1),"pigauto_mi_error")))
errors <- t(vapply(seq_len(m),function(i) validate_one(fits[[i]],mi$datasets[[i]]),numeric(2)))
pooled <- pool_mi(fits)
ex <- lapply(fits,extract_one)
Q <- t(vapply(ex,function(z) z$q[names(truth)],numeric(length(truth))))
U <- t(vapply(ex,function(z) z$se^2,numeric(length(truth))))
qbar <- colMeans(Q); W <- colMeans(U); B <- colSums(sweep(Q,2,qbar)^2)/(m-1)
Tvar <- W+(1+1/m)*B; r <- (1+1/m)*B/W
nu <- (m-1)*(1+1/r)^2; crit <- qt(.975,nu)
oracle <- data.frame(term=names(truth),estimate=qbar,W=W,B=B,totalvar=Tvar,df=nu,
                     conf.low=qbar-crit*sqrt(Tvar),conf.high=qbar+crit*sqrt(Tvar),truth=truth)
idx <- match(names(truth),pooled$term)
designated <- if (backend=="drm") "mu:x" else "traity1:x"
stopifnot(!anyNA(idx),B[[designated]]>0)
arith <- c(estimate=close_enough(pooled$estimate[idx],qbar,1e-8,"Rubin mean"),
           variance=close_enough(pooled$std.error[idx]^2,Tvar,1e-8,"Rubin variance"),
           df=close_enough(pooled$df[idx],nu,1e-8,"Rubin df"),
           low=close_enough(pooled$conf.low[idx],oracle$conf.low,1e-8,"Rubin lower"),
           high=close_enough(pooled$conf.high[idx],oracle$conf.high,1e-8,"Rubin upper"))
result <- list(backend=backend,seed=seed,streams=streams,n=n,m=m,truth=truth,full=full,
 masked=masked,tree=tree,mi=mi,fits=fits,fullfit=fullfit,pooled=pooled,oracle=oracle,
 oracle_errors=errors,full_oracle=full_oracle,arithmetic_errors=arith,
 wall_s=proc.time()[[3L]]-t0,versions=sapply(c("pigauto","drmTMB","gllvmTMB"),function(p) as.character(packageVersion(p))),
 session=capture.output(sessionInfo()))
path <- file.path(outdir,paste0(backend,"-",seed,".rds"));saveRDS(result,path)
reloaded <- readRDS(path);close_enough(pool_mi(reloaded$fits)$estimate,pooled$estimate,1e-8,"saved-fit pooling")
write.csv(oracle,file.path(outdir,paste0(backend,"-",seed,".csv")),row.names=FALSE)
cat("RECOVERY_SEED_OK",backend,seed,"wall_s",result$wall_s,"\n")
