# Fixed paired route comparison at n=300 under current baseline defaults.
args <- commandArgs(TRUE)
seed_count <- as.integer(args[[1L]])
out <- args[[2L]]
library(pigauto)
rows <- list(); k <- 0L; t0 <- proc.time()[[3L]]
for (s in seq_len(seed_count)) for (lambda in c(.3,1)) {
  set.seed(20261000L+s); tree <- ape::rtree(300L)
  R <- cov2cor(ape::vcv(tree)); V <- lambda*R+(1-lambda)*diag(300L)
  L <- chol(V+diag(1e-10,300L))
  Z <- t(L)%*%matrix(rnorm(300L*5L),300L,5L)
  u <- .7*Z[,1L]+sqrt(1-.7^2)*Z[,2L]
  catliab <- cbind(.6*Z[,1L]+Z[,3L], -.6*Z[,1L]+Z[,4L],Z[,5L])
  df <- data.frame(a=Z[,1L],b=u,
    binary=factor(ifelse(u>0,"yes","no"),levels=c("no","yes")),
    categorical=factor(max.col(catliab,ties.method="first"),levels=1:3),row.names=tree$tip.label)
  pd <- preprocess_traits(df,tree,log_transform=FALSE)
  spl <- make_missing_splits(pd$X_scaled,seed=20262000L+s,trait_map=pd$trait_map)
  accuracy <- function(bl,tm) {
    ix <- unique(rbind(arrayInd(spl$val_idx,dim(pd$X_scaled)),arrayInd(spl$test_idx,dim(pd$X_scaled))))
    rr <- unique(ix[ix[,2L]%in%tm$latent_cols,1L])
    if (tm$type=="binary") {
      y <- pd$X_scaled[rr,tm$latent_cols]; p <- plogis(bl$mu[rr,tm$latent_cols]); pred <- as.numeric(p>.5)
      c(accuracy=mean(pred==y),brier=mean((p-y)^2),n=length(rr))
    } else {
      y <- max.col(pd$X_scaled[rr,tm$latent_cols,drop=FALSE],ties.method="first")
      lp <- bl$mu[rr,tm$latent_cols,drop=FALSE]; p <- exp(lp-apply(lp,1,max));p <- p/rowSums(p)
      pred <- max.col(p,ties.method="first"); yy <- matrix(0,length(rr),length(tm$latent_cols));yy[cbind(seq_along(y),y)]<-1
      c(accuracy=mean(pred==y),brier=mean(rowSums((p-yy)^2)),n=length(rr))
    }
  }
  arms <- lapply(c("exact","per_column"),function(route) fit_baseline(pd,tree,splits=spl,predict_method=route))
  for (tr in c("binary","categorical")) {
    e <- accuracy(arms[[1L]],pd$trait_map[[tr]]); p <- accuracy(arms[[2L]],pd$trait_map[[tr]])
    k<-k+1L;rows[[k]]<-data.frame(seed=s,lambda=lambda,trait=tr,n=300L,cells=e[["n"]],exact=e[["accuracy"]],per_column=p[["accuracy"]],difference=e[["accuracy"]]-p[["accuracy"]],brier_exact=e[["brier"]],brier_per_column=p[["brier"]])
  }
  cat("PAIR_OK",s,lambda,"elapsed_s",proc.time()[[3L]]-t0,"\n");flush.console()
}
res <- do.call(rbind,rows);write.csv(res,out,row.names=FALSE)
summary <- do.call(rbind,lapply(split(res,list(res$lambda,res$trait)),function(d) {
  se <- if(nrow(d)>1)sd(d$difference)/sqrt(nrow(d)) else NA_real_
  crit <- if(nrow(d)>1)qt(.975,nrow(d)-1) else NA_real_
  data.frame(lambda=d$lambda[1],trait=d$trait[1],pairs=nrow(d),exact=mean(d$exact),per_column=mean(d$per_column),difference=mean(d$difference),mcse=se,low=mean(d$difference)-crit*se,high=mean(d$difference)+crit*se,brier_difference=mean(d$brier_exact-d$brier_per_column))
}));write.csv(summary,sub(".csv$","-summary.csv",out),row.names=FALSE)
cat("PREDICTIVE_PAIRS_OK",seed_count,"wall_s",proc.time()[[3L]]-t0,"\n")
