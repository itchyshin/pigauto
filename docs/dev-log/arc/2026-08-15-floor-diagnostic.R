suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260429L); n <- 80L; tree <- ape::rtree(n); sp <- tree$tip.label
bm_draw <- function(seed) withr::with_seed(seed, {
  v <- as.numeric(ape::rTraitCont(tree, model="BM", sigma=1)); names(v) <- tree$tip.label; v[sp] })
v1<-bm_draw(11L); v2<-bm_draw(12L); v3<-bm_draw(13L); v4<-bm_draw(14L)
qs <- stats::quantile(v3, c(0,1/3,2/3,1), na.rm=TRUE); qs[1] <- qs[1]-1e-9
cat3 <- factor(c("a","b","c")[as.integer(cut(v3, qs, include.lowest=TRUE))], levels=c("a","b","c"))
df <- data.frame(row.names=sp, x_cont=v1,
                 x_bin=factor(ifelse(v2>0,"yes","no"), levels=c("no","yes")),
                 x_cat=cat3, x_cnt=as.integer(round(pmax(v4+5,0))))
for (r in 1:3) {
  msgs <- character(0)
  res <- withCallingHandlers(
    pigauto::impute(df, tree, epochs=30L, eval_every=10L, patience=5L,
                    missing_frac=0.30, verbose=TRUE, seed=20260429L),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  fit <- res$fit; lc <- fit$trait_map[[1]]$latent_cols[1]
  fired <- grep("Post-refine floor", msgs, value = TRUE)
  # measure B on the val+test-hidden surface (the internal one)
  splits <- res$splits; X <- res$data$X_scaled
  rows <- which(matrix(replace(logical(length(X)), splits$val_idx, TRUE), nrow(X), ncol(X))[, lc])
  truth <- X[rows, lc]; bm <- fit$baseline$mu[rows, lc]
  B2 <- predict(fit, return_se=FALSE,
                .mask_observed_idx = c(splits$val_idx, splits$test_idx))$imputed_latent[rows, lc]
  cat(sprintf("rep%d gate=%.2f  B2(val+test hidden)=%.4f  bm=%.4f  ratio=%.4f  floor_fired=%s\n",
              r, fit$r_cal_gnn[lc], mean((B2-truth)^2), mean((bm-truth)^2),
              mean((B2-truth)^2)/mean((bm-truth)^2),
              if (length(fired)) "YES" else "no"))
  if (length(fired)) cat("   ", fired[1], "\n")
}
