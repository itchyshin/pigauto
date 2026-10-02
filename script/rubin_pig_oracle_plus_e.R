# script/rubin_pig_oracle_plus_e.R: oracle MI plus pig_post's posterior Sigma_E noise on the imputed cells (mechanism test,
# docs/dev-log/arc/2026-10-02-rubin-pigauto-campaign.md). Run from the repo root after rubin_pig_cond_diag.R and
# rubin_pig_lambda_diag.R outputs are under .unlazy/pigauto-rubin-arm/.
here <- "script"
suppressMessages({ source(file.path(here, "campaign_gnn_off_lib.R")); source(file.path(here, "rubin_lib.R")) })
src <- readLines(file.path(here, "rubin_pig_oracle.R")); i0 <- grep("^oracle_sets <- function", src); i1 <- grep("^rows <- list\\(\\)", src) - 1
eval(parse(text = src[i0:i1]))      # oracle_sets(), score()
RNGkind("L'Ecuyer-CMRG")
dc <- ".unlazy/pigauto-rubin-arm/diag_cond"; dl <- ".unlazy/pigauto-rubin-arm/diag_lambda1"
rows <- list()
for (s in 1:6) {
  x <- readRDS(file.path(dc, sprintf("cond_diag_n1000_l1_r0.5_s%d.rds", s)))
  L <- readRDS(file.path(dl, sprintf("lambda_diag_n1000_l1_r0.5_s%d.rds", s)))
  SE <- L$params$Sigma_E; SEm <- apply(SE, c(1, 2), mean)          # posterior mean Sigma_E (latent z-scale of pigauto)
  # pigauto works on z-scored latents; put Sigma_E back on the data scale with the observed SDs of the block
  Yt <- as.matrix(x$truth); Yt[, "prp"] <- stats::qlogis(Yt[, "prp"]); Yo <- Yt; Yo[as.matrix(x$mask)] <- NA
  sdv <- apply(Yo, 2, sd, na.rm = TRUE); SEd <- SEm * outer(sdv, sdv)
  eig <- pagel_eigen(x$tree, rownames(x$truth)); comp <- est_pgls_slope_fast(x$truth, x$tree, eig = eig)
  set.seed(s + 4242L); os <- oracle_sets(x)
  set.seed(s + 777L)
  oe <- lapply(os, function(d) {
    E <- MASS::mvrnorm(nrow(d), rep(0, 4), SEd); colnames(E) <- x$block_traits
    Y <- as.matrix(d); Y[, "prp"] <- stats::qlogis(Y[, "prp"]); mk <- as.matrix(x$mask)
    Y[mk] <- Y[mk] + E[mk]; d2 <- as.data.frame(Y); d2$prp <- stats::plogis(d2$prp); d2 })
  r <- function(sets) score(sets, x, eig)[["estimate"]] - comp$estimate
  rows[[s]] <- data.frame(seed = s, SE_c1_data = SEd[1, 1], SE_c2_data = SEd[2, 2], SE_cor = SEd[1, 2] / sqrt(SEd[1, 1] * SEd[2, 2]),
                          oracle = r(os), oracle_plus_pigE = r(oe), pig_post = r(x$pig_post))
}
d <- do.call(rbind, rows); print(signif(d, 3), row.names = FALSE); cat("\nmeans:\n"); print(signif(colMeans(d[, -1]), 3))
