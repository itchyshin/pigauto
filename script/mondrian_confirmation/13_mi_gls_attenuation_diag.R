# Diagnostic (2026-09-23): why pigauto conformal MI halves a phylogenetic GLS slope.
# One tree, n = 400, bivariate BM (rho = 0.7), 120 x cells missing completely at random,
# epochs 150. Compares point imputation, conformal MI draws (OLS vs GLS), BM-correlated
# noise of equal size, and an oracle proper imputation under the true model.
# Run from the repo root: OPENBLAS_NUM_THREADS=1 Rscript script/mondrian_confirmation/13_mi_gls_attenuation_diag.R

suppressMessages({devtools::load_all(quiet = TRUE); library(ape)})
set.seed(20260823); n <- 400
tree <- rtree(n); V <- vcv(tree); V <- V / max(V); L <- chol(V + 1e-8 * diag(n))
Sig <- matrix(0.7, 2, 2); diag(Sig) <- 1
Z <- t(L) %*% matrix(rnorm(n * 2), n, 2) %*% chol(Sig)
truth <- data.frame(row.names = tree$tip.label, x = Z[, 1], y = Z[, 2])
miss <- sample(n, 120); df <- truth; df$x[miss] <- NA
t0 <- Sys.time()
res <- impute(df, tree, epochs = 150, verbose = FALSE, seed = 1)
xi <- res$completed$x[miss]
cat("cor(truth x, y) among missing:", round(cor(truth$x[miss], truth$y[miss]), 3), "\n")
cat("cor(imputed x, y) among missing:", round(cor(xi, truth$y[miss]), 3), "\n")
cat("cor(imputed x, truth x):", round(cor(xi, truth$x[miss]), 3), "\n")
cat("var ratio imputed/true (missing):", round(var(xi) / var(truth$x[miss]), 3), "\n")
cat("baseline path:", paste(names(res$baseline %||% list()), collapse = ","), "\n")
d <- res$completed; d$species <- rownames(d)
f <- nlme::gls(y ~ x, correlation = corBrownian(phy = tree, form = ~species), data = d, method = "ML")
cat("single-imputation gls slope:", round(coef(f)[2], 3), " (truth 0.70)\n")
cat("elapsed", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
mi <- multi_impute(df, tree, m = 5, draws_method = "conformal", epochs = 150, verbose = FALSE, seed = 1)
D <- sapply(mi$datasets, function(d) d$x[miss])
cat("draw SD across m (mean over cells):", round(mean(apply(D, 1, sd)), 3), "\n")
cat("residual SD of point imputation:", round(sd(xi - truth$x[miss]), 3), "\n")
cat("var(one draw)/var(truth) among missing:", round(var(D[, 1]) / var(truth$x[miss]), 3), "\n")
cat("cor(one draw, y):", round(cor(D[, 1], truth$y[miss]), 3), "\n")
sl <- sapply(mi$datasets, function(d) { d$species <- rownames(d); coef(nlme::gls(y ~ x, correlation = corBrownian(phy = tree, form = ~species), data = d, method = "ML"))[2] })
cat("per-draw gls slopes:", round(sl, 3), "\n")
cat("conformal score x:", round(res$fit$conformal_scores["x"], 3), " trait sd:", round(sd(df$x, na.rm = TRUE), 3), "\n")
d1 <- mi$datasets[[1]]
cat("rownames match tips:", identical(rownames(d1), rownames(df)), " head:", head(rownames(d1), 3), "\n")
cat("observed x unchanged:", isTRUE(all.equal(d1$x[-miss], df$x[-miss])), " max abs diff:", signif(max(abs(d1$x[-miss] - df$x[-miss])), 3), "\n")
cat("y unchanged:", isTRUE(all.equal(d1$y, df$y)), " max abs diff:", signif(max(abs(d1$y - df$y)), 3), "\n")
cat("cor(d1$y, truth y):", round(cor(d1$y, df$y), 3), "\n")
ols_draw <- sapply(mi$datasets, function(d) coef(lm(y ~ x, data = d))[2])
cat("per-draw OLS slopes:", round(ols_draw, 3), "\n")
cat("OLS slope on truth:", round(coef(lm(y ~ x, data = truth))[2], 3), " on single imputation:", round(coef(lm(y ~ x, data = res$completed))[2], 3), "\n")
# correlated (BM-conditional) noise of the same marginal size, as a contrast
set.seed(2); Vm <- V[miss, miss] - V[miss, -miss] %*% solve(V[-miss, -miss], V[-miss, miss])
e <- as.numeric(t(chol(Vm + 1e-8 * diag(length(miss)))) %*% rnorm(length(miss))); e <- e / sd(e) * sd(D[, 1] - xi)
dc <- res$completed; dc$x[miss] <- xi + e; dc$species <- rownames(dc)
cat("gls slope, BM-correlated noise of equal size:", round(coef(nlme::gls(y ~ x, correlation = corBrownian(phy = tree, form = ~species), data = dc, method = "ML"))[2], 3), "\n")
# Oracle proper imputation: (x, y) ~ N(0, Sig %x% V) up to trait scale; condition x_mis on x_obs and all y
s <- c(sd(truth$x), sd(truth$y)); S <- diag(s) %*% Sig %*% diag(s)
C <- kronecker(S, V)                        # order: x(1..n), y(1..n)
im <- miss; io <- c(setdiff(1:n, miss), n + 1:n)
obs <- c(truth$x[-miss], truth$y)
mu_c <- C[im, io] %*% solve(C[io, io], obs)
V_c <- C[im, im] - C[im, io] %*% solve(C[io, io], C[io, im])
Lc <- t(chol(V_c + 1e-9 * diag(length(im))))
or <- sapply(1:5, function(k) { d <- truth; d$x[miss] <- as.numeric(mu_c + Lc %*% rnorm(length(im))); d$species <- rownames(d)
  coef(nlme::gls(y ~ x, correlation = corBrownian(phy = tree, form = ~species), data = d, method = "ML"))[2] })
cat("oracle proper-MI gls slopes:", round(or, 3), "\n")
cat("oracle conditional SD (mean):", round(mean(sqrt(diag(V_c))), 3), " vs pigauto draw SD:", round(mean(apply(D, 1, sd)), 3), "\n")
cat("gls slope on truth:", round(coef(nlme::gls(y ~ x, correlation = corBrownian(phy = tree, form = ~species), data = transform(truth, species = rownames(truth)), method = "ML"))[2], 3), "\n")
