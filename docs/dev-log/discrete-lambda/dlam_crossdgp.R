# Screen 4 (cross-DGP): pigauto variants on the four-arm factorial cells (paired with stored BACE where present) and on
# bace_dgp, bm_mixed, ou_mixed, AVONET300. Variants: default (gate + floor, discrete lambda off) = "default";
# floor only + discrete lambda = "floor,est"; neither + lambda off = "none,l1"; neither + discrete lambda = "none,est".
# Question for the gate verdict: is "none,est" (or "floor,est") ever worse than "default" by more than 2 MCSE?
setwd("~/pigauto_sim")
get <- function(p) { r <- readRDS(p)$results; r[r$metric %in% c("accuracy", "zRMSE", "brier"), c("arm", "trait", "metric", "value", "coverage")] }
disc_traits <- c("bin", "ord", "cat3", "y", "Trophic.Level", "Primary.Lifestyle", "Migration", "x_bin", "b1", "d_bin", "z_cat")
rows <- list(); bad <- 0
for (f in list.files("dlam/a2_off", "\\.rds$")) {
  fo <- file.path("dlam/a2_on", f); if (!file.exists(fo)) { bad <- bad + 1; next }
  a <- tryCatch(get(file.path("dlam/a2_off", f)), error = function(e) NULL); b <- tryCatch(get(fo), error = function(e) NULL)
  if (is.null(a) || is.null(b)) { bad <- bad + 1; next }
  a <- a[a$arm %in% c("gnn_off", "gnn_off_pure", "floor"), ]; a$arm <- c(gnn_off = "default", gnn_off_pure = "none,l1", floor = "mode")[a$arm]
  b <- b[b$arm %in% c("gnn_off_nogate", "gnn_off_pure"), ]; b$arm <- c(gnn_off_nogate = "floor,est", gnn_off_pure = "none,est")[b$arm]
  fb <- file.path("results/factorial_bace", f)
  bb <- if (file.exists(fb)) { z <- get(fb); z <- z[z$arm == "bace", ]; if (nrow(z)) z$arm <- "BACE"; z } else NULL
  dgp <- sub("_(BM|OU)_l.*", "", f); evo <- regmatches(f, regexpr("_(BM|OU)_", f)); evo <- gsub("_", "", evo)
  m <- regmatches(f, regexec("_l([0-9.]+)_r([0-9.]+)_([a-z]+)([0-9.]+)_n([0-9]+)_s([0-9]+)", f))[[1]]
  grp <- if (dgp == "types_mixed") sprintf("tm %s %s%s l%s", evo, m[4], m[5], m[2]) else dgp
  rows[[f]] <- cbind(group = grp, n = as.integer(m[6]), rho = m[3], seed = as.integer(m[7]), rbind(a, b, bb))
}
x <- do.call(rbind, rows)
cat("files paired off/on:", length(rows), " unusable:", bad, "\n")
x$kind <- ifelse(x$metric == "zRMSE", "cont", ifelse(x$metric %in% c("accuracy", "brier"), "disc", NA))
lab <- c("mode", "default", "floor,est", "none,l1", "none,est", "BACE")
# dataset-level means within a dataset (pooled over traits of a kind), then mean and MCSE over datasets
summ <- function(metric, kind, val = "value") {
  d <- x[x$metric == metric & x$kind == kind & is.finite(x[[val]]), ]
  d <- aggregate(as.formula(paste(val, "~ group + n + rho + seed + arm")), d, mean)
  w <- reshape(d, idvar = c("group", "n", "rho", "seed"), timevar = "arm", direction = "wide"); names(w) <- sub(paste0(val, "."), "", names(w), fixed = TRUE)
  out <- do.call(rbind, lapply(split(w, list(w$group, w$n), drop = TRUE), function(g) {
    r <- data.frame(group = g$group[1], n = g$n[1], k = nrow(g))
    for (a in intersect(lab, names(g))) r[[a]] <- mean(g[[a]], na.rm = TRUE)
    for (a in c("none,est", "floor,est")) if (all(c(a, "default") %in% names(g))) {
      dd <- g[[a]] - g$default; r[[paste0("d[", a, "]")]] <- mean(dd, na.rm = TRUE); r[[paste0("se[", a, "]")]] <- sd(dd, na.rm = TRUE) / sqrt(sum(is.finite(dd))) }
    r }))
  out[order(out$group, out$n), ] }
fmt <- function(t, title) { cat("\n", title, "\n", sep = ""); num <- vapply(t, is.numeric, TRUE); t[num] <- lapply(t[num], round, 3); print(t, row.names = FALSE) }
acc <- summ("accuracy", "disc"); fmt(acc, "Discrete accuracy (higher better); d[.] = variant - default, paired, se = MCSE")
bri <- summ("brier", "disc"); fmt(bri, "Discrete Brier (lower better)")
zr <- summ("zRMSE", "cont"); fmt(zr, "Continuous zRMSE (lower better)")
cv <- summ("zRMSE", "cont", val = "coverage"); fmt(cv[, intersect(c("group", "n", "k", "default", "floor,est", "none,l1", "none,est"), names(cv))], "Continuous conformal coverage (nominal 0.95)")
worse <- function(t, higher_better, a) { d <- t[[paste0("d[", a, "]")]]; s <- t[[paste0("se[", a, "]")]]
  idx <- if (higher_better) which(d < -2 * s) else which(d > 2 * s); if (length(idx)) paste(t$group[idx], paste0("n", t$n[idx]), sprintf("%+.3f", d[idx]), collapse = "; ") else "none" }
cat("\nCells where the variant is WORSE than default by > 2 MCSE:\n")
for (a in c("none,est", "floor,est")) {
  cat(a, " accuracy:", worse(acc, TRUE, a), "\n"); cat(a, " Brier:   ", worse(bri, FALSE, a), "\n"); cat(a, " zRMSE:   ", worse(zr, FALSE, a), "\n") }
ok <- length(rows) >= 0.95 * 3785
cat(if (ok) "\nCROSSDGP_OK\n" else "\nCROSSDGP_INCOMPLETE\n")
saveRDS(x, "dlam/screen4_long.rds")
