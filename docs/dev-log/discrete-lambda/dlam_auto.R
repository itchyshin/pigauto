# Screen 7: options(pigauto.discrete_lambda = "auto") (per binary/categorical trait, lambda fixed vs estimated chosen
# by validation Brier; branch feat/discrete-lambda-auto, 6969a20), gate and floor off ("none,auto"), against the
# earlier runs on the same datasets: "default" (gate + floor, lambda not estimated for discrete), "none,est" and BACE.
setwd("~/pigauto_sim")
get <- function(p) { r <- readRDS(p)$results; r[r$metric %in% c("accuracy", "zRMSE", "brier"), c("arm", "trait", "metric", "value", "coverage")] }
pairs <- list(c("dlam/out_off", "dlam/out_on"), c("dlam/a2_off", "dlam/a2_on"), c("dlam/c_off", "dlam/c_on"))
rows <- list(); miss <- 0
for (f in list.files("dlam/x_auto", "\\.rds$")) {
  src <- Filter(function(p) file.exists(file.path(p[1], f)) && file.exists(file.path(p[2], f)), pairs)
  if (!length(src)) { miss <- miss + 1; next }
  p <- src[[1]]
  a <- get(file.path(p[1], f)); a <- a[a$arm %in% c("gnn_off", "floor"), ]; a$arm <- c(gnn_off = "default", floor = "mode")[a$arm]
  b <- get(file.path(p[2], f)); b <- b[b$arm == "gnn_off_pure", ]; b$arm <- "none,est"
  x <- get(file.path("dlam/x_auto", f)); x$arm <- c(gnn_off_pure = "none,auto", gnn_off_nogate = "floor,auto")[x$arm]; x <- x[!is.na(x$arm), ]
  bb <- NULL; for (d in c("results/core_bace", "results/factorial_bace")) if (file.exists(file.path(d, f))) { z <- get(file.path(d, f)); bb <- z[z$arm == "bace", ]; if (nrow(bb)) bb$arm <- "BACE" }
  dgp <- sub("_(BM|OU)_l.*", "", f); evo <- gsub("_", "", regmatches(f, regexpr("_(BM|OU)_", f)))
  m <- regmatches(f, regexec("_l([0-9.]+)_r([0-9.]+)_([a-z]+)([0-9.]+)_n([0-9]+)_s([0-9]+)", f))[[1]]
  grp <- if (dgp == "types_mixed") sprintf("tm %s %s%s l%s", evo, m[4], m[5], m[2]) else dgp
  rows[[f]] <- cbind(group = grp, n = as.integer(m[6]), rho = m[3], seed = as.integer(m[7]), rbind(a, b, x, bb))
}
x <- do.call(rbind, rows); cat("datasets:", length(rows), " without an earlier run:", miss, "\n")
x$kind <- ifelse(x$metric == "zRMSE", "cont", "disc")
lab <- c("mode", "default", "none,est", "none,auto", "floor,auto", "BACE")
summ <- function(metric, kind) {
  d <- x[x$metric == metric & x$kind == kind & is.finite(x$value), ]
  d <- aggregate(value ~ group + n + rho + seed + arm, d, mean)
  w <- reshape(d, idvar = c("group", "n", "rho", "seed"), timevar = "arm", direction = "wide"); names(w) <- sub("value.", "", names(w), fixed = TRUE)
  out <- do.call(rbind, lapply(split(w, list(w$group, w$n), drop = TRUE), function(g) {
    r <- data.frame(group = g$group[1], n = g$n[1], k = nrow(g)); for (a in intersect(lab, names(g))) r[[a]] <- mean(g[[a]], na.rm = TRUE)
    for (a in c("none,auto", "none,est")) { dd <- g[[a]] - g$default; r[[paste0("d[", a, "]")]] <- mean(dd, na.rm = TRUE); r[[paste0("se[", a, "]")]] <- sd(dd, na.rm = TRUE) / sqrt(sum(is.finite(dd))) }
    r })); out[order(out$group, out$n), ] }
fmt <- function(t, title) { cat("\n", title, "\n", sep = ""); num <- vapply(t, is.numeric, TRUE); t[num] <- lapply(t[num], round, 3); print(t, row.names = FALSE) }
acc <- summ("accuracy", "disc"); bri <- summ("brier", "disc"); zr <- summ("zRMSE", "cont")
fmt(acc, "Discrete accuracy (higher better)"); fmt(bri, "Discrete Brier (lower better)"); fmt(zr, "Continuous zRMSE (lower better)")
worse <- function(t, hb, a) { d <- t[[paste0("d[", a, "]")]]; s <- t[[paste0("se[", a, "]")]]
  i <- if (hb) which(d < -2 * s) else which(d > 2 * s); if (length(i)) paste(t$group[i], paste0("n", t$n[i]), sprintf("%+.3f", d[i]), collapse = "; ") else "none" }
cat("\nWorse than default by > 2 MCSE:\n")
for (a in c("none,auto", "none,est")) { cat(a, "accuracy:", worse(acc, TRUE, a), "\n"); cat(a, "Brier:   ", worse(bri, FALSE, a), "\n"); cat(a, "zRMSE:   ", worse(zr, FALSE, a), "\n") }
cat("\nAUTO_ROWS", length(rows), "\n")
