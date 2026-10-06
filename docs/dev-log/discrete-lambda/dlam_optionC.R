# Option-C tests (Fisher review): (i) very weak signal, lambda = 0 and 0.1; (ii) OU continuous (ou_mixed, evo OU,
# replacing the duplicated run); (iii) all-discrete data (label-propagation path); (iv) n = 3000, lambda = 1.
# Same variant labels as screen 4. Paired per dataset; MCSE over datasets.
setwd("~/pigauto_sim")
get <- function(p) { r <- readRDS(p)$results; r[r$metric %in% c("accuracy", "zRMSE", "brier"), c("arm", "trait", "metric", "value", "coverage")] }
rows <- list()
for (f in list.files("dlam/c_off", "\\.rds$")) {
  fo <- file.path("dlam/c_on", f); if (!file.exists(fo)) next
  a <- get(file.path("dlam/c_off", f)); a <- a[a$arm %in% c("gnn_off", "gnn_off_pure", "floor"), ]; a$arm <- c(gnn_off = "default", gnn_off_pure = "none,l1", floor = "mode")[a$arm]
  b <- get(fo); b <- b[b$arm %in% c("gnn_off_nogate", "gnn_off_pure"), ]; b$arm <- c(gnn_off_nogate = "floor,est", gnn_off_pure = "none,est")[b$arm]
  dgp <- sub("_(BM|OU)_l.*", "", f); m <- regmatches(f, regexec("_l([0-9.]+)_r([0-9.]+)_[a-z]+[0-9.]+_n([0-9]+)_s([0-9]+)", f))[[1]]
  grp <- if (dgp == "types_mixed") sprintf("tm l%s", m[2]) else paste(dgp, "OU")
  rows[[f]] <- cbind(group = grp, n = as.integer(m[4]), rho = m[3], seed = as.integer(m[5]), rbind(a, b))
}
for (f in list.files("dlam/c_ad", "\\.rds$")) {
  r <- readRDS(file.path("dlam/c_ad", f)); a <- r$results; a <- a[a$metric %in% c("accuracy", "brier"), c("arm", "trait", "metric", "value", "coverage")]
  a$arm <- c(gnn_off = "default", gnn_off_nogate = "floor,est", gnn_off_pure = "none,est", floor = "mode")[a$arm]
  rows[[f]] <- cbind(group = sprintf("alldisc l%s", format(r$lambda)), n = r$n, rho = format(r$rho), seed = r$seed, a)
}
x <- do.call(rbind, rows); x$kind <- ifelse(x$metric == "zRMSE", "cont", "disc")
lab <- c("mode", "default", "floor,est", "none,l1", "none,est")
summ <- function(metric, kind, val = "value") {
  d <- x[x$metric == metric & x$kind == kind & is.finite(x[[val]]), ]; if (!nrow(d)) return(NULL)
  d <- aggregate(as.formula(paste(val, "~ group + n + rho + seed + arm")), d, mean)
  w <- reshape(d, idvar = c("group", "n", "rho", "seed"), timevar = "arm", direction = "wide"); names(w) <- sub(paste0(val, "."), "", names(w), fixed = TRUE)
  out <- do.call(rbind, lapply(split(w, list(w$group, w$n), drop = TRUE), function(g) {
    r <- data.frame(group = g$group[1], n = g$n[1], k = nrow(g)); for (a in intersect(lab, names(g))) r[[a]] <- mean(g[[a]], na.rm = TRUE)
    for (a in c("none,est", "floor,est")) if (all(c(a, "default") %in% names(g))) { dd <- g[[a]] - g$default
      r[[paste0("d[", a, "]")]] <- mean(dd, na.rm = TRUE); r[[paste0("se[", a, "]")]] <- sd(dd, na.rm = TRUE) / sqrt(sum(is.finite(dd))) }
    r })); out[order(out$group, out$n), ] }
fmt <- function(t, title) { if (is.null(t)) return(); cat("\n", title, "\n", sep = ""); num <- vapply(t, is.numeric, TRUE); t[num] <- lapply(t[num], round, 3); print(t, row.names = FALSE) }
fmt(summ("accuracy", "disc"), "Discrete accuracy (higher better)")
fmt(summ("brier", "disc"), "Discrete Brier (lower better)")
fmt(summ("zRMSE", "cont"), "Continuous zRMSE (lower better)")
cv <- summ("zRMSE", "cont", val = "coverage"); if (!is.null(cv)) fmt(cv[, intersect(c("group", "n", "k", "default", "floor,est", "none,l1", "none,est"), names(cv))], "Continuous conformal coverage")
cat("\nOPTIONC_ROWS", length(rows), "\n")
