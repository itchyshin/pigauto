# Screen 3 (attribution): which part of the safety machinery blocks discrete lambda?
# Variants on the same 1,197 four-arm core datasets, paired with BACE:
#   gate+floor (default gnn_off) | floor only (gnn_off_nogate) | neither (gnn_off_pure), each with discrete lambda off/on.
# Floor-off-gate-on is not a separate variant: with safety_floor = FALSE the gate falls back to the BM baseline.
setwd("~/pigauto_sim")
get <- function(p) { r <- readRDS(p)$results; r[r$metric %in% c("accuracy", "zRMSE", "brier"), c("arm", "trait", "metric", "value", "coverage")] }
rows <- list()
for (f in list.files("dlam/out_off", "\\.rds$")) {
  p <- c(off = "dlam/out_off", on = "dlam/out_on", ngoff = "dlam/out_ng_off", ngon = "dlam/out_ng_on")
  fb <- file.path("results/core_bace", f); if (!file.exists(fb) || !all(file.exists(file.path(p, f)))) next
  m <- regmatches(f, regexec("_l([0-9.]+)_r([0-9.]+)_mcar0.3_n([0-9]+)_s([0-9]+)", f))[[1]]
  a <- get(file.path(p["off"], f)); a <- a[a$arm %in% c("gnn_off", "gnn_off_pure", "floor"), ]; a$arm <- paste0(a$arm, "|off")
  b <- get(file.path(p["on"], f)); b <- b[b$arm %in% c("gnn_off", "gnn_off_pure"), ]; b$arm <- paste0(b$arm, "|on")
  c1 <- get(file.path(p["ngoff"], f)); c1 <- c1[c1$arm == "gnn_off_nogate", ]; c1$arm <- "gnn_off_nogate|off"
  c2 <- get(file.path(p["ngon"], f)); c2 <- c2[c2$arm == "gnn_off_nogate", ]; c2$arm <- "gnn_off_nogate|on"
  bb <- get(fb); bb <- bb[bb$arm == "bace", ]
  rows[[f]] <- cbind(lambda = as.numeric(m[2]), rho = as.numeric(m[3]), n = as.integer(m[4]), seed = as.integer(m[5]), rbind(a, b, c1, c2, bb))
}
x <- do.call(rbind, rows)
lab <- c("floor|off" = "mode", "gnn_off|off" = "gate+floor,l1", "gnn_off|on" = "gate+floor,est",
         "gnn_off_nogate|off" = "floor,l1", "gnn_off_nogate|on" = "floor,est",
         "gnn_off_pure|off" = "none,l1", "gnn_off_pure|on" = "none,est", "bace" = "BACE")
x$arm <- lab[x$arm]
np <- length(unique(paste(x$n, x$lambda, x$rho, x$seed)))
tab <- function(d, val = "value", by = c("n", "lambda")) {
  a <- aggregate(as.formula(paste(val, "~", paste(c(by, "arm"), collapse = "+"))), d, mean)
  w <- reshape(a, idvar = by, timevar = "arm", direction = "wide"); names(w) <- sub(paste0(val, "."), "", names(w), fixed = TRUE)
  w[do.call(order, w[by]), c(by, intersect(unname(lab), names(w)))] }
disc <- x[x$metric == "accuracy" & x$trait %in% c("bin", "ord", "cat3"), ]
cat("Discrete accuracy pooled over bin/ord/cat3 and rho:\n"); print(tab(disc), digits = 3, row.names = FALSE)
cat("\nBy trait:\n"); print(tab(disc, by = c("n", "lambda", "trait")), digits = 3, row.names = FALSE)
cat("\nDiscrete Brier:\n"); print(tab(x[x$metric == "brier" & x$trait %in% c("bin", "ord", "cat3"), ]), digits = 3, row.names = FALSE)
cont <- x[x$metric == "zRMSE" & x$trait %in% c("c1", "c2", "cnt", "prp"), ]
cat("\nContinuous zRMSE:\n"); print(tab(cont), digits = 3, row.names = FALSE)
cc <- cont[!(cont$arm %in% c("BACE", "mode")) & is.finite(cont$coverage), ]
cat("\nContinuous conformal coverage:\n"); print(tab(cc, val = "coverage"), digits = 3, row.names = FALSE)
k <- aggregate(value ~ n + lambda + rho + seed + arm, disc, mean)
kw <- reshape(k, idvar = c("n", "lambda", "rho", "seed"), timevar = "arm", direction = "wide")
cat("\nPaired difference from BACE, discrete accuracy pooled over traits, mean (SE):\n")
for (a in c("gate+floor,l1", "floor,est", "none,est")) { kw$d <- kw[[paste0("value.", a)]] - kw$value.BACE
  s <- aggregate(d ~ n + lambda, kw, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
  cat(a, ":", paste0("n", s$n, " l", s$lambda, " ", sprintf("%+.3f (%.3f)", s$d[, "m"], s$d[, "se"]), collapse = "; "), "\n") }
cat("\nPAIRED", np, "\n")
saveRDS(x, "dlam/screen3_long.rds")
