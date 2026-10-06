# Screen 2: as dlam_compare.R, with the "on" runs from out_on2 (ordinal BM candidate also estimates lambda),
# plus continuous zRMSE and conformal coverage (c1, c2, cnt, prp) and discrete Brier, to price the gate change.
setwd("~/pigauto_sim")
get <- function(path) { r <- readRDS(path)$results; r[r$metric %in% c("accuracy", "zRMSE", "brier") , c("arm", "trait", "metric", "value", "coverage")] }
rows <- list()
for (f in list.files("dlam/out_off", "\\.rds$")) {
  fb <- file.path("results/core_bace", f); fn <- file.path("dlam/out_on2", f)
  if (!file.exists(fb) || !file.exists(fn)) next
  m <- regmatches(f, regexec("_l([0-9.]+)_r([0-9.]+)_mcar0.3_n([0-9]+)_s([0-9]+)", f))[[1]]
  off <- get(file.path("dlam/out_off", f)); on <- get(fn); b <- get(fb)
  off$arm <- paste0(off$arm, "|off"); on$arm <- paste0(on$arm, "|on"); b <- b[b$arm == "bace", ]
  rows[[f]] <- cbind(lambda = as.numeric(m[2]), rho = as.numeric(m[3]), n = as.integer(m[4]), seed = as.integer(m[5]),
                     rbind(off, on[on$arm != "floor|on", ], b))
}
x <- do.call(rbind, rows)
lab <- c("floor|off" = "floor", "gnn_off|off" = "gate,fixed1", "gnn_off_pure|off" = "nogate,fixed1",
         "gnn_off|on" = "gate,est", "gnn_off_pure|on" = "nogate,est", "bace" = "BACE")
x$arm <- lab[x$arm]; x <- x[!is.na(x$arm), ]
tab <- function(d, val = "value", by = c("n", "lambda")) {
  a <- aggregate(as.formula(paste(val, "~", paste(c(by, "arm"), collapse = "+"))), d, mean)
  w <- reshape(a, idvar = by, timevar = "arm", direction = "wide"); names(w) <- sub(paste0(val, "."), "", names(w), fixed = TRUE)
  w[do.call(order, w[by]), c(by, intersect(unname(lab), names(w)))] }
disc <- x[x$metric == "accuracy" & x$trait %in% c("bin", "ord", "cat3"), ]
cat("Discrete accuracy by trait (pooled over rho):\n"); print(tab(disc, by = c("n", "lambda", "trait")), digits = 3, row.names = FALSE)
cat("\nDiscrete accuracy pooled over traits and rho:\n"); print(tab(disc), digits = 3, row.names = FALSE)
k <- aggregate(value ~ n + lambda + rho + seed + arm, disc, mean)
kw <- reshape(k, idvar = c("n", "lambda", "rho", "seed"), timevar = "arm", direction = "wide")
cat("\nPaired difference from BACE (pooled traits), mean (SE):\n")
for (a in c("gate,fixed1", "nogate,est")) { kw$d <- kw[[paste0("value.", a)]] - kw$value.BACE
  s <- aggregate(d ~ n + lambda, kw, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
  cat(a, ":", paste0("n", s$n, " l", s$lambda, " ", sprintf("%+.3f (%.3f)", s$d[, "m"], s$d[, "se"]), collapse = "; "), "\n") }
cat("\nDiscrete Brier (lower is better):\n"); print(tab(x[x$metric == "brier" & x$trait %in% c("bin", "ord", "cat3"), ]), digits = 3, row.names = FALSE)
cont <- x[x$metric == "zRMSE" & x$trait %in% c("c1", "c2", "cnt", "prp"), ]
cat("\nContinuous zRMSE (lower is better), pooled over c1 c2 cnt prp and rho:\n"); print(tab(cont), digits = 3, row.names = FALSE)
cc <- cont[cont$arm != "BACE" & cont$arm != "floor" & is.finite(cont$coverage), ]
cat("\nContinuous conformal coverage (nominal 0.95):\n"); print(tab(cc, val = "coverage"), digits = 3, row.names = FALSE)
saveRDS(x, "dlam/screen2_long.rds")
