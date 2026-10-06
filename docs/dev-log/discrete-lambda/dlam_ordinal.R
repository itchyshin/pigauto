# Screen 5 (ordinal): cumulative "class >= k" decomposition (options pigauto.ordinal_method = "cumulative",
# with discrete lambda estimated) vs the current ordinal route, BACE and the mode floor, on the 1,197 paired core datasets.
setwd("~/pigauto_sim")
get <- function(p) { r <- readRDS(p)$results; r[r$metric %in% c("accuracy", "brier", "macroF1"), c("arm", "trait", "metric", "value")] }
rows <- list()
for (f in list.files("dlam/out_ord", "\\.rds$")) {
  fb <- file.path("results/core_bace", f); fo <- file.path("dlam/out_on", f); fs <- file.path("results/core", f)
  if (!file.exists(fb) || !file.exists(fo)) next
  m <- regmatches(f, regexec("_l([0-9.]+)_r([0-9.]+)_mcar0.3_n([0-9]+)_s([0-9]+)", f))[[1]]
  a <- get(file.path("dlam/out_ord", f)); a <- a[a$arm %in% c("gnn_off_pure", "gnn_off_nogate", "floor"), ]
  a$arm <- c(gnn_off_pure = "cumul,none", gnn_off_nogate = "cumul,floor", floor = "mode")[a$arm]
  b <- get(fo); b <- b[b$arm %in% c("gnn_off_pure", "gnn_off"), ]; b$arm <- c(gnn_off_pure = "route,none", gnn_off = "route,gate+floor")[b$arm]
  bb <- get(fb); bb <- bb[bb$arm == "bace", ]; bb$arm <- "BACE"
  fr <- if (file.exists(fs)) { z <- get(fs); z <- z[z$arm == "freq_lambda", ]; if (nrow(z)) z$arm <- "freq"; z } else NULL
  rows[[f]] <- cbind(lambda = as.numeric(m[2]), rho = as.numeric(m[3]), n = as.integer(m[4]), seed = as.integer(m[5]), rbind(a, b, bb, fr))
}
x <- do.call(rbind, rows)
np <- length(unique(paste(x$n, x$lambda, x$rho, x$seed)))
ord <- c("mode", "route,gate+floor", "route,none", "cumul,floor", "cumul,none", "BACE", "freq")
tab <- function(d, by = c("n", "lambda")) { a <- aggregate(as.formula(paste("value ~", paste(c(by, "arm"), collapse = "+"))), d, mean)
  w <- reshape(a, idvar = by, timevar = "arm", direction = "wide"); names(w) <- sub("value.", "", names(w), fixed = TRUE)
  w[do.call(order, w[by]), c(by, intersect(ord, names(w)))] }
for (met in c("accuracy", "macroF1", "brier")) { cat("\nOrdinal", met, "(pooled over rho):\n"); print(tab(x[x$trait == "ord" & x$metric == met, ]), digits = 3, row.names = FALSE) }
cat("\nAll discrete traits, accuracy (bin, cat3 unchanged by the ordinal option):\n"); print(tab(x[x$trait %in% c("bin", "ord", "cat3") & x$metric == "accuracy", ], by = c("n", "lambda", "trait")), digits = 3, row.names = FALSE)
k <- x[x$trait == "ord" & x$metric == "accuracy", ]
kw <- reshape(k[, c("n", "lambda", "rho", "seed", "arm", "value")], idvar = c("n", "lambda", "rho", "seed"), timevar = "arm", direction = "wide")
cat("\nOrdinal accuracy, paired difference, mean (SE):\n")
for (a in c("cumul,none", "route,none")) { kw$d <- kw[[paste0("value.", a)]] - kw$value.BACE
  s <- aggregate(d ~ n + lambda, kw, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
  cat(a, "- BACE:", paste0("n", s$n, " l", s$lambda, " ", sprintf("%+.3f (%.3f)", s$d[, "m"], s$d[, "se"]), collapse = "; "), "\n") }
kw$d <- kw[["value.cumul,none"]] - kw[["value.route,none"]]
s <- aggregate(d ~ n + lambda, kw, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
cat("cumul,none - route,none:", paste0("n", s$n, " l", s$lambda, " ", sprintf("%+.3f (%.3f)", s$d[, "m"], s$d[, "se"]), collapse = "; "), "\n")
cat("\nORD_SCREEN", np, "\n")
