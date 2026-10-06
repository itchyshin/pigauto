# Discrete-lambda screen: pigauto (switch off / on; safety gate on = gnn_off, off = gnn_off_pure) vs BACE and the
# mode floor on the four-arm core datasets (types_mixed, BM, MCAR 30%), paired per dataset (file name).
setwd("~/pigauto_sim")
acc_of <- function(path) { r <- readRDS(path)$results; r <- r[r$metric == "accuracy" & r$trait %in% c("bin", "ord", "cat3"), ]
  r[, c("arm", "trait", "value")] }
rows <- list()
for (f in list.files("dlam/out_off", "\\.rds$")) {
  fb <- file.path("results/core_bace", f); fn <- file.path("dlam/out_on", f)
  if (!file.exists(fb) || !file.exists(fn)) next
  m <- regmatches(f, regexec("_l([0-9.]+)_r([0-9.]+)_mcar0.3_n([0-9]+)_s([0-9]+)", f))[[1]]
  off <- acc_of(file.path("dlam/out_off", f)); on <- acc_of(fn); b <- acc_of(fb)
  off$arm <- paste0(off$arm, "|off"); on$arm <- paste0(on$arm, "|on"); b <- b[b$arm == "bace", ]
  d <- rbind(off, on[on$arm != "floor|on", ], b)
  rows[[f]] <- cbind(lambda = as.numeric(m[2]), rho = as.numeric(m[3]), n = as.integer(m[4]), seed = as.integer(m[5]), d)
}
x <- do.call(rbind, rows)
arm_lab <- c("floor|off" = "floor", "gnn_off|off" = "gate,fixed1", "gnn_off_pure|off" = "nogate,fixed1",
             "gnn_off|on" = "gate,est", "gnn_off_pure|on" = "nogate,est", "bace" = "BACE")
x$arm <- arm_lab[x$arm]
cat("datasets paired with BACE:", length(unique(paste(x$n, x$lambda, x$rho, x$seed))), "\n")
# mean accuracy pooled over rho, by n, lambda, trait
w <- aggregate(value ~ n + lambda + trait + arm, x, mean)
ww <- reshape(w, idvar = c("n", "lambda", "trait"), timevar = "arm", direction = "wide"); names(ww) <- sub("value.", "", names(ww))
ww <- ww[order(ww$n, ww$lambda, ww$trait), c("n", "lambda", "trait", unname(arm_lab))]
print(ww, digits = 3, row.names = FALSE)
# pooled over the three trait types
p <- aggregate(value ~ n + lambda + arm, x, mean)
pp <- reshape(p, idvar = c("n", "lambda"), timevar = "arm", direction = "wide"); names(pp) <- sub("value.", "", names(pp))
cat("\nPooled over bin/ord/cat3 and rho:\n"); print(pp[order(pp$n, pp$lambda), c("n", "lambda", unname(arm_lab))], digits = 3, row.names = FALSE)
# paired difference vs BACE (dataset-level mean over traits), SE over datasets
k <- aggregate(value ~ n + lambda + rho + seed + arm, x, mean)
kw <- reshape(k, idvar = c("n", "lambda", "rho", "seed"), timevar = "arm", direction = "wide")
cat("\nPaired difference from BACE (pooled traits), mean (SE), per n x lambda:\n")
for (a in c("gate,fixed1", "nogate,est", "gate,est")) {
  kw$d <- kw[[paste0("value.", a)]] - kw$value.BACE
  s <- aggregate(d ~ n + lambda, kw, function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v))))
  cat(a, ":", paste0("n", s$n, " l", s$lambda, " ", sprintf("%+.3f (%.3f)", s$d[, "mean"], s$d[, "se"]), collapse = "; "), "\n")
}
saveRDS(x, "dlam/screen_long.rds")
