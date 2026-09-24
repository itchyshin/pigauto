# Real-data check for the lambda default (lane feat/joint-lambda-default, brain D-278).
#
# One job = (dataset, size, seed, lambda_mode). For each job: subsample `size` species (seeded),
# mask 20% of the observed cells of every continuous trait MCAR, impute with gnn = FALSE under the
# given lambda_mode, and score z-RMSE on the masked cells (RMSE / SD of the observed truth, on the log
# scale for positive skewed traits). Discrete traits present in the data are kept in the fit so the
# mixed-type joint path runs; their masked-cell accuracy is scored too. Writes one small rds per job.
#
# Usage:
#   Rscript script/bench_lambda_datasets.R --list                     # print the job grid
#   Rscript script/bench_lambda_datasets.R --job <i> --out <dir> [--data <snapshot dir>]
#   Rscript script/bench_lambda_datasets.R --collect <dir>            # summary tables
suppressMessages({ library(pigauto); library(ape) })
args <- commandArgs(trailingOnly = TRUE)
arg <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1L] }
snap <- arg("--data", "useful/bace_data_snapshot/data")

specs <- list(
  avonet    = list(cont = c("mass_g", "wing_length_mm", "beak_length_culmen_mm", "tarsus_length_mm",
                            "tail_length_mm"),
                   log = c("mass_g", "wing_length_mm", "beak_length_culmen_mm", "tarsus_length_mm",
                           "tail_length_mm"),
                   disc = c("trophic_level", "migration")),
  amphibio  = list(cont = c("body_size_mm", "body_mass_g"), log = c("body_size_mm", "body_mass_g"),
                   disc = c("habitat")),
  bien      = list(cont = c("height_m", "leaf_area", "sla", "seed_mass", "wood_density"),
                   log = c("height_m", "leaf_area", "sla", "seed_mass"), disc = character(0)),
  globtherm = list(cont = c("Tmax", "Tmin", "elevation_max"), log = character(0), disc = character(0)),
  leptraits = list(cont = c("wingspan_lower", "flight_duration", "n_hostplant_families"),
                   log = c("wingspan_lower", "flight_duration", "n_hostplant_families"),
                   disc = character(0)),
  pantheria = list(cont = c("body_mass_g", "head_body_length_mm", "gestation_d", "max_longevity_m"),
                   log = c("body_mass_g", "head_body_length_mm", "gestation_d", "max_longevity_m"),
                   disc = c("terrestriality"))
)
sizes <- c(small = 300L, large = 2000L)
modes <- c("fixed_1", "estimate", "cv", "bayes")
grid <- rbind(
  expand.grid(dataset = c(names(specs)), size = names(sizes), seed = 1:5, mode = modes,
              stringsAsFactors = FALSE),
  expand.grid(dataset = "avonet300", size = "bundled", seed = 1:5, mode = modes, stringsAsFactors = FALSE))
grid <- grid[!(grid$dataset == "globtherm" & grid$size == "large"), ]   # globtherm has 1969 species: "large" = all
grid$size[grid$dataset == "globtherm" & grid$size == "small"] <- "small"
grid <- rbind(grid, expand.grid(dataset = "globtherm", size = "full", seed = 1:5, mode = modes,
                                stringsAsFactors = FALSE))
rownames(grid) <- NULL

load_data <- function(nm) {
  if (nm == "avonet300") {
    utils::data("avonet300", "tree300", package = "pigauto", envir = environment())
    df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
    return(list(df = df, tree = tree300,
                spec = list(cont = c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length"),
                            log = c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length"),
                            disc = c("Trophic.Level", "Primary.Lifestyle", "Migration"))))
  }
  e <- new.env()
  load(file.path(snap, paste0(nm, "_traits.rda")), e); load(file.path(snap, paste0(nm, "_tree.rda")), e)
  obs <- mget(ls(e), e)
  tree <- Filter(function(x) inherits(x, "phylo"), obs)[[1]]
  df <- Filter(is.data.frame, obs)[[1]]
  if (is.null(rownames(df)) || all(rownames(df) == seq_len(nrow(df)))) {
    sc <- intersect(c("species", "Species", "tip_label"), colnames(df))
    if (length(sc)) { rownames(df) <- df[[sc[1]]] }
  }
  list(df = df, tree = tree, spec = specs[[nm]])
}

if ("--grid" %in% args) { utils::write.csv(cbind(job = seq_len(nrow(grid)), grid), stdout(), row.names = FALSE); quit(save = "no") }
if ("--list" %in% args) {
  cat("jobs:", nrow(grid), "\n"); print(table(grid$dataset, grid$size)); quit(save = "no")
}

if ("--collect" %in% args) {
  dir <- arg("--collect")
  rs <- lapply(list.files(dir, "\\.rds$", full.names = TRUE), readRDS)
  ok <- Filter(function(r) is.null(r$error), rs)
  cat("jobs read:", length(rs), " ok:", length(ok), " failed:", length(rs) - length(ok), "\n")
  if (length(rs) > length(ok)) print(table(sapply(Filter(function(r) !is.null(r$error), rs), `[[`, "error")))
  z <- do.call(rbind, lapply(ok, function(r) r$scores))
  saveRDS(z, file.path(dirname(normalizePath(dir)), paste0(basename(dir), "_scores.rds")))
  cont <- z[z$kind == "continuous", ]
  a <- aggregate(value ~ dataset + size + n + mode, cont, mean)
  w <- reshape(a, idvar = c("dataset", "size", "n"), timevar = "mode", direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w$est_vs_fixed_pct <- round(100 * (w$estimate - w$fixed_1) / w$fixed_1, 1)
  w$best <- modes[apply(w[, modes], 1, which.min)]
  w <- w[order(w$dataset, w$n), ]
  cat("\n== mean z-RMSE over continuous traits (masked cells, 5 seeds), gnn = FALSE ==\n")
  print(w, digits = 3, row.names = FALSE)
  lam <- aggregate(lambda_hat ~ dataset + size + trait, cont[cont$mode == "estimate", ], mean)
  cat("\n== estimated lambda per trait (mean over seeds) ==\n"); print(lam, digits = 2, row.names = FALSE)
  pt <- aggregate(value ~ dataset + size + trait + mode, cont, mean)
  pt <- reshape(pt, idvar = c("dataset", "size", "trait"), timevar = "mode", direction = "wide")
  names(pt) <- sub("^value\\.", "", names(pt))
  pt$est_vs_fixed_pct <- round(100 * (pt$estimate - pt$fixed_1) / pt$fixed_1, 1)
  cat("\n== per trait z-RMSE ==\n"); print(pt[order(pt$dataset, pt$size, pt$trait), ], digits = 3, row.names = FALSE)
  d <- z[z$kind == "discrete", ]
  if (nrow(d)) {
    da <- aggregate(value ~ dataset + size + trait + mode, d, mean)
    da <- reshape(da, idvar = c("dataset", "size", "trait"), timevar = "mode", direction = "wide")
    names(da) <- sub("^value\\.", "", names(da))
    cat("\n== discrete accuracy (masked cells) ==\n"); print(da, digits = 3, row.names = FALSE)
  }
  tm <- aggregate(secs ~ dataset + size + mode, do.call(rbind, lapply(ok, function(r) r$timing)), median)
  tm <- reshape(tm, idvar = c("dataset", "size"), timevar = "mode", direction = "wide")
  names(tm) <- sub("^secs\\.", "", names(tm))
  cat("\n== median seconds per fit ==\n"); print(tm, digits = 3, row.names = FALSE)
  write.csv(w, file.path(dirname(normalizePath(dir)), paste0(basename(dir), "_summary.csv")), row.names = FALSE)
  write.csv(pt, file.path(dirname(normalizePath(dir)), paste0(basename(dir), "_per_trait.csv")), row.names = FALSE)
  quit(save = "no")
}

i <- as.integer(arg("--job")); out <- arg("--out"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
job <- grid[i, ]
file <- file.path(out, sprintf("%s_%s_s%d_%s.rds", job$dataset, job$size, job$seed, job$mode))
if (file.exists(file)) quit(save = "no")
res <- tryCatch({
  d <- load_data(job$dataset); spec <- d$spec
  keep <- intersect(c(spec$cont, spec$disc), colnames(d$df))
  df <- d$df[intersect(rownames(d$df), d$tree$tip.label), keep, drop = FALSE]
  for (v in spec$log) if (v %in% names(df)) { x <- df[[v]]; x[!is.finite(x) | x <= 0] <- NA; df[[v]] <- x }
  set.seed(1000L * job$seed + nchar(job$dataset))
  n_target <- switch(job$size, small = 300L, large = 2000L, full = nrow(df), bundled = nrow(df))
  # species with at least one observed continuous value, so every subsample carries signal
  has <- rowSums(!is.na(df[intersect(spec$cont, names(df))])) > 0
  pool <- rownames(df)[has]
  spp <- if (length(pool) > n_target) sample(pool, n_target) else pool
  df <- df[spp, , drop = FALSE]; tree <- ape::keep.tip(d$tree, spp)
  for (v in spec$disc) if (v %in% names(df)) df[[v]] <- droplevels(df[[v]])
  truth <- df; mask <- matrix(FALSE, nrow(df), ncol(df), dimnames = list(NULL, names(df)))
  for (v in names(df)) {
    o <- which(!is.na(df[[v]])); if (length(o) < 20L) next
    m <- sample(o, max(1L, round(0.2 * length(o)))); mask[m, v] <- TRUE; df[m, v] <- NA
  }
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressWarnings(impute(df, tree, gnn = FALSE, lambda_mode = job$mode))
  secs <- proc.time()[["elapsed"]] - t0
  lam <- fit$fit$model_config$lambda_per_trait
  sc <- list()
  for (v in names(df)) {
    sel <- mask[, v]; if (!any(sel)) next
    if (v %in% spec$cont) {
      tv <- truth[[v]]; pv <- fit$completed[[v]]
      if (v %in% spec$log) { tv <- log(tv); pv <- log(pmax(pv, .Machine$double.eps)) }
      val <- sqrt(mean((pv[sel] - tv[sel])^2)) / stats::sd(tv[!is.na(tv)])
      lh <- if (!is.null(lam) && v %in% names(lam)) unname(lam[[v]]) else NA_real_
      sc[[v]] <- data.frame(trait = v, kind = "continuous", value = val, lambda_hat = lh)
    } else {
      val <- mean(as.character(fit$completed[[v]][sel]) == as.character(truth[[v]][sel]))
      sc[[v]] <- data.frame(trait = v, kind = "discrete", value = val, lambda_hat = NA_real_)
    }
  }
  sc <- do.call(rbind, sc)
  sc <- cbind(dataset = job$dataset, size = job$size, n = nrow(df), seed = job$seed, mode = job$mode, sc)
  list(scores = sc, timing = data.frame(dataset = job$dataset, size = job$size, mode = job$mode, secs = secs))
}, error = function(e) list(error = conditionMessage(e), job = job))
saveRDS(res, file)
