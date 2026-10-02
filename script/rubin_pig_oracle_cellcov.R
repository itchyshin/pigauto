here <- "script"; suppressMessages({ source(file.path(here, "campaign_gnn_off_lib.R")); source(file.path(here, "rubin_lib.R")) })
src <- readLines(file.path(here, "rubin_pig_oracle.R")); eval(parse(text = src[grep("^oracle_sets <- function", src):(grep("^rows <- list[(][)]", src) - 1)]))
s2 <- readLines(file.path(here, "rubin_pig_prior_trial_summary.R")); eval(parse(text = s2[grep("^cell_cov <- function", s2):(grep("^files <- ", s2) - 1)]))
RNGkind("L'Ecuyer-CMRG")
for (l in c("1", "0.7")) { v <- c()
  for (f in Sys.glob(sprintf(".unlazy/pigauto-rubin-arm/prior_trial/prior_base_n1000_l%s_r0.5_s*.rds", l))) { x <- readRDS(f); set.seed(x$seed + 4242L); v <- c(v, cell_cov(oracle_sets(x), x)) }
  cat("lambda", l, "oracle 20-draw cell coverage:", round(mean(v), 3), "over", length(v), "datasets\n") }
