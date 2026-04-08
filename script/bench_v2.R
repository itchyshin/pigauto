suppressPackageStartupMessages({ devtools::load_all(quiet = TRUE) })

cat("\n=== pigauto benchmark v2: direct delta training ===\n\n")

bench <- simulate_benchmark(
  n_species = 150L,
  n_traits = 3L,
  scenarios = c("BM", "OU", "nonlinear", "regime_shift", "mixed"),
  missing_frac = 0.25,
  n_reps = 2L,
  epochs = 600L,
  verbose = TRUE
)
saveRDS(bench, file = "script/bench_v2.rds")

cat("\n\n=== Summary ===\n")
print(summary(bench))
