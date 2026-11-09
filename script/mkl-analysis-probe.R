pkgload::load_all('.', quiet = TRUE)
set.seed(1202L)
n <- 180L
z <- rnorm(n); x <- 0.4 * z + rnorm(n, sd = 0.8)
y <- 0.5 + 0.7 * x - 0.35 * z + rnorm(n, sd = 0.9)
data <- data.frame(y = y, x = x, z = z)
data$x[seq(3L, n, by = 4L)] <- NA_real_
run <- function(seed) {
  v <- multi_impute_analysis(data, y ~ x + z, 'x', model = 'lm', m = 3L, seed = seed)
  list(data = v$datasets, rng = .Random.seed)
}
a <- run(9L)
for (i in seq_len(15L)) {
  b <- run(9L)
  vals <- unlist(Map(function(q,r) q$x - r$x, a$data, b$data))
  cat('repeat', i, 'identical=', identical(a$data, b$data), 'rng_identical=', identical(a$rng,b$rng), 'max_abs=', max(abs(vals)), 'n_diff=', sum(vals != 0), '\n')
}
c <- run(10L)
cat('different_seed_max_abs=', max(abs(unlist(Map(function(q,r) q$x-r$x,a$data,c$data)))), '\n')
cat('different_seed_rng_identical=', identical(a$rng,c$rng), '\n')
