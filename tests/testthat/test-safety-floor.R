# tests/testthat/test-safety-floor.R
# Smoke canary for safety-floor (three-way mean-gate) calibration.
# See specs/2026-04-23-safety-floor-mean-gate-design.md.

test_that("simplex_grid(0.05) returns 231 rows on the simplex", {
  g <- pigauto:::simplex_grid(step = 0.05)
  expect_equal(dim(g), c(231L, 3L))
  expect_equal(colnames(g), c("r_bm", "r_gnn", "r_mean"))
  expect_true(all(g >= 0 & g <= 1))
  expect_equal(rowSums(g), rep(1, 231L), tolerance = 1e-10)
})

test_that("simplex_grid(0.25) returns 15 rows (coarse grid for tests)", {
  g <- pigauto:::simplex_grid(step = 0.25)
  expect_equal(nrow(g), 15L)
  expect_equal(rowSums(g), rep(1, 15L), tolerance = 1e-10)
})

test_that("simplex_grid always contains the three corners", {
  g <- pigauto:::simplex_grid(step = 0.05)
  corners <- matrix(c(1,0,0, 0,1,0, 0,0,1), nrow = 3L, byrow = TRUE,
                    dimnames = list(NULL, c("r_bm", "r_gnn", "r_mean")))
  for (k in seq_len(nrow(corners))) {
    matched <- any(apply(g, 1L, function(row) all(row == corners[k, ])))
    expect_true(matched,
                info = sprintf("corner (%g,%g,%g) missing",
                               corners[k,1], corners[k,2], corners[k,3]))
  }
})
