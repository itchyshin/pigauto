test_that("missingness plots use traits, not expanded latent columns", {
  x <- structure(list(X_scaled = matrix(c(1, NA, 1, NA, 0, NA), 2, 3),
    trait_map = list(list(name = "mass", type = "continuous", latent_cols = 1),
                     list(name = "diet", type = "categorical", latent_cols = 2:3)),
    species_names = c("a", "b"), obs_to_species = NULL), class = "pigauto_data")
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off())
  m <- plot_missingness(x, "matrix")
  expect_equal(unname(m), matrix(c(FALSE, TRUE, FALSE, TRUE), 2))
  expect_equal(colnames(m), c("mass", "diet"))
  expect_silent(plot_missingness(x, "co_missing"))
  expect_error(plot_missingness(x, "phylo"), "phylo")
})

test_that("baseline report figures show finite metrics and offline SVG", {
  metrics <- list(list(name = "mass & size", rmse_gnn = 0.5),
                  list(name = "diet", accuracy_gnn = 0.75))
  fig <- pigauto:::.report_svg_metrics(metrics, "rmse_gnn", "RMSE")
  expect_match(fig, "<svg", fixed = TRUE)
  expect_match(fig, "mass &amp; size", fixed = TRUE)
  expect_match(fig, "0.500", fixed = TRUE)
  expect_match(pigauto:::.report_svg_metrics(metrics, "conformal_score", "Scores"),
               "No data available", fixed = TRUE)
})

test_that("phylogenetic map orders repeated observations by tree traversal", {
  tree <- ape::read.tree(text = "((a,c),b);")
  x <- structure(list(X_scaled = matrix(c(NA, 1, NA, 1), 4, 1),
    trait_map = list(list(name = "mass", type = "continuous", latent_cols = 1)),
    species_names = c("b", "a", "c"), obs_to_species = c(1L, 2L, 3L, 2L)),
    class = "pigauto_data")
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off())
  ordered <- plot_missingness(x, "phylo", tree)
  expect_equal(as.vector(ordered), c(FALSE, FALSE, TRUE, TRUE))
})

test_that("single-trait matrices plot and graphics parameters are restored", {
  x <- structure(list(X_scaled = matrix(c(1, NA), 2, 1),
    trait_map = list(list(name = "mass", type = "continuous", latent_cols = 1)),
    species_names = c("a", "b"), obs_to_species = NULL), class = "pigauto_data")
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 8, height = 6)
  on.exit(grDevices::dev.off())
  before <- graphics::par("mar")
  expect_silent(plot_missingness(x, "matrix"))
  expect_silent(plot_missingness(x, "co_missing"))
  expect_equal(graphics::par("mar"), before)
})

test_that("bare fit reports the same held-out metrics as its result", {
  set.seed(7)
  tree <- ape::rtree(30)
  traits <- data.frame(mass = abs(as.numeric(ape::rTraitCont(tree))) + 1,
                       wing = abs(as.numeric(ape::rTraitCont(tree))) + 1,
                       row.names = tree$tip.label)
  traits$mass[1:3] <- NA_real_
  result <- suppressWarnings(impute(traits, tree, gnn = FALSE, verbose = FALSE, seed = 7))
  paths <- c(tempfile(fileext = ".html"), tempfile(fileext = ".html"))
  on.exit(unlink(paths), add = TRUE)
  pigauto_report(result, output_path = paths[1], open = FALSE)
  pigauto_report(result$fit, output_path = paths[2], open = FALSE)
  html <- lapply(paths, function(p) paste(readLines(p, warn = FALSE), collapse = "\n"))
  table <- lapply(html, function(h) regmatches(h, regexpr("<tbody>.*?</tbody>", h)))
  expect_identical(table[[1]], table[[2]])
  expect_match(table[[1]], "[0-9]+\\.[0-9]{4}")
  expect_match(html[[2]], "Held-out RMSE (latent scale)", fixed = TRUE)
  rmse_panel <- regmatches(html[[2]], regexpr("<h3>Held-out RMSE \\(latent scale\\)</h3>.*?</div>", html[[2]]))
  expect_match(rmse_panel, "<svg", fixed = TRUE)
})
