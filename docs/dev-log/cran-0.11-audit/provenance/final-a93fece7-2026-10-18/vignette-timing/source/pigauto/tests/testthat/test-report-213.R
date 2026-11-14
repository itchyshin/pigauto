# Issue #213: pigauto_report() validates data/splits/title/output_path/open
# before writing, and HTML-escapes title so markup cannot inject into the file.

tiny_report_fit <- function() {
  set.seed(213L)
  tree <- ape::rtree(16)
  df <- data.frame(
    mass  = abs(as.numeric(ape::rTraitCont(tree))) + 0.5,
    range = abs(as.numeric(ape::rTraitCont(tree))) + 0.5,
    row.names = tree$tip.label
  )
  df$mass[seq_len(2L)] <- NA_real_
  pd <- preprocess_traits(df, tree)
  spl <- make_missing_splits(
    pd$X_scaled, missing_frac = 0.25, seed = 213L, trait_map = pd$trait_map
  )
  fit <- fit_pigauto(
    pd, tree, splits = spl, gnn = FALSE, verbose = FALSE, seed = 213L
  )
  list(fit = fit, pd = pd, spl = spl, df = df)
}

other_report_data <- function() {
  set.seed(214L)
  tree <- ape::rtree(20)
  df <- data.frame(
    mass  = abs(as.numeric(ape::rTraitCont(tree))) + 0.5,
    range = abs(as.numeric(ape::rTraitCont(tree))) + 0.5,
    row.names = tree$tip.label
  )
  pd <- preprocess_traits(df, tree)
  spl <- make_missing_splits(
    pd$X_scaled, missing_frac = 0.25, seed = 214L, trait_map = pd$trait_map
  )
  list(pd = pd, spl = spl)
}

test_that("pigauto_report rejects bad data, splits, title, output_path, and open", {
  obj <- tiny_report_fit()
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)
  rpt <- function(...) {
    pigauto_report(obj$fit, open = FALSE, output_path = tmp, ...)
  }

  expect_error(rpt(data = obj$df, splits = obj$spl), "pigauto_data")
  expect_error(rpt(data = list(a = 1)), "pigauto_data")
  expect_error(rpt(data = obj$pd, splits = list()), "splits")
  expect_error(rpt(data = obj$pd, splits = matrix(1, 2, 2)), "splits")
  expect_error(rpt(data = obj$pd, splits = "x"), "splits")

  other <- other_report_data()
  expect_error(
    rpt(data = other$pd, splits = obj$spl),
    "does not match fit"
  )
  expect_error(
    rpt(data = obj$pd, splits = other$spl),
    "do not match fit"
  )

  expect_error(rpt(title = c("a", "b")), "title")
  expect_error(rpt(title = NULL), "title")
  expect_error(
    pigauto_report(obj$fit, open = FALSE, output_path = NULL),
    "output_path"
  )
  expect_error(
    pigauto_report(obj$fit, open = FALSE, output_path = c("a.html", "b.html")),
    "output_path"
  )
  expect_error(
    pigauto_report(obj$fit, open = FALSE, output_path = ""),
    "output_path"
  )
  expect_error(
    pigauto_report(obj$fit, open = NA, output_path = tmp),
    "open"
  )
})

test_that("pigauto_report HTML-escapes a title that contains markup", {
  obj <- tiny_report_fit()
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)
  nasty <- "My <b>mammals</b> & <script>x</script>"
  pigauto_report(obj$fit, open = FALSE, output_path = tmp, title = nasty)
  html <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  expect_false(grepl("<script>x</script>", html, fixed = TRUE))
  expect_true(grepl("&lt;script&gt;x&lt;/script&gt;", html, fixed = TRUE))
  expect_true(grepl("&lt;b&gt;mammals&lt;/b&gt;", html, fixed = TRUE))
  expect_true(grepl("My &lt;b&gt;mammals&lt;/b&gt; &amp; &lt;script&gt;x&lt;/script&gt;",
                    html, fixed = TRUE))
})
