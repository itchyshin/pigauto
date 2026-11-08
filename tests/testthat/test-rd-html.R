test_that("help tables render as valid HTML paragraphs", {
  tidy <- Sys.which("tidy")
  if (!nzchar(tidy)) skip("HTML Tidy is not installed")
  root <- testthat::test_path("..", "..")
  pages <- file.path(root, "man", c("predict.pigauto_fit.Rd", "preprocess_traits.Rd"))
  if (!all(file.exists(pages))) skip("source Rd files are unavailable")
  check_html <- function(path) {
    suppressWarnings(system2(tidy, c("-errors", "-quiet", shQuote(path)),
                             stdout = TRUE, stderr = TRUE))
  }
  invalid <- tempfile(fileext = ".html")
  writeLines('<!DOCTYPE html><html><head><title>test</title></head><body><p>x<table><tr><td>y</td></tr></table></p></body></html>', invalid)
  expect_true(any(grepl("inserting implicit <p>|trimming empty <p>",
                        check_html(invalid))))
  for (page in pages) {
    html <- tempfile(fileext = ".html")
    tools::Rd2HTML(page, out = html)
    expect_false(any(grepl("inserting implicit <p>|trimming empty <p>",
                           check_html(html))), info = basename(page))
  }
})
