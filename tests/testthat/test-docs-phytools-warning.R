# Issues #196 / #197 / #198 / #215: docs match the prediction object;
# phytools-missing warning names an install step.

pkg_src <- function(...) {
  file.path(testthat::test_path("../.."), ...)
}

skip_if_no_pkg_src <- function(...) {
  path <- pkg_src(...)
  if (!file.exists(path)) {
    skip(paste("source file not installed during R CMD check:", path))
  }
  path
}

test_that("getting-started Step 6 uses matrix indexing and evaluate() for coverage", {
  gs <- readLines(skip_if_no_pkg_src("vignettes", "getting-started.Rmd"))
  text <- paste(gs, collapse = "\n")
  expect_false(grepl('conformal_lower\\[\\["Mass"\\]\\]', text))
  expect_true(grepl('conformal_lower\\[,\\s*"Mass"\\]', text))
  expect_false(grepl("pred\\$conformal_coverage", text))
  expect_true(grepl("evaluate\\(\\)", text))
})

test_that("install text does not send a CRAN 0.10.0 user into check_pigauto()", {
  rd <- paste(readLines(skip_if_no_pkg_src("README.md")), collapse = "\n")
  gs <- paste(readLines(skip_if_no_pkg_src("vignettes", "getting-started.Rmd")),
              collapse = "\n")
  expect_true(grepl("0\\.10\\.0", rd))
  expect_true(grepl("0\\.10\\.0", gs))
  expect_true(grepl("does not include `check_pigauto", rd))
  expect_true(grepl("does not have them", gs))
})

test_that("phytools-missing warning names an install step", {
  src <- paste(readLines(skip_if_no_pkg_src("R", "phylo_signal.R")), collapse = "\n")
  expect_true(grepl("phytools", src, ignore.case = TRUE))
  expect_true(grepl("install", src, ignore.case = TRUE))
  expect_true(grepl("install.packages", src, fixed = TRUE))
})

test_that("docs state the split-conformal clade assumption", {
  rd <- paste(readLines(skip_if_no_pkg_src("README.md")), collapse = "\n")
  gs <- paste(readLines(skip_if_no_pkg_src("vignettes", "getting-started.Rmd")),
              collapse = "\n")
  expect_true(grepl("clade", rd, ignore.case = TRUE))
  expect_true(grepl("clade", gs, ignore.case = TRUE))
  expect_true(grepl("resemble", rd, ignore.case = TRUE))
  expect_true(grepl("resemble", gs, ignore.case = TRUE))
})

test_that("Wave 2 coverage figures are not quoted as a new measurement", {
  rd <- paste(readLines(skip_if_no_pkg_src("README.md")), collapse = "\n")
  gs <- paste(readLines(skip_if_no_pkg_src("vignettes", "getting-started.Rmd")),
              collapse = "\n")
  blob <- paste(rd, gs, collapse = "\n")
  isolated <- function(fig) {
    grepl(paste0("(^|[^0-9.])", gsub("\\.", "\\\\.", fig), "([^0-9]|$)"), blob)
  }
  for (fig in c("0.94", "0.88", "0.79", "0.18")) {
    if (isolated(fig)) {
      expect_true(grepl("#215", blob) || grepl("Wave 2", blob) ||
                    grepl("07abad3", blob))
    }
  }
  succeed()
})
