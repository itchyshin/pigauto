# Exercise the production dependency-parser expressions without loading packages.
expressions <- parse("script/cran-0.11-platform/check.R")
parse_dependencies <- function(description) {
  env <- new.env()
  env$desc <- description
  for (expression in expressions) {
    if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        as.character(expression[[2]]) %in% c("fields", "required")) {
      eval(expression, env)
    }
  }
  env$required
}
desc <- read.dcf("frozen-deps/DESCRIPTION")
actual <- parse_dependencies(desc)
stopifnot(all(nzchar(actual)), !"R" %in% actual)
db <- desc
for (field in setdiff(c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances"), colnames(db))) {
  db <- cbind(db, matrix(NA_character_, 1, 1, dimnames = list(NULL, field)))
}
rownames(db) <- desc[1, "Package"]
expected <- tools::package_dependencies(rownames(db), db = db,
  which = c("Depends", "Imports", "LinkingTo", "Suggests"), recursive = FALSE)[[1]]
stopifnot(setequal(actual, setdiff(expected, "R")))
fixture <- matrix(c("R (>= 4.1)", "ape,", "testthat,"), 1,
                  dimnames = list(NULL, c("Depends", "Imports", "Suggests")))
stopifnot(setequal(parse_dependencies(fixture), c("ape", "testthat")))
cat("DEPENDENCY_PARSER_CONTROLS_OK packages=", length(actual), "\n", sep = "")
