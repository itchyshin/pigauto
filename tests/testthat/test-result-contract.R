test_that("completed_data returns the stored completed object exactly", {
  completed <- data.frame(x = c(1, 2), row.names = c("a", "b"))
  result <- structure(list(completed = completed), class = "pigauto_result")
  expect_identical(completed_data(result), result$completed)
})

test_that("result completion preserves input rows, classes, observed cells, and data-only gaps", {
  original <- data.frame(species = c("ghost", "a", "ghost", "c"),
                         value = c(NA_real_, 1, NA_real_, NA_real_),
                         group = factor(c("g", "a", "g", "c")),
                         stringsAsFactors = FALSE)
  imputed <- data.frame(value = c(10, 30), group = factor(c("a", "c")))
  info <- build_completed(original, imputed, species_col = "species",
                          input_row_order = c(2L, 4L))
  result <- structure(list(completed = info$completed, imputed_mask = info$imputed_mask),
                      class = "pigauto_result")
  expect_identical(completed_data(result), result$completed)
  expect_identical(result$completed$species, original$species)
  expect_identical(result$completed$value, c(NA_real_, 1, NA_real_, 30))
  expect_s3_class(result$completed$group, "factor")
  expect_identical(unname(result$imputed_mask[, "value"]), c(FALSE, FALSE, FALSE, TRUE))
})

test_that("character completion keeps character values rather than factor codes", {
  original <- data.frame(x = c("A", NA_character_, "B"), stringsAsFactors = FALSE)
  imputed <- data.frame(x = factor(c("A", "B", "B")))
  out <- build_completed(original, imputed, input_row_order = 1:3)$completed
  expect_identical(out$x, c("A", "B", "B"))
  expect_type(out$x, "character")
})

test_that("result summary exposes the shared output-role contract", {
  completed <- data.frame(x = c(1, 2), row.names = c("a", "b"))
  result <- structure(
    list(completed = completed,
         imputed_mask = matrix(c(FALSE, TRUE), ncol = 1,
                               dimnames = list(c("a", "b"), "x")),
         check = list(status = "ready", ready = TRUE),
         data = list(species_names = c("a", "b"),
                     trait_map = list(x = list(name = "x", type = "continuous"))),
         evaluation = NULL),
    class = "pigauto_result"
  )
  out <- summary(result)
  expect_s3_class(out, "summary_pigauto_result")
  expect_true(is.data.frame(out$output_roles))
  expect_true(any(out$output_roles$role == "completed"))
  expect_true(any(grepl("does not authorize inference", out$output_roles$text,
                        fixed = TRUE)))
})

test_that("summary partitions target cells into observed, filled, and unresolved", {
  completed <- data.frame(x = c(1, 9, NA_real_, NA_real_), species = c("a", "b", "c", "ghost"))
  mask <- matrix(c(FALSE, TRUE, FALSE, FALSE), ncol = 1, dimnames = list(NULL, "x"))
  check <- structure(list(species = list(species_col = "species", data_only = list(names = "ghost"))), class = c("pigauto_check", "list"))
  out <- summary(structure(list(completed = completed, imputed_mask = mask, check = check, evaluation = NULL), class = "pigauto_result"))
  expect_identical(out$counts, list(filled = 1L, observed = 1L, unresolved = 1L))
  expect_output(print(out), "observed cells.*1")
  expect_output(print(out), "inference=not authorized")
})

test_that("summary treats an attempted but missing fill as unresolved", {
  completed <- data.frame(x = c(1, NA_real_), row.names = c("a", "b"))
  mask <- matrix(c(FALSE, TRUE), ncol = 1, dimnames = list(c("a", "b"), "x"))
  out <- summary(structure(
    list(completed = completed, imputed_mask = mask, check = NULL, evaluation = NULL),
    class = "pigauto_result"
  ))
  expect_identical(out$counts, list(filled = 0L, observed = 1L, unresolved = 1L))
  expect_identical(out$per_trait$filled, 0L)
  expect_identical(sum(unlist(out$counts)), 2L)
})

test_that("report contract handles legacy checks and escapes input messages", {
  result <- structure(list(completed = data.frame(x = 1),
                           check = list(status = "ready<script>",
                                        messages = data.frame(code = "<code>", field = "x",
                                                              text = "<unsafe>", action = "use 'x'",
                                                              severity = "warning"))),
                      class = "pigauto_result")
  html <- .report_contract_html(result$check, result)
  expect_match(html, "&lt;unsafe&gt;", fixed = TRUE)
  expect_false(grepl("<unsafe>", html, fixed = TRUE))
  legacy <- .report_contract_html(NULL, NULL)
  expect_match(legacy, "Preflight check unavailable", fixed = TRUE)
})
