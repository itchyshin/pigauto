# Deliberately perturb an in-memory function, then require an existing assertion to fail.
# No package source file is changed.
devtools::load_all(quiet = TRUE)
original <- get(".plot_gates", envir = asNamespace("pigauto"))
mutated <- original
code <- paste(deparse(body(original)), collapse = "\n")
stopifnot(grepl("mean(gates[lc])", code, fixed = TRUE))
body(mutated) <- parse(text = sub("mean(gates[lc])", "sum(gates[lc])", code, fixed = TRUE))[[1]]
assignInNamespace(".plot_gates", mutated, ns = "pigauto")
r <- testthat::test_file("tests/testthat/test-plot-contract.R", reporter = "silent",
                         stop_on_failure = FALSE)
counts <- as.data.frame(r)
stopifnot(sum(counts$failed) > 0L, any(grepl("gate plots average", counts$test)))
assignInNamespace(".plot_gates", original, ns = "pigauto")
cat("MUTATION_CONTROL_PASS gate mean-to-sum was detected\n")
