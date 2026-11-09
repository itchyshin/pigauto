# Probe that module registration adds counters without clearing existing ones.
root <- normalizePath(".")
pkgload::load_all(root, quiet = TRUE, helpers = FALSE, export_all = FALSE)
torch::torch_set_num_threads(1L)
torch::torch_set_num_interop_threads(1L)
ns <- asNamespace("pigauto")
module_bodies <- function() lapply(c("GraphTransformerBlock", "ResidualPhyloDAE"),
  function(name) lapply(attr(get(name, ns), "module")$public_methods[c("initialize", "forward")],
                        function(fun) paste(deparse(body(fun)), collapse="\n")))
before_bodies <- module_bodies()
probe <- tempfile(fileext = ".R")
on.exit(unlink(probe), add = TRUE)
writeLines(c(
  "stopifnot(trimws(format_metric(.5)) == '0.500')",
  "prior <- as.list(get('.counters', asNamespace('covr')))",
  "before <- sum(vapply(prior, function(x) x$value, numeric(1)))",
  sprintf("sys.source(%s, envir = new.env(parent = globalenv()))", encodeString(file.path(root,'tests/testthat/helper-coverage-torch.R'), quote='"')),
  "after <- as.list(get('.counters', asNamespace('covr')))",
  "stopifnot(all(names(prior) %in% names(after)))",
  "stopifnot(identical(before, sum(vapply(after[names(prior)], function(x) x$value, numeric(1)))))",
  "b <- pigauto:::GraphTransformerBlock(hidden_dim=8L,n_heads=2L,dropout=0)",
  "h <- torch::torch_ones(3,8); adj <- torch::torch_eye(3)",
  "stopifnot(all(is.finite(as.numeric(b(h,adj)))))"), probe)
e <- new.env(parent = globalenv())
e$format_metric <- get("format_metric", envir = ns)
cov <- covr::environment_coverage(e, test_files = probe)
tally <- covr::tally_coverage(cov, by = "line")
files <- basename(tally$filename)
stopifnot(any(files == "evaluate.R" & tally$value > 0),
          any(files == "graph_transformer_block.R" & tally$value > 0),
          any(files == "model_residual_dae.R" & tally$value == 0))
stopifnot(identical(before_bodies, module_bodies()))
cat("REGISTRATION_SELFTEST_PASS retained primary counters and added both modules\n")
