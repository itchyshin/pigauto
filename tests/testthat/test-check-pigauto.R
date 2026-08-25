make_p2_tree <- function() {
  ape::read.tree(text = "((a:1,b:1):1,c:2);")
}

test_that("preprocess retains original row ids across data-only filtering", {
  tree <- make_p2_tree()
  traits <- data.frame(x = c(NA_real_, 1, NA_real_, NA_real_),
                       row.names = c("ghost_before", "a", "ghost_middle", "c"))

  pd <- NULL
  expect_warning(pd <- preprocess_traits(traits, tree), "not found in tree")
  expect_equal(pd$input_row_order, c(2L, NA_integer_, 4L))

  imputed <- data.frame(x = c(11, 22, 33), row.names = tree$tip.label)
  completed <- build_completed(traits, imputed,
                               input_row_order = pd$input_row_order)
  expect_identical(completed$completed$x, c(NA_real_, 1, NA_real_, 33))
  expect_identical(unname(completed$imputed_mask[, "x"]),
                   c(FALSE, FALSE, FALSE, TRUE))
})

test_that("preprocess maps species_col rows around data-only rows back to input", {
  tree <- make_p2_tree()
  traits <- data.frame(species = c("ghost", "c", "a", "ghost", "b"),
                       x = c(NA_real_, NA_real_, 1, NA_real_, NA_real_),
                       stringsAsFactors = FALSE)

  pd <- NULL
  expect_warning(pd <- preprocess_traits(traits, tree, species_col = "species"),
                 "not found in tree")
  expect_equal(pd$input_row_order, c(3L, 5L, 2L))
  imputed <- data.frame(x = c(10, 20, 30))
  completed <- build_completed(traits, imputed, species_col = "species",
                               input_row_order = pd$input_row_order)
  expect_identical(completed$completed$x, c(NA_real_, 30, 1, NA_real_, 20))
  expect_identical(unname(completed$imputed_mask[, "x"]),
                   c(FALSE, TRUE, FALSE, FALSE, TRUE))
})

test_that("check_pigauto returns a structured no-fit preflight", {
  tree <- make_p2_tree()
  traits <- data.frame(x = c(1, NA_real_, 3), row.names = tree$tip.label)

  check <- check_pigauto(traits, tree)
  expect_s3_class(check, "pigauto_check")
  expect_identical(check$schema_version, 1L)
  expect_true(check$status %in% c("ready", "ready_with_warnings"))
  expect_true(check$ready)
  expect_identical(check$fingerprint$algorithm, "md5-xdr-v1")
  expect_match(check$fingerprint$value, "^[0-9a-f]{32}$")
  expect_named(check, c("schema_version", "status", "ready", "fingerprint",
                        "species", "tree", "traits", "runtime", "size", "messages"))
})

test_that("check_pigauto identifies a fully observed matched input as not needed", {
  tree <- make_p2_tree()
  traits <- data.frame(x = c(1, 2, 3), row.names = tree$tip.label)
  check <- check_pigauto(traits, tree)
  expect_identical(check$status, "not_needed")
  expect_false(check$ready)
  expect_true(any(check$messages$code == "fully_observed"))
})

test_that("fingerprints are stable, input-sensitive, and survive RDS", {
  tree <- make_p2_tree()
  dat <- data.frame(x = c(1, NA, 3), row.names = tree$tip.label)
  a <- check_pigauto(dat, tree)
  b <- check_pigauto(dat, tree)
  changed <- dat; changed$x[[1]] <- 9
  c <- check_pigauto(changed, tree)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(a, path)
  expect_identical(a$fingerprint$value, b$fingerprint$value)
  expect_false(identical(a$fingerprint$value, c$fingerprint$value))
  expect_identical(readRDS(path)$fingerprint$value, a$fingerprint$value)
})

test_that("tree and species diagnostics retain explicit consequences", {
  tree <- make_p2_tree()
  dat <- data.frame(species = c("a", "a", "", "ghost"),
                    x = c(1, 2, 3, NA_real_))
  chk <- check_pigauto(dat, tree, species_col = "species")
  expect_identical(chk$status, "error")
  expect_true(any(chk$messages$code == "missing_species_id"))
  expect_true(all(c("matched", "data_only", "tree_only", "repeated") %in% names(chk$species)))
  expect_match(chk$species$data_only$consequence, "completed")

  zero <- tree; zero$edge.length[] <- 0
  branch <- check_pigauto(data.frame(x = c(1, NA, 3), row.names = tree$tip.label), zero)
  expect_true(any(branch$messages$code == "branch_length_all_zero"))
  some_zero <- tree; some_zero$edge.length[[1]] <- 0
  branch2 <- check_pigauto(data.frame(x = c(1, NA, 3), row.names = tree$tip.label), some_zero)
  expect_true(any(branch2$messages$code == "branch_length_zero"))
})

test_that("trait, runtime, size, and no-graph diagnostics are structured", {
  tree <- make_p2_tree()
  dat <- data.frame(count = as.integer(c(1, NA, 3)),
                    ordered = ordered(c("low", NA, "high"), levels = c("low", "high")),
                    row.names = tree$tip.label)
  chk <- check_pigauto(dat, tree, trait_types = c(count = "continuous"))
  expect_true(all(c("trait", "input_columns", "type", "n_missing", "constant") %in% names(chk$traits)))
  expect_identical(chk$traits$type[chk$traits$trait == "count"], "continuous")
  expect_identical(chk$traits$correction[chk$traits$trait == "count"],
                   "")
  expect_true(is.list(chk$runtime))
  expect_true(chk$size$dense_lower_bound_bytes > 0)
  expect_true(chk$size$memory_class %in% c("low", "moderate", "high", "extreme"))

  testthat::local_mocked_bindings(
    build_phylo_graph = function(...) stop("graph must not run"),
    .package = "pigauto"
  )
  expect_no_error(check_pigauto(dat, tree, trait_types = c(count = "continuous")))
})

test_that("automatic integer traits carry a copy-ready correction", {
  tree <- make_p2_tree()
  dat <- data.frame(measure = as.integer(c(1, NA, 3)), row.names = tree$tip.label)
  chk <- check_pigauto(dat, tree)
  expect_identical(chk$traits$correction[[1L]],
                   "trait_types = c(measure = \"continuous\")")
  expect_no_error(parse(text = chk$traits$correction[[1L]]))
  expect_true(any(chk$messages$code == "integer_auto_ambiguous"))
  expect_true(any(chk$messages$action == "trait_types = c(measure = \"continuous\")"))
  weird <- data.frame(`body mass` = as.integer(c(1, NA, 3)), row.names = tree$tip.label, check.names = FALSE)
  expect_identical(check_pigauto(weird, tree)$traits$correction[[1L]],
                   "trait_types = stats::setNames(\"continuous\", \"body mass\")")
  strange_name <- "body`\\mass"
  strange <- data.frame(value = as.integer(c(1, NA, 3)), row.names = tree$tip.label,
                        check.names = FALSE)
  names(strange) <- strange_name
  strange_action <- check_pigauto(strange, tree)$traits$correction[[1L]]
  expect_match(strange_action, "stats::setNames")
  expect_no_error(parse(text = strange_action))
})

test_that("bounded numeric and zero-bearing integer hints are conditional", {
  tree <- make_p2_tree()
  prop <- data.frame(p = c(0, NA_real_, 1), row.names = tree$tip.label)
  cp <- check_pigauto(prop, tree)
  expect_true(any(cp$messages$code == "bounded_numeric_may_be_proportion"))
  expect_true(any(cp$messages$action == "trait_types = c(p = \"proportion\")"))
  expect_no_error(parse(text = cp$messages$action[cp$messages$code == "bounded_numeric_may_be_proportion"]))
  zi <- data.frame(z = as.integer(c(0, NA, 2)), row.names = tree$tip.label)
  cz <- check_pigauto(zi, tree)
  expect_true(any(cz$messages$code == "zero_integer_may_be_zi_count"))
  expect_true(any(cz$messages$action == "trait_types = c(z = \"zi_count\")"))
  expect_no_error(parse(text = cz$messages$action[cz$messages$code == "zero_integer_may_be_zi_count"]))

  named_prop <- data.frame(`body prop` = c(.1, NA_real_, .9), row.names = tree$tip.label, check.names = FALSE)
  declared <- check_pigauto(named_prop, tree, trait_types = c(`body prop` = "proportion"))
  expect_identical(declared$traits$correction[[1L]], "trait_types = stats::setNames(\"proportion\", \"body prop\")")
  expect_no_error(parse(text = declared$traits$correction[[1L]]))

  comp <- data.frame(`red blue` = c(.2, NA, .3), green = c(.8, NA, .7), row.names = tree$tip.label, check.names = FALSE)
  grouped <- check_pigauto(comp, tree, multi_proportion_groups = list(`body comp` = c("red blue", "green")))
  expect_identical(grouped$traits$correction[[1L]], "multi_proportion_groups = stats::setNames(list(c(\"red blue\", \"green\")), \"body comp\")")
  expect_no_error(parse(text = grouped$traits$correction[[1L]]))

  odd_comp <- data.frame(one = c(.2, NA, .3), two = c(.8, NA, .7),
                         row.names = tree$tip.label, check.names = FALSE)
  names(odd_comp) <- c("one`\\", "two\\`")
  odd_group <- "body`\\ group"
  odd <- check_pigauto(odd_comp, tree,
                       multi_proportion_groups = stats::setNames(list(names(odd_comp)), odd_group))
  odd_action <- odd$traits$correction[[1L]]
  expect_match(odd_action, "stats::setNames")
  expect_no_error(parse(text = odd_action))
})

test_that("ambiguity hints use matched target rows and actions parse", {
  tree <- make_p2_tree()
  prop <- data.frame(p = c(.1, NA_real_, .9, 5),
                     row.names = c("a", "b", "c", "ghost"))
  cp <- check_pigauto(prop, tree)
  prop_action <- cp$messages$action[cp$messages$code == "bounded_numeric_may_be_proportion"]
  expect_length(prop_action, 1L)
  expect_no_error(parse(text = prop_action))

  zi <- data.frame(z = as.integer(c(1, NA, 2, 0)),
                   row.names = c("a", "b", "c", "ghost"))
  cz <- check_pigauto(zi, tree)
  expect_false(any(cz$messages$code == "zero_integer_may_be_zi_count"))
})

test_that("branch status is independent of species errors and transforms are semantic", {
  tree <- make_p2_tree()
  dat <- data.frame(x = c(1, NA, 3), row.names = c("a", "ghost", "c"))
  chk <- check_pigauto(dat, tree)
  expect_identical(chk$tree$branch$status, "ok")
  count <- data.frame(x = as.integer(c(1, NA, 3)), row.names = tree$tip.label)
  expect_identical(check_pigauto(count, tree)$traits$transform[[1L]], "log1p")
})

test_that("all-missing, constant, composition, and runtime failures are explicit", {
  tree <- make_p2_tree()
  all_missing <- data.frame(x = c(NA_real_, NA_real_, NA_real_), row.names = tree$tip.label)
  expect_identical(check_pigauto(all_missing, tree)$status, "error")
  constant <- data.frame(x = c(1, NA_real_, 1), row.names = tree$tip.label)
  expect_true(any(check_pigauto(constant, tree)$messages$code == "trait_constant"))
  comp <- data.frame(a = c(.2, NA, .3), b = c(.8, .8, .7), row.names = tree$tip.label)
  cc <- check_pigauto(comp, tree, multi_proportion_groups = list(comp = c("a", "b")))
  expect_identical(cc$traits$n_missing[[1]], 1L)

  testthat::local_mocked_bindings(
    .check_pigauto_runtime = function() list(torch_package = TRUE, torch_is_installed = FALSE,
                                               device = "unavailable", unavailable = TRUE,
                                               detail = "mock runtime unavailable"),
    .package = "pigauto"
  )
  runtime <- check_pigauto(data.frame(x = c(1, NA, 3), row.names = tree$tip.label), tree)
  expect_identical(runtime$status, "error")
  expect_true(any(runtime$messages$code == "torch_unavailable"))
})

test_that("runtime probes the selected accelerator and records failures", {
  expect_identical(.check_pigauto_select_device(cuda_available = FALSE, mps_available = TRUE), "mps")
  expect_identical(.check_pigauto_select_device(cuda_available = TRUE, mps_available = TRUE), "cuda")
  selected <- character()
  mps_probe <- .check_pigauto_probe_device("mps", function(x, device) {
    selected <<- device
    list(x = x, device = device)
  })
  expect_true(mps_probe$ok)
  expect_identical(selected, "mps")
  cuda_selected <- character()
  cuda_ok <- .check_pigauto_probe_device("cuda", function(x, device) {
    cuda_selected <<- device
    list(x = x, device = device)
  })
  expect_true(cuda_ok$ok)
  expect_identical(cuda_selected, "cuda")
  cuda_probe <- .check_pigauto_probe_device("cuda", function(...) stop("broken CUDA"))
  expect_false(cuda_probe$ok)
  expect_match(cuda_probe$detail, "broken CUDA")

  testthat::local_mocked_bindings(
    .check_pigauto_cuda_available = function() TRUE,
    .check_pigauto_mps_available = function() FALSE,
    .check_pigauto_probe_device = function(device, ...) {
      expect_identical(device, "cuda")
      list(ok = FALSE, detail = "cuda scalar tensor probe failed: mocked failure")
    },
    .package = "pigauto"
  )
  broken <- check_pigauto(data.frame(x = c(1, NA_real_, 3), row.names = make_p2_tree()$tip.label), make_p2_tree())
  expect_true(broken$runtime$unavailable)
  expect_identical(broken$runtime$device, "cuda")
  expect_identical(broken$status, "error")
})

test_that("check captures preprocessing messages and impute forwards warnings once", {
  tree <- make_p2_tree()
  tree_only <- data.frame(x = c(1, NA_real_), row.names = c("a", "c"))
  check <- NULL
  expect_silent(check <- check_pigauto(tree_only, tree))
  expect_true(any(check$messages$code == "preprocess_message"))

  dummy_check <- structure(list(
    status = "ready_with_warnings",
    messages = data.frame(code = "branch_length_zero", field = "tree",
                          text = "forward this warning once", action = "inspect",
                          severity = "warning", stringsAsFactors = FALSE)
  ), class = c("pigauto_check", "list"))
  testthat::local_mocked_bindings(
    .check_pigauto_internal = function(...) list(
      check = dummy_check,
      prepared = list(X_scaled = matrix(NA_real_, 1L, 1L), trait_map = list())
    ),
    make_missing_splits = function(...) stop("after warning"),
    .package = "pigauto"
  )
  warnings <- character()
  expect_error(withCallingHandlers(
    impute(data.frame(x = NA_real_, row.names = "a"), tree, missing_frac = 1),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ), "after warning")
  expect_identical(warnings, "forward this warning once")
})

test_that("a broken selected MPS runtime is not ready", {
  tree <- make_p2_tree()
  testthat::local_mocked_bindings(
    .check_pigauto_cuda_available = function() FALSE,
    .check_pigauto_mps_available = function() TRUE,
    .check_pigauto_probe_device = function(device, ...) {
      expect_identical(device, "mps")
      list(ok = FALSE, detail = "mps scalar tensor probe failed: mocked failure")
    },
    .package = "pigauto"
  )
  broken <- check_pigauto(data.frame(x = c(1, NA_real_, 3), row.names = tree$tip.label), tree)
  expect_identical(broken$runtime$device, "mps")
  expect_true(broken$runtime$unavailable)
  expect_identical(broken$status, "error")
})

test_that("zi-count and composition trait descriptions state their semantics", {
  tree <- make_p2_tree()
  zi <- data.frame(z = as.integer(c(0, NA, 2)), row.names = tree$tip.label)
  zcheck <- check_pigauto(zi, tree, trait_types = c(z = "zi_count"))
  expect_identical(zcheck$traits$transform[[1L]], "binary gate + log1p magnitude")

  comp <- data.frame(a = c(.2, NA, .2), b = c(.8, NA, .8), row.names = tree$tip.label)
  ccheck <- check_pigauto(comp, tree, multi_proportion_groups = list(colour = c("a", "b")))
  expect_true(ccheck$traits$constant[[1L]])
})

test_that("impute retains only the small preflight check on a result", {
  skip_if_no_libtorch()
  set.seed(20260825)
  tree <- ape::rtree(12L)
  dat <- data.frame(x = seq_len(12L) + 1, row.names = tree$tip.label)
  dat$x[[4L]] <- NA_real_
  result <- NULL
  expect_warning(result <- impute(dat, tree, epochs = 2L, verbose = FALSE, seed = 1L),
                 "Small validation set")
  expect_s3_class(result$check, "pigauto_check")
  expect_false(any(c("prepared", "raw_traits", "matrices") %in% names(result$check)))
  expect_false(inherits(result$check$tree, "phylo"))
  expect_identical(completed_data(result), result$completed)
})
