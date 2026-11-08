# ordinal_levels on impute() / multi_impute() / multi_impute_trees() /
# check_pigauto() / preprocess_traits(): a declared level order must behave
# identically to a user-constructed ordered factor on every entry point, and
# invalid declarations must fail closed with the read_traits() messages.

ordinal_fixture <- function(n = 24L, seed = 2026L) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  df_chr <- data.frame(
    cont = stats::rnorm(n, 10, 2),
    status = rep(c("low", "mid", "high"), length.out = n),
    row.names = sp,
    stringsAsFactors = FALSE
  )
  df_chr$cont[c(2L, 9L)] <- NA_real_
  df_chr$status[c(4L, n - 2L)] <- NA_character_
  df_ord <- df_chr
  df_ord$status <- ordered(df_ord$status, levels = c("low", "mid", "high"))
  list(tree = tree, chr = df_chr, ord = df_ord,
       levels = c("low", "mid", "high"))
}

test_that("preprocess_traits(ordinal_levels) matches the ordered-factor route", {
  fx <- ordinal_fixture()
  pd_decl <- preprocess_traits(fx$chr, fx$tree,
                               ordinal_levels = list(status = fx$levels))
  pd_manual <- preprocess_traits(fx$ord, fx$tree)

  expect_identical(pd_decl$trait_map$status$type, "ordinal")
  # The declared order, not R's alphabetical sort ("high", "low", "mid").
  expect_identical(pd_decl$trait_map$status$levels, fx$levels)
  expect_identical(pd_decl$trait_map, pd_manual$trait_map)
  expect_identical(pd_decl$X_scaled, pd_manual$X_scaled)
})

test_that("impute(ordinal_levels) equals manual ordered-factor imputation", {
  fx <- ordinal_fixture()
  res_decl <- suppressWarnings(impute(
    fx$chr, fx$tree, ordinal_levels = list(status = fx$levels),
    verbose = FALSE, seed = 11L))
  res_manual <- suppressWarnings(impute(
    fx$ord, fx$tree, verbose = FALSE, seed = 11L))

  expect_s3_class(res_decl$completed$status, "ordered")
  expect_identical(levels(res_decl$completed$status), fx$levels)
  expect_identical(res_decl$completed$status, res_manual$completed$status)
  expect_equal(res_decl$completed$cont, res_manual$completed$cont)
  expect_identical(res_decl$imputed_mask, res_manual$imputed_mask)
  expect_false(anyNA(res_decl$completed$status))
  # Observed cells untouched.
  obs <- !is.na(fx$chr$status)
  expect_identical(as.character(res_decl$completed$status[obs]),
                   fx$chr$status[obs])
})

test_that("ordinal_levels combines with an explicit trait_types declaration", {
  fx <- ordinal_fixture()
  pd <- preprocess_traits(fx$chr, fx$tree,
                          trait_types = c(status = "ordinal"),
                          ordinal_levels = list(status = fx$levels))
  expect_identical(pd$trait_map$status$type, "ordinal")
  expect_identical(pd$trait_map$status$levels, fx$levels)
})

test_that("ordinal_levels declarations fail closed on every entry point", {
  fx <- ordinal_fixture()
  expect_error(
    impute(fx$chr, fx$tree, ordinal_levels = list(status = c("low", "mid")),
           verbose = FALSE),
    "not included in ordinal_levels")
  expect_error(
    preprocess_traits(fx$chr, fx$tree,
                      ordinal_levels = list(nope = fx$levels)),
    "Ordinal columns not found")
  expect_error(
    multi_impute(fx$chr, fx$tree, m = 2L,
                 ordinal_levels = list(c("low", "mid", "high")),
                 verbose = FALSE),
    "uniquely named list")
  expect_error(
    preprocess_traits(fx$chr, fx$tree,
                      ordinal_levels = list(status = c("low", "low", "high"))),
    "unique, nonmissing strings")
  expect_error(
    multi_impute_trees(fx$chr, list(fx$tree, fx$tree), m_per_tree = 1L,
                       ordinal_levels = list(status = c("low", "mid")),
                       verbose = FALSE),
    "not included in ordinal_levels")
})

test_that("declaring ordinal without an order points at ordinal_levels", {
  fx <- ordinal_fixture()
  # Character column: type label alone cannot supply the order.
  expect_error(
    preprocess_traits(fx$chr, fx$tree, trait_types = c(status = "ordinal")),
    "ordinal_levels", fixed = TRUE)
  # Unordered factor: same guidance.
  f <- fx$chr
  f$status <- factor(f$status)
  expect_error(
    preprocess_traits(f, fx$tree, trait_types = c(status = "ordinal")),
    "ordinal_levels", fixed = TRUE)
})

test_that("safety advisories fire identically on both ordinal routes", {
  # The Phase B3 imbalance advisory keys on is.ordered(); the declared route
  # must trigger it exactly like a user-constructed ordered factor would.
  set.seed(99L)
  n <- 30L
  tree <- ape::rtree(n)
  status <- c(rep("low", 24L), rep("mid", 3L), rep("high", 3L))
  df_chr <- data.frame(status = status, row.names = tree$tip.label,
                       stringsAsFactors = FALSE)
  df_chr$status[c(3L, 17L)] <- NA_character_
  df_ord <- df_chr
  df_ord$status <- ordered(df_ord$status, levels = c("low", "mid", "high"))

  grab_warnings <- function(expr) {
    w <- character()
    withCallingHandlers(expr, warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
    w
  }
  w_decl <- grab_warnings(impute(
    df_chr, tree, ordinal_levels = list(status = c("low", "mid", "high")),
    verbose = FALSE, seed = 3L))
  w_manual <- grab_warnings(impute(df_ord, tree, verbose = FALSE, seed = 3L))

  expect_true(any(grepl("collapse to the majority class", w_decl)))
  expect_true(any(grepl("collapse to the majority class", w_manual)))
})

test_that("the missing-order error suggests valid code for non-syntactic names", {
  fx <- ordinal_fixture()
  df <- fx$chr
  names(df)[names(df) == "status"] <- "body size"
  expect_error(
    preprocess_traits(df, fx$tree,
                      trait_types = stats::setNames("ordinal", "body size")),
    "list(`body size` = c(...))", fixed = TRUE)
})

test_that("a two-level declaration detects as binary with the order preserved", {
  fx <- ordinal_fixture()
  df <- fx$chr
  df$status <- ifelse(is.na(df$status), NA_character_,
                      ifelse(df$status == "low", "small", "large"))
  pd <- preprocess_traits(df, fx$tree,
                          ordinal_levels = list(status = c("small", "large")))
  expect_identical(pd$trait_map$status$type, "binary")
  expect_identical(pd$trait_map$status$levels, c("small", "large"))
})

test_that("check_pigauto reports declared ordinals and bad declarations", {
  fx <- ordinal_fixture()
  chk <- suppressWarnings(check_pigauto(
    fx$chr, fx$tree, ordinal_levels = list(status = fx$levels)))
  expect_true(chk$ready)
  expect_identical(chk$traits$type[chk$traits$trait == "status"], "ordinal")

  # An invalid declaration is reported, not thrown: check_pigauto()'s
  # report-not-throw contract holds for this argument too.
  bad <- suppressWarnings(check_pigauto(
    fx$chr, fx$tree, ordinal_levels = list(status = c("low", "mid"))))
  expect_identical(bad$status, "error")
  expect_true(any(bad$messages$code == "ordinal_levels_invalid"))
  expect_false(bad$ready)
})

test_that("multi_impute(ordinal_levels) returns ordered factors in datasets", {
  fx <- ordinal_fixture()
  mi <- suppressWarnings(suppressMessages(multi_impute(
    fx$chr, fx$tree, m = 2L, draws_method = "conformal",
    ordinal_levels = list(status = fx$levels),
    verbose = FALSE, seed = 7L)))
  datasets <- get_imputed_datasets(mi)
  expect_length(datasets, 2L)
  for (d in datasets) {
    expect_s3_class(d$status, "ordered")
    expect_identical(levels(d$status), fx$levels)
    expect_false(anyNA(d$status))
  }
})

test_that("multi_impute_trees(ordinal_levels) coerces before each tree fit", {
  fx <- ordinal_fixture(n = 12L, seed = 5L)
  tree2 <- ape::rtree(12L)
  tree2$tip.label <- fx$tree$tip.label
  out <- suppressWarnings(multi_impute_trees(
    fx$chr, list(fx$tree, tree2), m_per_tree = 1L,
    ordinal_levels = list(status = fx$levels),
    verbose = FALSE, seed = 5L))
  datasets <- get_imputed_datasets(out)
  expect_s3_class(datasets[[1L]]$status, "ordered")
  expect_identical(levels(datasets[[1L]]$status), fx$levels)
})

test_that("new argument preserves the existing positional API", {
  # Calls written before ordinal_levels existed must continue to bind the
  # composition declaration to multi_proportion_groups on every entry point.
  legacy_calls <- list(
    list(fun = impute, call = quote(impute(d, tree, NULL, NULL, groups))),
    list(fun = preprocess_traits, call = quote(preprocess_traits(d, tree, NULL, NULL, groups))),
    list(fun = check_pigauto, call = quote(check_pigauto(d, tree, NULL, NULL, groups))),
    list(fun = multi_impute, call = quote(multi_impute(d, tree, 2L, "conformal", NULL, NULL, groups))),
    list(fun = multi_impute_trees, call = quote(multi_impute_trees(d, trees, 1L, NULL, NULL, groups)))
  )
  for (case in legacy_calls) {
    matched <- match.call(definition = case$fun, call = case$call)
    expect_identical(matched$multi_proportion_groups, quote(groups))
    expect_null(matched$ordinal_levels)
  }

  set.seed(42)
  tree <- ape::rtree(6)
  d <- data.frame(a = c(.2,.3,.4,.5,.6,NA),
                  b = c(.8,.7,.6,.5,.4,NA), row.names = tree$tip.label)
  groups <- list(comp = c("a", "b"))
  named <- suppressWarnings(preprocess_traits(d, tree, multi_proportion_groups = groups))
  positional <- suppressWarnings(preprocess_traits(d, tree, NULL, NULL, groups))
  expect_identical(positional$trait_map, named$trait_map)
  expect_identical(positional$X_scaled, named$X_scaled)
})

test_that("suggested ordinal declaration parses for difficult names", {
  fx <- ordinal_fixture()
  for (tricky in c("body`size", "body\\size")) {
    df <- fx$chr
    names(df)[names(df) == "status"] <- tricky
    msg <- tryCatch(preprocess_traits(df, fx$tree,
      trait_types = stats::setNames("ordinal", tricky)), error = conditionMessage)
    code <- sub(".*ordinal_levels = (.*) listing the levels.*", "\\1", msg)
    code <- sub("c(...)", "c('low', 'mid', 'high')", code, fixed = TRUE)
    declaration <- eval(parse(text = code))
    expect_identical(names(declaration), tricky)
    expect_identical(declaration[[1]], c("low", "mid", "high"))
  }
})
