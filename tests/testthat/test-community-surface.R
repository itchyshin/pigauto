test_that("canonical novice journey is present in the public sources", {
  root <- testthat::test_path("..", "..")
  journey <- c(
    'traits <- read_traits("traits.csv")',
    'tree <- read_tree("tree.nwk")',
    "check_pigauto(traits, tree)",
    "result <- impute(traits, tree)",
    "completed <- completed_data(result)",
    "pigauto_report(result)"
  )
  for (path in c("README.md", file.path("vignettes", "getting-started.Rmd"))) {
    text <- paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
    expect_true(all(vapply(journey, grepl, logical(1), x = text, fixed = TRUE)), info = path)
  }
})

test_that("installed novice fixture and script parse without fitting", {
  fixture_dir <- testthat::test_path("..", "..", "inst", "extdata")
  traits_path <- file.path(fixture_dir, "novice_traits.csv")
  tree_path <- file.path(fixture_dir, "novice_tree.nwk")
  script_path <- testthat::test_path("..", "..", "inst", "examples", "novice-workflow.R")
  expect_true(file.exists(traits_path))
  expect_true(file.exists(tree_path))
  expect_true(file.exists(script_path))
  exprs <- parse(script_path)
  expect_length(exprs, 6L)
  expression_text <- vapply(exprs, function(x) paste(deparse(x), collapse = ""), character(1))
  expect_match(expression_text[[1L]], "^traits <- pigauto::read_traits\\(")
  expect_match(expression_text[[2L]], "^tree <- pigauto::read_tree\\(")
  expect_match(expression_text[[3L]], "^pigauto::check_pigauto\\(traits, tree\\)$")
  expect_match(expression_text[[4L]], "^result <- pigauto::impute\\(traits, tree\\)$")
  expect_match(expression_text[[5L]], "^completed <- pigauto::completed_data\\(result\\)$")
  expect_match(expression_text[[6L]], "^pigauto::pigauto_report\\(result\\)$")
  traits <- read_traits(traits_path)
  tree <- read_tree(tree_path)
  expect_s3_class(traits, "data.frame")
  expect_s3_class(tree, "phylo")
  expect_true(anyNA(traits))
  check <- check_pigauto(traits, tree)
  expect_true(check$status %in% c("ready", "ready_with_warnings"))
  expect_no_error(parse(text = 'traits$state <- ordered(traits$state, levels = c("low", "high"))'))
  expect_no_error(parse(text = 'trait_types = c(length_mm = "continuous", cover = "proportion", parasites = "zi_count")'))
  expect_no_error(parse(text = 'multi_proportion_groups = list(diet = c("diet_insect", "diet_fruit", "diet_seed"))'))
})

test_that("recipes, advanced controls, and bounded claims remain visible", {
  root <- testthat::test_path("..", "..")
  readme <- paste(readLines(file.path(root, "README.md"), warn = FALSE), collapse = "\n")
  getting <- paste(readLines(file.path(root, "vignettes", "getting-started.Rmd"), warn = FALSE), collapse = "\n")
  mixed <- paste(readLines(file.path(root, "vignettes", "mixed-types.Rmd"), warn = FALSE), collapse = "\n")
  expect_match(readme, "Ordered ecological states")
  expect_match(readme, "Integer-valued continuous")
  expect_match(readme, "Repeated observations")
  expect_match(readme, "check\\$species")
  expect_match(readme, "## Advanced controls", fixed = TRUE)
  expect_match(readme, 'predict_method = "exact"', fixed = TRUE)
  expect_match(getting, "nominal held-out diagnostic", fixed = TRUE)
  expect_match(mixed, "multi_impute_analysis", fixed = TRUE)
  expect_false(grepl("maximally reduce imputation uncertainty|Why this is novel|certified coverage", readme, ignore.case = TRUE))
})

test_that("current public claim surfaces retain their stated boundaries", {
  root <- testthat::test_path("..", "..")
  read_text <- function(...) paste(readLines(file.path(root, ...), warn = FALSE), collapse = "\n")
  expect_false(grepl("ensures the network only|ensures improvement", read_text("DESCRIPTION"), ignore.case = TRUE))
  expect_false(grepl("maximise imputation precision|maximally reduce|Why this is novel|appears to be new", read_text("R", "active_impute.R"), ignore.case = TRUE))
  expect_false(grepl("0\\.14-1\\.27|lower z-RMSE|~13%", read_text("R", "fit_baseline.R")))
  expect_false(grepl("\\+6\\.6|will outperform|materially outperforms|typically helps|λ > 0\\.7", read_text("vignettes", "common-pitfalls.Rmd"), ignore.case = TRUE))
  expect_false(grepl("confirmed .* beats|z-RMSE 1\\.038|giving \\>= 95", read_text("vignettes", "gnn-architecture.Rmd"), ignore.case = TRUE))
  expect_match(read_text("R", "active_impute.R"), "model-based BM/label-propagation proxy", fixed = TRUE)
  expect_match(read_text("vignettes", "getting-started.Rmd"), "nominal\\s+held-out diagnostics")
})

test_that("tree benchmark and interval help retain current boundaries", {
  root <- testthat::test_path("..", "..")
  read_text <- function(...) paste(readLines(file.path(root, ...), warn = FALSE), collapse = "\n")
  pkgdown <- read_text("_pkgdown.yml")
  tombstone <- read_text("pkgdown", "assets", "dev", "bench_tree_uncertainty.html")
  tree_vignette <- read_text("vignettes", "tree-uncertainty.Rmd")
  predict_source <- read_text("R", "predict_pigauto.R")
  predict_rd <- read_text("man", "predict.pigauto_fit.Rd")
  fit_source <- read_text("R", "fit_pigauto.R")
  fit_rd <- read_text("man", "fit_pigauto.Rd")

  expect_false(grepl("dev/bench_tree_uncertainty\\.html", pkgdown))
  expect_match(tombstone, "Historical tree-sensitivity benchmark", fixed = TRUE)
  expect_match(tombstone, "articles/tree-uncertainty.html", fixed = TRUE)
  expect_match(tombstone, "not MI, Rubin pooling, SE, or FMI", fixed = TRUE)
  expect_false(grepl("Pooled estimates|Rubin's rules", tombstone, fixed = TRUE))
  expect_match(tree_vignette, "descriptive prediction sensitivity only", fixed = TRUE)
  expect_match(tree_vignette, "do not pass them to `pool_mi()`", fixed = TRUE)

  for (text in list(predict_source, predict_rd)) {
    expect_false(grepl("safe for interval arithmetic", text, fixed = TRUE))
    expect_match(text, "model-dependent BM/joint conditional SD", fixed = TRUE)
    expect_match(text, "not total final blended-prediction uncertainty", fixed = TRUE)
  }
  for (text in list(fit_source, fit_rd)) {
    plain <- gsub("[[:space:]]+", " ", text)
    expect_false(grepl("\\b(near[- ]exact|guarantee[sd]?)\\b.{0,100}\\b(coverage|interval)", plain,
      ignore.case = TRUE, perl = TRUE
    ))
    expect_false(grepl("\\b(restore[sd]?|fix(?:es)?)\\b.{0,100}\\bexchangeability\\b", plain,
      ignore.case = TRUE, perl = TRUE
    ))
    expect_false(grepl("use this when.{0,100}\\b(coverage|interval)", plain,
      ignore.case = TRUE, perl = TRUE
    ))
    expect_false(grepl("\\b(coverage|interval).{0,80}\\b(SD|range|ranging|seed|simulation|benchmark)\\b|\\b(mean|empirical).{0,80}\\bcoverage\\b", plain,
      ignore.case = TRUE, perl = TRUE
    ))
    expect_match(plain, "nominal held-out diagnostic", fixed = TRUE)
    expect_match(plain, "does not certify package-wide coverage", fixed = TRUE)
    expect_match(plain, "cross_validate", fixed = TRUE)
  }
})

test_that("NEWS labels superseded public claims in their historical sections", {
  root <- testthat::test_path("..", "..")
  news <- paste(readLines(file.path(root, "NEWS.md"), warn = FALSE), collapse = "\n")
  news_plain <- gsub("[[:space:]]+", " ", news)
  forbidden_proximity <- paste(
    "bace_final_imp.{0,240}(proper multiple-imputation|Rubin \\(1987\\))",
    "PMM.{0,240}(pool_mi|Rubin|honest standard errors|properly calibrated SEs)",
    "tree uncertainty.{0,120}downstream standard errors.{0,120}Rubin",
    "pooled via Rubin's rules.{0,120}default path|Rubin's rules stable",
    sep = "|"
  )
  has_forbidden_proximity <- function(text) {
    grepl(forbidden_proximity, gsub("[[:space:]]+", " ", text),
      ignore.case = TRUE, perl = TRUE
    )
  }

  expect_false(grepl("net win-loss vs BACE|outperformed BACE|competitive with BACE|beating BACE|pigauto is uniquely positioned|Why this is novel|Coverage is guaranteed under exchangeability", news_plain, ignore.case = TRUE))
  expect_false(grepl("accuracy switch|0\\.14-1\\.27 z-RMSE", news_plain, ignore.case = TRUE))
  expect_false(grepl("Where PMM genuinely helps", news_plain, ignore.case = TRUE))
  expect_false(has_forbidden_proximity(news))
  expect_true(has_forbidden_proximity(
    "bace_final_imp() produces\nproper multiple-imputation draws for a historical route"
  ))
  expect_true(has_forbidden_proximity(
    "PMM draws can be passed to\npool_mi() for Rubin inference"
  ))
  expect_true(has_forbidden_proximity(
    "tree uncertainty can propagate into\ndownstream standard errors through Rubin pooling"
  ))
  expect_true(has_forbidden_proximity(
    "datasets are pooled via Rubin's rules\nand become the default path"
  ))
  expect_match(news, "Historical comparator snapshot — withdrawn as current evidence", fixed = TRUE)
  expect_match(news, "Historical rationale — novelty not established", fixed = TRUE)
  expect_match(news, "Historical tree workflow — now unsupported", fixed = TRUE)
  expect_match(news, "Historical conformal description — superseded", fixed = TRUE)
  expect_match(news, "external comparator gate remains open", fixed = TRUE)
  expect_match(news, "opt-in routing with fallback", fixed = TRUE)
  expect_match(news, "Historical BACE bridge note — removed and not inference-validated", fixed = TRUE)
  expect_match(news, "Historical PMM rationale — superseded", fixed = TRUE)
  expect_match(news, "multi_impute_analysis()", fixed = TRUE)
})
