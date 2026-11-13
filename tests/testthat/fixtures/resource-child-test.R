test_that("bounded real allocation preserves exact conditional safety", {
  exact <- getFromNamespace("exact_conditional_mvn", "pigauto")
  henderson <- getFromNamespace("build_henderson_S_inv", "pigauto")
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  h <- henderson(tree)
  small <- cbind(x = c(1, NA, -1, .5), y = c(NA, .2, .3, .4))
  rownames(small) <- tree$tip.label
  control <- exact(small, diag(2), h, eps = 0)
  expect_false(is.null(control))
  expect_equal(unname(control$mu[is.na(small)]), c(.5, .1), tolerance = 1e-8)
  expect_equal(unname(control$var[is.na(small)]), c(.75, .75), tolerance = 1e-8)

  K <- 1000L
  L <- matrix(0, 4L, K); L[2, 1] <- NA_real_
  S <- matrix(.1, K, K); diag(S) <- 1.1
  expected_var <- .75 / (1 - .1 / (1 + .1 * K))
  check <- function(x) {
    expect_false(is.null(x))
    expect_identical(dim(x$mu), dim(L)); expect_identical(dim(x$var), dim(L))
    expect_true(all(is.finite(x$mu)) && all(is.finite(x$var)))
    expect_equal(x$mu[2, 1], 0, tolerance = 1e-8)
    expect_equal(x$var[2, 1], expected_var, tolerance = 1e-7)
    expect_true(all(x$mu[!is.na(L)] == 0) && all(x$var[!is.na(L)] == 0))
  }
  rss <- function() as.numeric(system2("ps", c("-o", "rss=", "-p", Sys.getpid()), stdout = TRUE))
  bounded <- function() {
    old <- mem.maxVSize()
    on.exit({
      restored <- mem.maxVSize(old)
      cat("RESOURCE_LIMIT_RESTORED", restored, "EXPECTED", old, "\n")
      if (!identical(restored, old)) stop("vector heap limit restoration failed")
    }, add = TRUE)
    g <- gc()
    used <- g["Vcells", "used"] * 8 / 1024^2
    trigger <- g["Vcells", "gc trigger"] * 8 / 1024^2
    target <- ceiling(trigger) + 128
    if (is.finite(old) && old < target) skip("resource child vector heap limit unavailable")
    active <- withCallingHandlers(mem.maxVSize(target), warning = function(w) {
      if (grepl("a limit lower than current usage, so ignored", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
    cat("RESOURCE_HEAP_MIB used", used, "trigger", trigger, "old", old,
        "target", target, "active", active, "\n")
    if (!identical(active, as.double(target))) skip("resource child vector heap limit unavailable")
    cat("RESOURCE_BOUNDED_RSS_KIB_BEFORE", rss(), "\n")
    out <- exact(L, S, h, eps = 0, max_cells = 5000L)
    cat("RESOURCE_BOUNDED_RSS_KIB_AFTER", rss(), "\n")
    if (!is.null(out)) check(out)
    cat("RESOURCE_BOUNDED", if (is.null(out)) "NULL" else "FINITE", "\n")
  }
  bounded()
  cat("RESOURCE_HEALTHY_RSS_KIB_BEFORE", rss(), "\n")
  healthy <- exact(L, S, h, eps = 0, max_cells = 5000L)
  cat("RESOURCE_HEALTHY_RSS_KIB_AFTER", rss(), "\n")
  check(healthy)
})
