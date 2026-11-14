test_that("resource boundary runs in an owned source-bound child", {
  skip_if(!Sys.info()[["sysname"]] %in% c("Darwin", "Linux"),
          "resource child requires Linux or macOS")
  skip_if(!file.exists("/usr/bin/time"), "resource child requires /usr/bin/time")
  skip_if(!nzchar(Sys.which("ps")), "resource child requires ps")
  skip_if(!file.exists(file.path(R.home("bin"), "Rscript")), "resource child requires Rscript")
  skip_if(!exists("mem.maxVSize", envir = baseenv(), mode = "function", inherits = FALSE),
          "resource child requires mem.maxVSize")
  skip_if(!"timeout" %in% names(formals(system2)), "resource child requires system2 timeout")
  skip_if(identical(Sys.info()[["sysname"]], "Darwin") &&
            (!nzchar(Sys.which("vm_stat")) || !nzchar(Sys.which("sysctl"))),
          "resource child requires macOS memory tools")
  root <- normalizePath(getNamespaceInfo(asNamespace("pigauto"), "path"))
  fixture <- testthat::test_path("fixtures", "resource-child-test.R")
  driver <- testthat::test_path("fixtures", "resource-child-driver.R")
  receipt <- tempfile(fileext = ".rds"); expected_file <- tempfile(fileext = ".rds")
  on.exit(unlink(c(receipt, expected_file)), add = TRUE)
  # Conservative availability check before allocating the fixture, in either process.
  headroom_unknown <- function() stop(errorCondition("memory headroom unavailable", class = "resource_headroom_unknown"))
  safe_read <- function(path) {
    if (!file.exists(path)) headroom_unknown()
    tryCatch(readLines(path, warn = FALSE), error = function(e) headroom_unknown())
  }
  available <- tryCatch({ if (identical(Sys.info()[["sysname"]], "Linux")) {
    x <- safe_read("/proc/meminfo")
    host <- suppressWarnings(as.numeric(sub(".*:[[:space:]]*([0-9]+).*", "\\1", grep("^MemAvailable:", x, value = TRUE)))) * 1024
    if (length(host) != 1L || !is.finite(host) || host < 0) headroom_unknown()
    groups <- safe_read("/proc/self/cgroup")
    v2 <- grep("^0::", groups, value = TRUE)
    v1 <- grep("^[0-9]+:[^:]*memory[^:]*:", groups, value = TRUE)
    if (length(v2)) {
      base <- "/sys/fs/cgroup"
      rel <- sub("^0::", "", v2[[1]])
      limit_name <- "memory.max"; usage_name <- "memory.current"
    } else if (length(v1)) {
      base <- "/sys/fs/cgroup/memory"
      rel <- sub("^[^:]*:[^:]*:", "", v1[[1]])
      limit_name <- "memory.limit_in_bytes"; usage_name <- "memory.usage_in_bytes"
    } else {
      base <- NULL
    }
    if (!is.null(base)) {
      if (grepl("(^|/)\\.\\.(/|$)", rel)) headroom_unknown()
      dir <- file.path(base, sub("^/", "", rel))
      base <- tryCatch(normalizePath(base, mustWork = TRUE), error = function(e) headroom_unknown())
      dir <- tryCatch(normalizePath(dir, mustWork = TRUE), error = function(e) headroom_unknown())
      if (!(identical(dir, base) || startsWith(dir, paste0(base, "/")))) headroom_unknown()
      repeat {
        limit_path <- file.path(dir, limit_name); usage_path <- file.path(dir, usage_name)
        # The cgroup v2 root may omit controller limit files; its bound is
        # physical memory, already represented by MemAvailable above.
        if (identical(dir, base) && length(v2) &&
            !file.exists(limit_path) && !file.exists(usage_path)) break
        if (!file.exists(limit_path) || !file.exists(usage_path)) headroom_unknown()
        limit_text <- safe_read(limit_path)[1]
        usage <- as.numeric(safe_read(usage_path)[1])
        limit <- if (identical(limit_text, "max")) Inf else as.numeric(limit_text)
        if (is.na(limit) || !is.finite(usage)) headroom_unknown()
        host <- min(host, limit - usage)
        if (identical(dir, base)) break
        dir <- dirname(dir)
      }
    }
    host
  } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
    x <- system2("vm_stat", stdout = TRUE); page <- as.numeric(system2("sysctl", c("-n", "hw.pagesize"), stdout = TRUE))
    n <- as.numeric(gsub("[^0-9]", "", grep("^Pages (free|inactive|speculative):", x, value = TRUE)))
    if (length(n) != 3L || any(!is.finite(n)) || length(page) != 1L ||
        !is.finite(page) || page <= 0) headroom_unknown()
    sum(n) * page
  } else headroom_unknown() }, resource_headroom_unknown = function(e) NA_real_)
  skip_if(length(available) != 1L || !is.finite(available), "resource child cannot verify memory headroom")
  cat("RESOURCE_ADMISSION_AVAILABLE_BYTES", available, "\n")
  skip_if(available < 2 * 1024^3, "resource child insufficient memory headroom")
  fn <- getFromNamespace("exact_conditional_mvn", "pigauto")
  src <- getSrcFilename(fn, full.names = TRUE)
  skip_if(length(src) != 1L || is.na(src) || !nzchar(src) || !file.exists(src),
          "resource child source fingerprint unavailable")
  body_file <- tempfile(); on.exit(unlink(body_file), add = TRUE)
  writeLines(deparse(body(fn), width.cutoff = 500L), body_file)
  expected <- list(namespace = root, source = normalizePath(src), source_md5 = unname(tools::md5sum(src)), body_md5 = unname(tools::md5sum(body_file)))
  saveRDS(expected, expected_file)
  parent_limit <- mem.maxVSize()
  thread_names <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
  parent_threads <- Sys.getenv(thread_names, unset = NA_character_)
  # The cap is on the child vector heap, not its total RSS. The R timeout
  # terminates the child process group; time records peak RSS.
  time_flag <- if (identical(Sys.info()[["sysname"]], "Darwin")) "-l" else "-v"
  output <- withr::with_envvar(setNames(rep("1", length(thread_names)), thread_names),
    system2("/usr/bin/time", c(time_flag, shQuote(file.path(R.home("bin"), "Rscript")),
      shQuote(driver), shQuote(root), shQuote(fixture), shQuote(receipt), shQuote(expected_file)),
      stdout = TRUE, stderr = TRUE, timeout = 60))
  expect_identical(Sys.getenv(thread_names, unset = NA_character_), parent_threads)
  expect_identical(mem.maxVSize(), parent_limit)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  expect_true(file.exists(receipt), info = paste(output, collapse = "\n"))
  if (!file.exists(receipt)) stop("resource child did not produce a receipt")
  got <- readRDS(receipt)
  expect_identical(got$identity, expected)
  if (identical(got$status, "prerequisite") &&
      identical(got$reason, "resource child vector heap limit unavailable"))
    skip("resource child vector heap limit unavailable")
  expect_identical(got$status, "pass")
  expect_gt(got$passed, 0)
  expect_true(any(grepl("RESOURCE_LIMIT_RESTORED", output, fixed = TRUE)))
  expect_true(any(grepl("RESOURCE_HEAP_MIB", output, fixed = TRUE)))
  expect_true(any(grepl("RESOURCE_BOUNDED_RSS_KIB_BEFORE", output, fixed = TRUE)))
  expect_true(any(grepl("RESOURCE_BOUNDED_RSS_KIB_AFTER", output, fixed = TRUE)))
  expect_true(any(grepl("RESOURCE_HEALTHY_RSS_KIB_BEFORE", output, fixed = TRUE)))
  expect_true(any(grepl("RESOURCE_HEALTHY_RSS_KIB_AFTER", output, fixed = TRUE)))
  expect_true(any(grepl("RESOURCE_THREADS", output, fixed = TRUE)))
  expect_true(any(grepl("maximum resident set size|Maximum resident set size", output)))
  cat(paste(output, collapse = "\n"), "\n")
})
