#!/usr/bin/env Rscript

# Bounded audit of current exported defaults and their effective paths.
# Run from the pigauto repository root. No posterior sampler or GNN is run.

if (!file.exists("DESCRIPTION") ||
    !identical(read.dcf("DESCRIPTION", fields = "Package")[[1L]],
               "pigauto")) {
  stop("Run this script from the pigauto repository root.", call. = FALSE)
}

results <- devtools::test(filter = "cran-audit-defaults",
                          stop_on_failure = TRUE)
counts <- as.data.frame(results)
if (nrow(counts) == 0L ||
    any(counts$failed > 0L | counts$error | counts$warning > 0L |
        counts$skipped)) {
  stop("The defaults audit did not finish with clean test results.",
       call. = FALSE)
}
cat("CRAN_DEFAULTS_AUDIT_OK\n")
