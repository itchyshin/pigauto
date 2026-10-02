# Byte-identity of converged campaign cells: new code vs 69670d4 (design.md 5e).
a <- commandArgs(trailingOnly = TRUE)   # old.rds new.rds
o <- readRDS(a[1]); n <- readRDS(a[2])
drop_t <- function(d) { if (is.data.frame(d)) d[, setdiff(names(d), c("wall_s")), drop = FALSE] else d }
cmp <- list(
  results = identical(drop_t(o$results), drop_t(n$results)),
  cell_detail = identical(o$cell_detail, n$cell_detail),
  cell_coverage = identical(drop_t(o$cell_coverage), drop_t(n$cell_coverage)),
  diagnostics_core = identical(o$diagnostics[, c("method","max_rhat","min_ess","converged")],
                               n$diagnostics[, c("method","max_rhat","min_ess","converged")]),
  n_missing = identical(c(o$n_missing_x, o$n_missing_y), c(n$n_missing_x, n$n_missing_y)))
print(unlist(cmp))
if (!isTRUE(cmp$results)) { print(all.equal(drop_t(o$results), drop_t(n$results))) }
cat(sprintf("new n_extensions: %s | code_sha old %s new %s\n", paste(n$diagnostics$n_extensions, collapse=","), substr(o$code_sha,1,10), substr(n$code_sha,1,10)))
cat(if (all(unlist(cmp))) "IDENTITY_OK\n" else "IDENTITY_DIFF\n")
