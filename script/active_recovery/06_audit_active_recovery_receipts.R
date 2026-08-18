#!/usr/bin/env Rscript
# Usage: Rscript 06_audit_active_recovery_receipts.R pilot_dir extension_dir audit.rds
#
# A deliberately separate receipt audit for the completed Stage-A campaign.
# It checks campaign provenance and treatment application before the combined
# estimator is interpreted. It is not the campaign summary or claim verdict.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("expected: pilot_dir extension_dir audit.rds", call. = FALSE)
}

expected <- list(
  pilot = list(dir = args[[1L]], receipts = 800L, reps = 1:100L,
               sha = "beee5df360759f53cd22b595366436bca50845a5"),
  extension = list(dir = args[[2L]], receipts = 7200L, reps = 101:1000L,
                   sha = "ad3990c5bea62be7536d4e0ab2c81ea8f9ef686e")
)
policies <- c("active", "random", "uncertainty", "no_acquisition")

read_campaign <- function(spec, label) {
  files <- list.files(file.path(spec$dir, "results"), pattern = "\\.rds$",
                      full.names = TRUE)
  if (length(files) != spec$receipts) {
    stop(label, ": expected ", spec$receipts, " receipts, found ", length(files),
         call. = FALSE)
  }
  raw <- do.call(rbind, lapply(files, readRDS))
  if (nrow(raw) != 4L * spec$receipts) {
    stop(label, ": expected four policy rows per receipt", call. = FALSE)
  }
  if (!all(raw$status == "ok")) stop(label, ": receipt errors retained", call. = FALSE)
  if (!identical(sort(unique(raw$source_sha)), spec$sha)) {
    stop(label, ": unexpected source SHA", call. = FALSE)
  }
  key <- interaction(raw$n, raw$lambda, raw$family, raw$replicate, drop = TRUE)
  per_replicate <- split(raw, key)
  if (!all(vapply(per_replicate, function(x) identical(sort(x$policy), sort(policies)), logical(1)))) {
    stop(label, ": a replicate lacks exactly one row per policy", call. = FALSE)
  }
  if (!all(vapply(per_replicate, function(x) length(unique(x$replicate)) == 1L, logical(1)))) {
    stop(label, ": malformed replicate key", call. = FALSE)
  }
  expected_keys <- expand.grid(n = c(100L, 300L), lambda = c(1, 0.2),
                               family = c("continuous", "binary"), replicate = spec$reps,
                               stringsAsFactors = FALSE)
  observed_keys <- unique(raw[c("n", "lambda", "family", "replicate")])
  key_string <- function(x) do.call(paste, c(x, sep = "\r"))
  if (!setequal(key_string(observed_keys), key_string(expected_keys))) {
    stop(label, ": observed replicate grid differs from registered grid", call. = FALSE)
  }
  no_acquisition <- raw[raw$policy == "no_acquisition", , drop = FALSE]
  treated <- raw[raw$policy != "no_acquisition", , drop = FALSE]
  if (!all(no_acquisition$changed_hash == no_acquisition$initial_hash)) {
    stop(label, ": no-acquisition treatment changed data", call. = FALSE)
  }
  if (!all(treated$changed_hash != treated$initial_hash)) {
    stop(label, ": acquisition policy did not change data", call. = FALSE)
  }
  if (anyNA(treated$selected_species) || any(!nzchar(treated$selected_species))) {
    stop(label, ": acquisition policy lacks a selected species", call. = FALSE)
  }
  if (!all(is.finite(raw$primary))) stop(label, ": non-finite primary score", call. = FALSE)
  if (!all(raw$n_initial + raw$n_candidate + raw$n_test == raw$n)) {
    stop(label, ": split sizes do not cover the simulated species", call. = FALSE)
  }
  list(receipts = length(files), policy_rows = nrow(raw), raw = raw)
}

pilot <- read_campaign(expected$pilot, "pilot")
extension <- read_campaign(expected$extension, "extension")
audit <- list(
  audited_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  pilot_receipts = pilot$receipts,
  extension_receipts = extension$receipts,
  total_receipts = pilot$receipts + extension$receipts,
  total_policy_rows = pilot$policy_rows + extension$policy_rows,
  result = "pass"
)
saveRDS(audit, args[[3L]])
print(audit)
