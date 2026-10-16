# After-task: retain candidate-site identity and checker evidence

## 1. Goal

Close the evidence provenance follow-up for the local site candidate reviewed for PR #234. This does not establish deployment or release readiness.

## 2. Implemented

- Retained exact source-commit and generated pkgdown configuration identity files beside the audit receipt.
- Updated the candidate-site receipt with hashes and kept the complete standalone pkgdown checker output referenced.
- Corrected the visual-review timestamp to the observed local date and timezone.

## 3a. Decisions and Rejected Alternatives

- Reused the existing candidate build and its identity markers; did not rebuild the site.
- Kept G7 open because the candidate site has not been deployed.
- Assumptions made without asking: none.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/site-search-index-candidate-2026-10-16.md`
- `docs/dev-log/cran-0.11-audit/provenance/SOURCE_COMMIT-6ebc676.txt`
- `docs/dev-log/cran-0.11-audit/provenance/pkgdown-config-6ebc676.yml`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-16-site-evidence-provenance.md`

## 5. Checks Run

- SHA-256 values for both retained identity files and the standalone pkgdown check log match the receipt.
- Unlazy status reports 8 of 11 gates met; G7, G8, and G9 remain unmet.
- `git diff --check` passed before staging.
- Independent reviewer Rose verified the source commit, pkgdown configuration, and checker-log hashes and approved this bounded evidence follow-up.

## 6. Tests of the Tests

No test code changed. Rose independently recomputed the retained-file hashes and checked their contents against the receipt.

## 7a. Issue Ledger

- Resolved: candidate build identity and separate pkgdown checker output are retained and hash-bound.
- Open: PR checks, merge, deployed-site verification, final exact-tarball checks, and independent release-artifact review.

## 8. Consistency Audit

The receipt, retained identity markers, checker log, and G6 local visual review agree on source commit `6ebc67636f31fa8ff324ddd4a327d5cd808fba99` and pkgdown 2.2.1. Deployment remains explicitly separate.

## 9. What Did Not Go Smoothly

The initial PR check run stayed queued for over an hour with no runner time. It was canceled and a new attempt was started; the new attempt is queued and not yet evidence of passing checks.

## 10. Known Residuals

The candidate source has not merged or deployed. No exact post-merge source tarball, checksum-bound platform results, or final independent artifact verdict exists yet. CRAN submission is not authorized.

## 11. Team Learning

Store generated build identity markers with the receipt, not only in a temporary build directory. A passing local website check does not establish deployment. Golden Set: not checked. No brain-vault notes were edited.

## 12. Cross-Product Coverage

Covers: provenance for one local candidate pkgdown build and its standalone checker output.

Does NOT cover: the public deployed site, all reader surfaces, the post-merge tarball, Windows/platform checks, or CRAN submission.
