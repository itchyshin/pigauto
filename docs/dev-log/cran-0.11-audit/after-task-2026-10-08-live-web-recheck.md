## 1. Goal

Record a fresh live-site and committed-source comparison that identifies which reader-facing discrepancies remain before freezing the CRAN 0.11 artifact.

## 2. Implemented

Added a cache-busted page inventory with source identities, page-level discrepancies, reviewer findings, and explicit limits on what the check proves.

## 3a. Decisions and Rejected Alternatives

Kept the rights gate open. Citation instructions and the `megatrees` package license do not establish redistribution terms for BirdTree derivatives in pigauto. Tree retrieval guidance should describe a user-initiated workflow and should not imply that pigauto downloads trees.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/live-web-recheck-2026-10-08.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-live-web-recheck.md`

## 5. Checks Run

- Ran lane preflight in the original checkout and the gate-evidence worktree. The original checkout had unrelated dirty files and was not edited.
- Verified the gate worktree branch and HEAD: `release/cran-0.11-gate` at `94512fcee26a205c466fe59d15cdba641b21d679`; its only pre-existing dirty path was another lane's `retirement-manifest.md`.
- Compared cache-busted deployed Home, Getting Started, Tree Sensitivity, Multiple Imputation, and CRAN package pages with `origin/main` source at `0b0f71fee838c6ed51ef832ed819270eeafaf29b`.
- Requested independent read-only checks from the site auditor and tree-provenance reviewer. Both confirmed material mismatches; their scoped findings are recorded in the evidence note.
- `slop_check.py` reported zero findings for both new files. `git diff --no-index --check /dev/null <file>` emitted no whitespace diagnostics for either ignored new file; exit 1 is expected because each file differs from `/dev/null`. The after-task validator passed its structural check, then exited nonzero because five unmet gates belong to the separate `.unlazy/imputation-sim` pipeline; I did not alter another lane's ledger.
- GitHub API/SSH access failed because this session could not resolve or connect to GitHub. No remote state was changed.

## 6. Tests of the Tests

No package code changed. The key control was comparing live pages to committed source rather than treating uncommitted candidate-worktree edits as deployed evidence. The rights check separately distinguishes citation, package-license metadata, and permission to redistribute the underlying data.

## 7a. Issue Ledger

- Fixed in this slice: created a source-linked record of current live documentation discrepancies and their limits.
- Open: correct `tree300`'s MCC wording; make candidate-version status explicit; add reviewed user-initiated BirdTree retrieval guidance; reconcile generated help and provenance across the source lane; rebuild and inspect the site.
- Open: establish the redistribution basis for bundled BirdTree-derived data.
- Open: final exact-tarball Windows/platform checks and release closure.

## 8. Consistency Audit

The live home page advertises a 0.11 development version, while the CRAN listing remains 0.10.0. The Getting Started page's MCC description conflicts with committed package data provenance. The live tree article has no retrieval example, consistent with current committed source. The MI page's `auto` description matches current resolver documentation. The release/audit source branch has uncommitted help changes and predates PRs #229 and #230, so neither those edits nor that branch alone can establish a final candidate.

## 9. What Did Not Go Smoothly

The in-app network path to GitHub was unavailable, so current PR #228 state could not be freshly verified. The cached public page reader initially served older pages; the independent browser audit used cache-busted URLs and saw the current 0.11.0 version badges. An initial lane-lease call could not write the protected registry and was not relied on; the exact-path lease succeeded after escalation.

## 10. Known Residuals

This slice does not correct the source pages, establish deployed commit identity, resolve data rights, validate a fresh pkgdown build, or validate a frozen artifact. The gate PR remains unmerged by user instruction.

## 11. Team Learning

Read the candidate badge, the CRAN package listing, the article source commit, and the deployed article together. A current home page can link to an old article, and a branch's uncommitted correction is not evidence of a deployed fix.

Memory receipt: loaded the pigauto route manifest and ran lane preflight; the exact-artifact and lane-ownership rules shaped this review. Brain content was not changed.

Golden Set: not checked; this was a reader-surface release audit, not a recurring product defect.

## 12. Cross-Product Coverage

Covers: selected live reader pages, their current committed vignette/source counterparts, the CRAN listing, and the mixed-backbone provenance discrepancy identified by two independent reviews.

Does NOT cover: complete inventory of all URLs/assets, a fresh site build, visual browser inspection by the primary reviewer, deployment commit verification, final source integration, package tests, platform checks, or submission readiness.
