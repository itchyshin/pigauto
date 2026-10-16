## 1. Goal

Refresh the bundled BirdTree data-rights evidence from current primary sources and reconcile the release summary with the recorded maintainer warranty.

## 2. Implemented

- Added a dated BirdTree and CRAN-policy source recheck to `provenance/rights-and-policy.md`.
- Clarified in the current release summary that the maintainer-warranty basis is recorded, while the final post-merge tarball must still reconcile the exact tree objects and `inst/NOTICE`.
- Updated the machine-readable release ledger to distinguish the historical pre-PR #234 tarball inspection from the required final-artifact recheck.

## 3a. Decisions and Rejected Alternatives

- Preserved the user's explicit maintainer confirmation that the bundled example trees may be included with citation and credit. The official BirdTree page supplies citation instructions but no separate redistribution licence; this is not described as a BirdTree-issued licence or independent legal opinion.
- Did not remove the bundled trees or add a downloader. The requested scope is an audit, and the recorded CRAN maintainer warranty is the current inclusion basis.
- Assumptions made without asking: none.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-16-rights-refresh.md`

## 5. Checks Run

- Opened the official BirdTree downloads page and verified its current citation instructions for full/partial tree use and web-tool use. It does not state a separate redistribution licence.
- Opened the current CRAN Repository Policy and verified the source-package requirements for clear rights in all components, including data, and the maintainer warranty for credited material.
- Parsed `release-ledger.json` with Python's standard JSON parser; it returned `RELEASE_LEDGER_JSON_OK`.
- `check-after-task.R` passed report structure; its embedded acceptance scan could not resolve Node from its PATH, so I ran Unlazy directly with the bundled Node. The scoped ledger reports 4 met and 3 unmet. `slop_check.py` found 0 hits in 634 words, and `git diff --check` passed.
- The active PR #234 checks were freshly inspected in Chrome; Ubuntu R-release, Ubuntu R-devel, and macOS R-release remain queued. The R-release job says it is waiting for a runner.

## 6. Tests of the Tests

No implementation tests changed. The rights claims are tied directly to the primary-source page content; the artifact-specific object and notice checks remain pending until the final post-merge tarball exists.

## 7a. Issue Ledger

- Resolved at the maintainer-warranty level: Shinichi's explicit inclusion confirmation is recorded and remains the release basis.
- Open: inspect the exact tree objects and `inst/NOTICE` in the final post-merge tarball and bind those checks to its checksum.
- Open: final source merge, deployed-site verification, exact artifact platform checks, and independent review.

## 8. Consistency Audit

The live summary, detailed rights note, and JSON release ledger now distinguish three things: BirdTree's citation instructions, the maintainer's inclusion warranty, and the final tarball's object/notice verification. The old tarball hash is explicitly historical because it predates PR #234. No public-facing package wording or tree object changed.

## 9. What Did Not Go Smoothly

The first lane-preflight invocation lacked the bundled Git path and falsely reported that the checkout was not a repository. With the fallback Git selected, preflight found the PR #234 Codex lane; all listed leases were expired. The lease registry lock could not be acquired under the default sandbox. A narrow two-hour lease was then granted using the approved escalated coordination write. No overlapping worktree files were edited.

## 10. Known Residuals

This recheck does not independently establish third-party rights. It records the maintainer warranty required by CRAN and the page's citation requirements. Exact shipped-object and NOTICE verification is still required on the final tarball. PR #234's checks are queued; no merge, deployment, submission, or external contact occurred.

## 11. Team Learning

Prior team findings already distinguished the `megatrees` software licence from BirdTree data rights. This continuation refreshed that evidence against the current official pages and corrected the live ledger summary. Memory receipt: the repository LOAD-FIRST manifest and prior pigauto release notes shaped the scope. Golden Set: not checked; no code changed. No brain-vault notes were edited.

## 12. Cross-Product Coverage

Covers: provenance and rights wording in the CRAN release ledger for the bundled BirdTree examples.

Does NOT cover: independent legal clearance, a final post-merge artifact, platform checks, deployment, final independent review, or CRAN submission.
