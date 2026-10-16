# After-task: candidate-site visual review

## 1. Goal

Visually review the current local site candidate for the pigauto 0.11 audit in Chrome and record whether the local-site gate can pass.

## 2. Implemented

- Opened the full candidate site built from PR #234 source commit `6ebc67636f31fa8ff324ddd4a327d5cd808fba99` and pkgdown 2.2.1 through a loopback-only HTTP server.
- Visually reviewed the home page, Getting Started, multiple-imputation, tree-uncertainty, and `tree300` reference pages.
- Marked local visual gate G6 met for this candidate build and kept deployed-site G7 open.

## 3a. Decisions and Rejected Alternatives

- Used the existing full candidate site and its recorded source/build identity. Did not rebuild it.
- Treated browser-rendered content and screenshots as local visual evidence, not deployment evidence.
- Assumptions made without asking: none.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-16-candidate-site-visual-review.md`

## 5. Checks Run

- Chrome rendered the candidate home page, Getting Started, multiple-imputation, tree-uncertainty, and `tree300` help pages from the local candidate site.
- Visual inspection found no clipped content or broken navigation. Code blocks, tables, headings, and the table of contents fit their rendered columns.
- Rendered content showed the seven original traits and 14 latent columns, partial-overlap behavior, BirdTree citations, optional drmTMB/gllvmTMB install instructions, fixed-effect-only pooling limits, and the inverse-Wishart prior marked as legacy-only.
- The retained candidate receipt binds the site to commit `6ebc67636f31fa8ff324ddd4a327d5cd808fba99`, pkgdown 2.2.1, and the structural verifier results.
- Unlazy reports G6 met; G7, G8, and G9 remain unmet.

## 6. Tests of the Tests

No automated test changed. The separate structural site verifier and crawler remain the machine-checkable checks; this review supplies the manual rendered-page inspection G6 requires.

## 7a. Issue Ledger

- Resolved: local visual review for the exact recorded candidate build.
- Open: verify the post-merge deployment and live retained, retired, search, and sitemap routes (G7).
- Open: freeze and run all required checks on one final post-merge tarball (G8).
- Open: obtain independent final-artifact review and prepare the unmerged release-evidence PR (G9).

## 8. Consistency Audit

The local screenshot review agrees with the candidate-site receipt and rendered content. The candidate site is explicitly distinguished from the public deployment, which still requires a later live check.

## 9. What Did Not Go Smoothly

The default sandbox prevented binding a localhost port. The approved local visual-review workflow was completed after requesting the scoped escalation to serve the already-built static candidate on loopback. No public network service was opened.

## 10. Known Residuals

The source PR remains open and its R checks were queued at the latest Chrome check. The visual review does not establish merge, Pages deployment, final tarball readiness, external platform checks, or CRAN readiness.

## 11. Team Learning

Using a loopback-only server allowed Chrome to inspect a local static candidate without enabling file-URL access. Prior independent review found the candidate verifier and 567-path count sound; it did not substitute for this visual check. Golden Set: not checked. No brain-vault notes were edited.

## 12. Cross-Product Coverage

Covers: local visual presentation of five primary reader pages from the candidate site built from PR #234.

Does NOT cover: every rendered page pixel-by-pixel, the public deployed site, post-merge source, final package tarball, platform checks, independent exact-artifact verdicts, or CRAN submission.
