## 1. Goal

Bind the local website checks and visual review to exact candidate source commit `77f858db920e3203e97f54dcf19baa1a098ea567` and record any reader-facing mismatches found during inspection.

## 2. Implemented

Removed a stale covariate-section placeholder from Getting Started and corrected two descriptions that called the seven-expression workflow six lines. Rebuilt the site from the exact candidate commit, reran structural retirement and link checks, and visually reviewed the updated Getting Started covariate section and `trees300` reference in the Codex In-app Browser. Updated the release gate ledger with this exact-head receipt.

## 3a. Decisions and Rejected Alternatives

Kept the website build and visual review separate from live deployment verification. Kept the BirdTree rights gate open because the displayed citations and upstream research-use guidance do not document the CRAN redistribution basis. Did not merge, deploy, freeze a new tarball, or submit to CRAN.

## 4. Files Touched

- `README.md`
- `vignettes/getting-started.Rmd`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- This after-task report

## 5. Checks Run

- `bash script/cran-0.11-site/build-release-site.sh`: exit 0 from exact source commit `77f858db920e3203e97f54dcf19baa1a098ea567`; fresh output at `/private/tmp/pigauto-cran011-site-77f858d.92JuKS`.
- Site receipt: 62 HTML pages, 613 search entries, 3,443 local references checked, 34 retired routes absent from output/search/sitemap, zero crawl errors, `SITE_CRAWL_OK`, and `pkgdown::check_pkgdown()` reported no problems.
- Build log SHA-256: `328158984d31fc3b2bf8c7ece727d501692e23ec509c587a3542fc9f703e2aac`; crawl log SHA-256: `2981f54a2056d62a71e37ea4f85eb8a8d6e22f242fc8e2d47d54daa1ab560ece`.
- Codex In-app Browser review of exact-head Getting Started and `trees300` rendered pages. Current covariate argument guidance is visible; the stale placeholder is absent. The tree page displays the mixed-backbone sample description, citations, and the licence boundary.
- `git diff --check`: passed after writing the report and ledger receipt.

## 6. Tests of the Tests

The retirement guard has a previously retained isolated negative control at `provenance/site-crawler-simulation-retirement-control.txt`. This continuation reran the positive build and crawler checks. No code test was changed.

## 7a. Issue Ledger

- Resolved: the stale Getting Started covariate placeholder and the incorrect six-line workflow descriptions.
- Resolved for this candidate source: fresh website build, structural crawl, retirement checks, and visual inspection of the two named pages.
- Open: G0 BirdTree redistribution basis; G7 live deployment and route verification; G8 exact post-merge tarball; G9 independent review of exact release evidence.

## 8. Consistency Audit

The visual review used files served from the build generated from the stated commit. Rendered Getting Started contains the corrected covariate guidance. Rendered `trees300` reports 28 Ericson and 22 Hackett trees and explicitly separates the `megatrees` software licence from rights to the BirdTree data. The ledger records 7 met and 4 open gates, with G0, G7, G8, and G9 still open.

## 9. What Did Not Go Smoothly

The first precise patch context did not match the long, historically accumulated gate ledger. No source was changed by that failed patch. I inspected the current ledger section and added a dated receipt at its end instead of rewriting prior evidence.

## 10. Known Residuals

This website review does not prove deployed state, exact tarball contents, platform checks for a frozen artifact, the BirdTree redistribution warranty, or CRAN acceptance. The local site builder emitted the known `interactive()` example-condition and Pandoc deprecation warnings. No layout issue was observed on the two reviewed pages; this is not a visual review of every generated page.

## 11. Team Learning

The single active pigauto lane was confirmed for this continuation. During exact-page review, check nearby headings and navigation prose as well as the targeted correction: that pass found both a stale placeholder and an inaccurate workflow line count. Preserve the distinction between right to use with attribution and the maintainer's separate CRAN redistribution warranty.

## 12. Cross-Product Coverage

Covers: exact-head local pkgdown build, structural links and retirement checks, rendered-content review, and visual inspection of Getting Started covariate guidance and the `trees300` reference.

Does NOT cover: all generated pages visually, deployed website state, exact post-merge tarball checks, unresolved data redistribution rights, CRAN acceptance, or new inferential claims.
