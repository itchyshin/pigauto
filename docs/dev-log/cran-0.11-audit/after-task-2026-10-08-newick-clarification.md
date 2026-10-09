# BirdTree file-format clarification: after-task report

## 1. Goal

Correct the last reviewed mismatch between BirdTree's download format and the file formats accepted by pigauto, then verify the rendered reader pages and route checks.

## 2. Implemented

Changed the README to state that BirdTree distributes Newick trees and that `read_tree()` also accepts NEXUS files from other sources. Updated the 0.11.0 NEWS entry to record the distinction. Added the fresh review and build evidence to the release gate ledger. G4 is now met.

## 3a. Decisions and Rejected Alternatives

- Kept `read_tree()`'s documented NEXUS support. It describes pigauto's parser, not BirdTree's download format.
- Did not add NEXUS download instructions to BirdTree's website guidance because the official BirdTree FAQ says its distributions are stored as Newick.
- Kept the rights, visual, deployment, and exact-artifact questions in their separate release gates.

## 4. Files Touched

- `README.md`
- `NEWS.md`
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-newick-clarification.md`

## 5. Checks Run

- Official [BirdTree FAQ](https://birdtree.org/faq/) says its distributions are stored as Newick trees.
- Independent documentation review confirmed the current README and NEWS source, and fresh rendered pages `_site/index.html` and `_site/news/index.html`, all distinguish Newick downloads from NEXUS parser support.
- Fresh offline pkgdown build: exit 0.
- `Rscript --vanilla pkgdown/clean-internal-pages.R`: exit 0.
- `python3 script/cran-0.11-site/crawl.py _site`: 62 HTML pages, 3,443 local references, 34 retired routes, zero errors; `SITE_CRAWL_OK`.
- `Rscript --vanilla -e 'pkgdown::check_pkgdown()'`: no problems found.
- Unlazy `--reverify --approve --timeout 600` on `GATES.md`: G5 and G5b both passed; 11 release gates parsed, 6 met and 5 unmet.
- Direct Chrome navigation to the local HTML site remains blocked by browser URL policy. No workaround was attempted.

## 6. Tests of the Tests

No code test was needed for this prose-only correction. The reviewer compared the source distinction with the official FAQ and checked the fresh rendered homepage and NEWS output. The route crawler also checked the entire generated site and reported zero errors.

## 7a. Issue Ledger

- Fixed: README conflated BirdTree's Newick downloads with `read_tree()`'s additional NEXUS support.
- Fixed: NEWS now records the corrected format distinction.
- Open: local visual inspection, merged deployment, BirdTree redistribution basis, final post-merge tarball checks, and independent exact-artifact review.

## 8. Consistency Audit

Searched README, NEWS, vignettes, R source, generated help, and NOTICE for BirdTree/NEXUS format claims. The generic README input instruction and `read_tree()` help already correctly describe both supported formats. The only provider-format conflation was in the BirdTree-specific README paragraph. The installed NOTICE continues to distinguish the `megatrees` software licence from rights to the underlying data.

## 9. What Did Not Go Smoothly

Chrome rejects local file URLs and prohibits local-server, alternate-browser, and indirect-access workarounds. The refreshed site was therefore checked structurally and by rendered-content inspection, while the visual gate remains open.

## 10. Known Residuals

G4, G5, and G5b are met. G0, G6, G7, G8, and G9 remain unmet. No merge, deployment, exact post-merge tarball check, or CRAN submission occurred. The BirdTree pages require attribution but do not document a redistribution grant for the bundled tree objects.

## 11. Team Learning

The independent review caught a provider-format claim that a generic parser description had made easy to miss. State the upstream format and the package's accepted formats separately, then verify both the source and rendered reader page.

Memory receipt: the pigauto LOAD-FIRST manifest informed the source-versus-artifact and evidence boundaries. Golden Set: not in scope for this documentation correction.

## 12. Cross-Product Coverage

- Covers: README source and rendered home page; NEWS source and rendered page; existing `read_tree()` source and generated help; local route, search, and sitemap checks.
- does NOT cover: BirdTree redistribution rights, local visual layout, deployed pages after merge, a frozen post-merge tarball, platform results for that tarball, or CRAN acceptance.
