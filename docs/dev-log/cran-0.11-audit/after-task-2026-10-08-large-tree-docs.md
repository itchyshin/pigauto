## 1. Goal

Correct the large-tree guidance in Getting Started so it matches the implemented graph-building paths and does not promise unsupported laptop performance.

## 2. Implemented

Updated the large-tree section to explain the optional Lanczos solver above 7,500 tips, its dense fallback, the dense matrices used during graph construction, and the limited meaning of the 32-feature `k_eigen = "auto"` cap. Added the candidate-source finding and local verification to the site-review record.

## 3a. Decisions and Rejected Alternatives

Kept the documented solver threshold and fallback behavior because they match `R/build_phylo_graph.R`. Removed the standard-laptop claim and the statement that sparse Lanczos is future work. Did not add a new performance benchmark or promise a runtime for any tree size.

## 4. Files Touched

- `vignettes/getting-started.Rmd`
- `docs/dev-log/cran-0.11-audit/site-review-2026-10-08.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-large-tree-docs.md`

## 5. Checks Run

- Candidate source inspection of `R/build_phylo_graph.R`, `DESCRIPTION`, and the sparse-versus-dense graph test: confirmed `RSpectra` is optional, the threshold is 7,500 tips, failures fall back to dense `eigen()`, and the graph path creates dense N-by-N intermediates.
- Fresh offline pkgdown build: completed and wrote the corrected Getting Started article. The shell wrapper then reported a zsh error because it assigned to the reserved variable `status`; the build log ends with `Finished building pkgdown site for package pigauto`. This shell error did not invalidate the completed build.
- `Rscript --vanilla pkgdown/clean-internal-pages.R` followed by `python3 script/cran-0.11-site/crawl.py _site`: 62 HTML pages, 3,443 local references, 34 retired routes, 0 errors; `SITE_CRAWL_OK`.
- Rendered Getting Started assertions: passed for the 7,500-tip optional Lanczos description, dense-matrix memory warning, and absence of the old laptop and future-Lanczos claims.
- `Rscript --vanilla -e 'pkgdown::check_pkgdown()'`: no problems found.
- Search of reader-facing source and generated help for the removed unsupported phrases: no remaining occurrences.
- `git diff --check`: passed.
- `closeout.py check` and the closeout-mode after-task validator do not pass while the wider approved audit ledger still has open gates; the structure check and evidence-promotion check pass. The release remains in progress, so this narrow documentation slice does not claim whole-goal completion.

## 6. Tests of the Tests

The crawler checked all generated local references and retired routes. Rendered-content assertions read the built HTML rather than relying on source inspection alone. The assertion was adjusted after the first attempt expected contiguous source wording where HTML formatting split the phrase; the final assertions passed.

## 7a. Issue Ledger

- Fixed: candidate Getting Started described an obsolete universal dense eigensolver and a laptop-scale 10,000-tip test claim.
- Open: visual review of the local rendered site, post-deployment page checks, BirdTree redistribution basis, and exact-artifact release gates.

## 8. Consistency Audit

Compared the vignette with the implementation, optional-dependency declaration, and existing graph tests. Searched README, vignettes, NEWS, installed help, and NOTICE for the removed laptop-performance and future-Lanczos statements. No further reader-facing copy using those claims was found. The deployed Getting Started page still needs checking after deployment.

## 9. What Did Not Go Smoothly

The shell wrapper used `status`, a reserved zsh variable, after pkgdown had already completed. I verified the completed build from its log and ran the cleanup, crawler, rendered assertions, and pkgdown consistency check separately. The first rendered assertion was too strict about inline HTML formatting and was replaced with checks that match the rendered structure.

## 10. Known Residuals

The fresh build used the approved offline system-font and metadata overrides. Chrome has not visually inspected this new local rendering because local file URLs are blocked by the browser policy. The live site still shows the old copy until deployment. The global closeout gate remains red while release gates are open. No source merge, deployment, frozen tarball, or CRAN submission occurred.

## 11. Team Learning

Memory receipt: loaded the pigauto project manifest and relevant CRAN audit memory. No team member was dispatched for this narrow copy correction; prior independent reviews remain recorded in the audit files. Golden Set: not in scope for this documentation correction.

## 12. Cross-Product Coverage

Covered: implementation details, optional `RSpectra` dependency behavior, Getting Started source, rendered article, site links, retired routes, and pkgdown consistency. Does NOT cover: measured large-tree performance, visual review, live deployment, BirdTree redistribution rights, or exact-tarball release checks.
