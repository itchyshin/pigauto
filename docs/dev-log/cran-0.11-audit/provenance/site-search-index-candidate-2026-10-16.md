# Candidate site search-index receipt, 2026-10-16

## Scope

This receipt checks a fresh full pkgdown site built from the exact source commit in PR #234. It resolves the candidate verifier's historical 568-path expectation and records the current deployment boundary. It is not evidence that the corrected site is deployed or that the final package tarball is ready.

## Inputs and identity

- Source commit: `6ebc67636f31fa8ff324ddd4a327d5cd808fba99` (`codex/cran-011-getting-started-correction`, PR #234).
- R: 4.6.0.
- pkgdown: 2.2.1, installed into `/private/tmp/pigauto-pkgdown-221-lib` from `pkgdown_2.2.1.tar.gz`; source tarball SHA-256 `9253a08bfbdddfa9f7f1404c5d1e5c9bc349845b078d06b68ef9288ece2db11c`.
- Candidate site: `/private/tmp/pigauto-cran011-site-6ebc676.OA2LOq/_site`.
- Retired-route manifest: `provenance/live-retired-routes-2026-10-09-v2.tsv`, SHA-256 `5414f6cb69e9d140c38d0a6308f8e0ddbe8952834eea01b6b861a78eecb807e1`.
- Candidate verifier and tests: `script/cran-0.11-site/verify_candidate_site.py`, SHA-256 `f922eef1b3f764d83afb751e92a438414bf916486ea103568271e951c1bcdd7a`; `test_candidate_site.py`, SHA-256 `cb05a428fa5dfc7840bac4a4bdfa8567faa5d65358af82cf0656fceacb146e59`.

## Commands and results

The exact PR #234 source snapshot was built with `PIGAUTO_SITE_CACHE_SOURCE=/private/tmp/pigauto-cran011-site-6ebc676.6hrGoz/_site R_LIBS=/private/tmp/pigauto-pkgdown-221-lib:/Users/z3437171/Library/R/arm64/4.6/library bash script/cran-0.11-site/build-release-site.sh`. The build completed, `pkgdown::check_pkgdown()` reported no problems, the crawl checked 3,584 local references with no errors, and the cleanup receipt reported 62 HTML pages and 613 search entries.

`python3 script/cran-0.11-site/test_candidate_site.py` passed all 16 tests. `python3 script/cran-0.11-site/verify_candidate_site.py --site-dir /private/tmp/pigauto-cran011-site-6ebc676.OA2LOq/_site --manifest docs/dev-log/cran-0.11-audit/provenance/live-retired-routes-2026-10-09-v2.tsv` returned `CANDIDATE_SITE_OK`, with 44 retired routes absent, 62 sitemap locations matching all 62 HTML files, 613 search records, and 567 non-empty search-path references. Every search path is present in the candidate sitemap, and no retired path appears in the output, sitemap, article index, or search index.

The passing test output is retained in `candidate-site-verifier-tests-2026-10-16.log`; the exact candidate output is in `candidate-site-verifier-2026-10-16-pkgdown-2.2.1.json.log`. The full build and crawl logs are retained beside this receipt. The earlier 2026-10-09 receipt reports 568 path references; that count was not reproduced with either pkgdown 2.2.0 or 2.2.1. The verifier and fixture now pin the observed 567 count. This count does not guarantee that every HTML page is independently searchable; five non-content or index pages have no search path.

The separately retained `pkgdown-check-2026-10-16.log` records `R=4.6.0`, `pkgdown=2.2.1`, and `✔ No problems found.` from a standalone `pkgdown::check_pkgdown()` run. The generated site's `SOURCE_COMMIT.txt` records the source commit above, and `_site/pkgdown.yml` records pkgdown 2.2.1. SHA-256 of the check log: `3d3dfaf7cc4dd264628f2b65a3b0953c815714aa2c43a91a1a21d179c857f946`.

Evidence-file SHA-256 values:

- build log: `a4643537f29fafa10a78386047a85000460a6b865079fa84236fc89ee55286c0`
- crawler log: `666b188ae633c31388a96b76c6b704ae575a69281e6606c220539e00f8f159d4`
- unit-test log: `4515cb569f2fd1cdb4b3e11c40a4c02697384c962de54d345466692c84723178`
- candidate-verifier output: `cff33162e02e900237f7a3110d6471956702e7333d14a2293875985840f908dd`
- standalone pkgdown check: `3d3dfaf7cc4dd264628f2b65a3b0953c815714aa2c43a91a1a21d179c857f946`

## Deployment state and limits

The Pages workflow log for run `38099234402` reports pkgdown 2.2.1 and deployed source `5e46c46b6970bf153a3b4bd00d55c9ace69225b0`. Chrome on 2026-10-16 still rendered the pre-correction Getting Started content (`X_scaled` shown as 300 × 4 although the example has seven traits and 14 latent columns). PR #234 remains open with three queued R checks. No local full-site visual inspection or post-correction live-site check has passed. This receipt supports the local candidate-site check only.
