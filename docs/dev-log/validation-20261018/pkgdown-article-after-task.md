## 1. Goal

Add a public pkgdown article that reports the measured results and limitations of the October 2026 pigauto validation campaign, with direct links to the retained evidence.

## 2. Implemented

Added a source-backed vignette and a menu entry. The article reports run accounting, measured prediction coverage, inferential classifications, package-check scope, and reproducibility identifiers. Links target the immutable validation commit.

## 3a. Decisions and Rejected Alternatives

Kept the current main-based pkgdown destination and navigation structure. The competing unmerged navbar changes included old benchmark pages or broader content that was not necessary for this article, so they were not merged. Used a vignette rather than a static HTML snapshot so pkgdown builds the page from source.

Assumptions made without asking: Shinichi's approval to proceed covers a public article, PR, merge after checks, and deployment verification; numerical claims remain limited to the retained campaign report and tables.

## 4. Files Touched

- `_pkgdown.yml`
- `vignettes/validation-results.Rmd`
- `docs/dev-log/validation-20261018/assessments/validation-results.assessment.json`
- `docs/dev-log/validation-20261018/pkgdown-article-after-task.md`

## 5. Checks Run

- `pkgdown::build_article("validation-results")`: passed and wrote `_site/articles/validation-results.html`.
- `Rscript pkgdown/clean-internal-pages.R`: passed; the retired simulation-study route remained absent.
- `pkgdown::check_pkgdown()`: passed, “No problems found.”
- `git diff --check`: passed.
- Exact article checks: generated page, articles index, navbar target, immutable commit links, and retired simulation-study absence all passed.
- `python3 ~/shinichi-brain/tools/slop_check.py <absolute article path>`: zero findings.
- `benchmark.py assessment-check ... --mode independent --version v1`: `ASSESSMENT_RECORDED` for SHA-256 `b59108ce5a2513d0a0fe642cc02329e878960c776295a63406af683de2514091`.
- A full `pkgdown::build_site()` attempt did not complete because the sandbox could not resolve `cloud.r-project.org`; a first retry also exposed an invalid `override = TRUE` invocation, corrected to the documented no-argument call. The article and source checks pass, but full-site generation remains for CI.

## 6. Tests of the Tests

This is documentation-only. No package test was added. The article build and pkgdown source checker verified that the new source and route resolve; an independent script checked the rendered page and retired-route handling. No deliberately broken vignette was introduced.

## 7a. Issue Ledger

- Fixed: initial local full-site command used an invalid `override = TRUE` argument; reran with `pkgdown::build_site()`.
- Deferred to CI: full-site build needs network access to resolve the CRAN package link.
- Open: pull request checks and main-branch deployment have not yet run.

## 8. Consistency Audit

The article's attempt counts, result-row totals, inferential classifications, and four prediction-coverage averages were recomputed from the retained CSV receipts. The posterior slope ranges and analysis-aware classifications match the same summary table. Existing experimental-status wording and the `_site` build destination remain intact. The public article menu points to the rendered route, and the detailed report, manifest, ledger, and summary links resolve to files at the validation commit.

## 9. What Did Not Go Smoothly

The shared checkout was dirty, so the change was prepared in an isolated temporary clone. A full local site build stopped at a network-dependent CRAN lookup. The existing validation branch worktree is outside the writable project roots, so the final change was prepared in the temporary clone instead of modifying that worktree.

## 10. Known Residuals

The full site has not completed locally, GitHub CI has not run for this change, and the public page is not live until merge and a successful Pages deployment. No new scientific validation was performed; this article summarizes the existing campaign.

## 11. Team Learning

Memory receipt: the pigauto `LOAD-FIRST` route, repository guidance, and pkgdown main-only deployment guidance shaped the scope. The current main configuration is the correct base for the article; unmerged navbar experiments were reviewed and left separate.

Golden Set: not run because this documentation-only change did not modify package behavior or a known code-mistake class.

## 12. Cross-Product Coverage

Covers: vignette source, pkgdown article index, curated article-menu route, and links to retained validation evidence.

Does NOT cover: package runtime behavior, unsupported inferential regimes, unavailable real-data cases, local full-site generation without network access, or live deployment before merge.
