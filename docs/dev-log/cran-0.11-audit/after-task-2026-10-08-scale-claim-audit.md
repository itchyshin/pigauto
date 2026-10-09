## 1. Goal

Align the package description and historical NEWS scaling claims with the saved run records before the CRAN 0.11 candidate advances.

## 2. Implemented

Removed the unqualified current-release claim that pigauto had been tested to 10,000 species. Corrected the v0.9.1 NEWS section to separate the completed 9,993-species missingness sweep from the v0.9.0 scaling curve that ended at 5,000 tips. Replaced the reference to a missing extended-scaling script with the existing driver, report, and log. A later independent review caught an inaccurate summary of the missingness results; the NEWS now reports the measured continuous and categorical comparisons and identifies them as one-run descriptive results without replicate-based uncertainty estimates.

## 3a. Decisions and Rejected Alternatives

Kept the 9,993-species result because `script/bench_avonet_missingness.md` and its log record all three missingness settings completing at that size. Kept the v0.9.0 timing result but scoped it to the measurements actually retained in `script/bench_scaling_v090.md`: a single replicate with four traits and 25% MCAR through n = 5,000. Removed the generic current-release ceiling instead of implying that a historical workload validates every 0.11.0 path. When review showed that NEWS said continuous pigauto RMSE beat BM at 80% missing and categorical accuracy stayed within about one point, compared every recorded trait and setting before correcting those claims.

## 4. Files Touched

- `DESCRIPTION`
- `NEWS.md`
- `docs/dev-log/cran-0.11-audit/site-review-2026-10-08.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-scale-claim-audit.md`

## 5. Checks Run

- Compared the 0.9.1 NEWS claim with `script/bench_scaling_v090.R`, `script/bench_scaling_v090.md`, and `script/bench_scaling_v090.log`. The report records commit `888d02052ff24513071283abb851d4ad64ebede7`, one run per size through n = 5,000, and 1,383.9 seconds total. The 5,000-tip graph stage records about 1,971 MB of R heap.
- Compared the separate AVONET sweep claim with `script/bench_avonet_missingness.md` and its log: 9,993 species, seven traits, three missingness settings, 273.9 minutes total.
- Rechecked every retained result row after independent review: all four continuous RMSE values equal BM at 80% missing; categorical accuracy differences across the six trait-setting comparisons range from −4.1 to +4.2 percentage points. The study has one run per setting and no replicate-based uncertainty estimates.
- `git diff --check`: passed after the review-driven NEWS correction.
- Fresh offline pkgdown build completed successfully.
- Site cleanup and crawler: 62 HTML pages, 3,443 local references, 34 retired routes, 0 errors; `SITE_CRAWL_OK`.
- Rendered assertions passed for removal of the blanket package-page scale claim, retention of the historical 9,993-species sweep, correct 5,000-tip scaling limit, removal of the absent script reference, and absence of the standard-laptop claim.
- `Rscript --vanilla -e 'pkgdown::check_pkgdown()'`: no problems found.
- `git diff --check`: passed.
- The after-task structure check passed. Its whole-workspace acceptance check remains red on five open `.unlazy/imputation-sim/gates/leaf-*.md` gates. Those simulation gates are outside this documentation correction, so I left them unchanged.

## 6. Tests of the Tests

Assertions inspect the generated package landing page and NEWS HTML, so stale source text or a build that does not render the correction fails the check. The site crawler separately verifies local references and all retired routes.

## 7a. Issue Ledger

- Fixed: `DESCRIPTION` claimed a current 10,000-species test ceiling without binding it to a current-version run.
- Fixed: v0.9.1 NEWS cited a nonexistent extended-scaling script and said the timing curve reached n = 10,000, while its retained v0.9.0 report stops at n = 5,000.
- Retained with scope: the separate historical 9,993-species AVONET missingness sweep.
- Fixed after review: the NEWS summary now reflects the one-run table rather than claiming continuous gains at 80% missing or categorical differences within about one point.
- Open: deployed-site refresh, BirdTree redistribution basis, exact-artifact checks, platform evidence on the frozen artifact, and independent final verdicts.

## 8. Consistency Audit

Checked candidate `DESCRIPTION`, the current 0.11.0 NEWS entry, the historical v0.9.1 scaling section, the saved v0.9.0 timing driver/report/log, the 9,993-species sweep report/log, and the freshly rendered home and NEWS pages. The scale statements now distinguish the specific data workflows and software-version evidence. After the independent review, corrected the NEWS comparison against each row of the saved sweep report; the supplemental verification below records the rebuilt site and rendered-text check. The deployed CRAN record still describes 0.10.0 and has not been changed.

## 9. What Did Not Go Smoothly

The first saved scaling driver was a local partial run, and the NEWS entry named a separate `bench_scaling_v090_extended.R` file that is not present. The actual completed output and its embedded commit were available, so no campaign rerun was needed. The CRAN page is necessarily still the prior release while this source PR remains unmerged.

## 10. Known Residuals

Chrome confirms CRAN currently serves 0.10.0. This audit does not establish performance beyond the stated historical runs or visual review of the local rendered site. The corrected NEWS is rebuilt and structurally checked, but still needs local visual review and a later deployed-site check. The wider acceptance check remains red on five open imputation-simulation gates, separate from this documentation slice. The full release remains in progress. No merge, deployment, tarball freeze, or CRAN submission occurred.

## 11. Team Learning

Memory receipt: loaded the pigauto LOAD-FIRST manifest and checked the existing release-audit record. The saved run's own commit, workload, output, and stopping size are the right basis for scaling prose; a nearby source comment or an old NEWS summary is not enough. No Golden Set regression class was in scope.

## 12. Cross-Product Coverage

Covered: package metadata, NEWS source and rendering, retained benchmark records, package landing page, site links, and retired routes. Does NOT cover: new performance campaigns, visual inspection, post-deployment site state, redistribution rights, exact-tarball checks, or CRAN acceptance.

## Supplemental verification, 2026-10-08

After the independent review, rechecked every retained result row: all four continuous RMSE values equal BM at 80% missing; categorical accuracy differences across the six trait-setting comparisons range from −4.1 to +4.2 percentage points. The study has one run per setting and no replicate-based uncertainty estimates. Built the edited candidate in an isolated source snapshot, ran `pkgdown/clean-internal-pages.R`, and checked the rendered NEWS page. It contains the corrected wording and neither superseded claim. The site crawler reported 62 HTML pages, 3,443 local references, 34 retired routes, zero errors, and `SITE_CRAWL_OK`; `pkgdown::check_pkgdown()` reported no problems. `git diff --check` passed. The live site has not been redeployed, and this check did not visually inspect rendered pages.
