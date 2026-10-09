# Exact-head CI and PR status refresh: after-task report

## 1. Goal

Record the completed exact-head candidate-source matrix for PR #231 and reconcile the visible status of both release PRs.

## 2. Implemented

Added the successful run #37888439948 receipt for source head `a48f61c748f72e5351c1a871c05cc52439e61858` to the audit ledger. Updated PR #231 to replace the pending run state with the completed result and refreshed PR #228 to show the current source head, completed run, local site review, and 8-of-11 gate count.

## 3a. Decisions and Rejected Alternatives

Candidate-source CI proves the three listed matrix jobs completed on the recorded source head. It does not prove the final artifact, force-Suggests checks, deployment, or independent review. G7, G8, and G9 remain open. No scientific default, package behavior, or website source changed in this receipt slice.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-ci-and-pr-refresh.md`

External review records updated in Chrome: PR #231 and PR #228 descriptions.

## 5. Checks Run

- Chrome Actions run #37888439948: Success, 3 of 3 jobs completed in 18m10s. Ubuntu R release, Ubuntu R-devel, and macOS arm64 R release succeeded.
- Chrome job pages: Ubuntu R release succeeded in 7m23s; Ubuntu R-devel succeeded in 12m45s; macOS arm64 R release succeeded in 18m02s.
- Chrome PR #231: description now identifies #37888439948 as completed on head `a48f61c`; PR remains Draft and unmerged.
- Chrome PR #228: description now reflects the completed source run and 8/11 audit tally; PR remains Draft and unmerged.
- Unlazy ledger status: 11 gates, 8 met and G7–G9 unmet.
- No code or website tests were rerun because this slice changes only release evidence.
- `Rscript ~/shinichi-brain/tools/check-after-task.R <report>`: the 12-section structure passed; the command exits 1 because the repo-wide Unlazy scan finds five unmet gates under `.unlazy/imputation-sim/gates/`.
- `slop_check.py <report>`: 0 findings.
- The five imputation-sim gates are outside this CRAN ledger slice and were left unchanged. They prevent whole-repository closeout, but do not change the separate CRAN ledger count of 8/11.

## 6. Tests of the Tests

GitHub's run summary confirms all three matrix jobs succeeded against the recorded PR head. This supports candidate-source CI only. It does not test a frozen tarball or the deployed website.

Golden Set: No implementation behavior changed in this evidence-recording slice.

## 7a. Issue Ledger

Resolved: the release PR descriptions no longer describe the successful #37888439948 matrix as pending, and the local ledger contains its exact source head and result.

Open: G7 deployed-site verification; G8 one exact post-merge tarball with local and platform results; G9 independent review of the exact artifact and deployed-site evidence.

## 8. Consistency Audit

Both PR descriptions state that the source changes remain unmerged and the website undeployed. The current source matrix passed, but the release tally remains 8/11 because its three outstanding gates require post-merge and exact-artifact evidence.

## 9. What Did Not Go Smoothly

GitHub CLI could not connect to `api.github.com`; Chrome supplied the authoritative run and PR state. Chrome editor updates completed and were verified in the rendered descriptions.

## 10. Known Residuals

This receipt does not satisfy deployed-site verification, exact post-merge artifact checks, Windows or other platform checks bound to that artifact, or independent final review. The repo-wide closeout validator also remains red on five imputation-sim gate files; this slice did not change them. Merge, deployment, and submission remain maintainer-controlled.

## 11. Team Learning

For each candidate-source matrix, record the exact head and completion duration. Keep source-CI results visibly separate from exact-artifact and deployment gates.

## 12. Cross-Product Coverage

This slice covers one completed three-platform candidate-source matrix, two PR description refreshes, and their local audit receipt. It does not cover merge, deployment, an exact final tarball, platform checks for that tarball, independent artifact review, or CRAN submission.

## Follow-up, 2026-10-09

Chrome verified run [#37897100223](https://github.com/itchyshin/pigauto/actions/runs/37897100223) on PR #231 head `8a00ce59382c7e7ba22845425299127d861d66e2`: the run succeeded in 17m48s and its Ubuntu R release, Ubuntu R-devel, and macOS arm64 R release jobs all completed successfully. The workflow sets `_R_CHECK_FORCE_SUGGESTS_=false`; this remains candidate-source evidence, not the exact-artifact force-Suggests gate. The run-page title identifies head `8a00ce5`, and the PR descriptions already reflect this receipt. The ledger now records it while preserving the 8-of-11 tally and open G7–G9 gates. No package source, user-facing documentation, bundled data, or website input changed in this receipt update.

Chrome verified run [#37899632824](https://github.com/itchyshin/pigauto/actions/runs/37899632824) on PR #231 head `2ebd5b97ab2a3e9b3c4c7507d568d365b478532c`: the run succeeded in 18m22s and Ubuntu R release, Ubuntu R-devel, and macOS arm64 R release all passed. The workflow sets `_R_CHECK_FORCE_SUGGESTS_=false`; this remains candidate-source evidence, not the exact-artifact force-Suggests gate. This head contains only the audit ledger and after-task records relative to the G0-reviewed source boundary. The ledger preserves the 8-of-11 tally and open G7–G9 gates. No package source, user-facing documentation, bundled data, or website input changed.

## Audit correction, 2026-10-09

The source ledger had mislabeled the October 7 predecessor archive as a post-merge tarball and described its older source/docs state as uncommitted. The archive still exists; its SHA-256 was rechecked as `00a0a323b46cb1d9a5d49c3542635ca529ad83ae95216db1ac38c3362913bf6e`, its size is 5,127,386 bytes, and its timestamp is October 7. It predates unmerged PR #231 and is now explicitly predecessor evidence. Its earlier macOS check remains a predecessor result; the CRAN-shaped incoming-feasibility step stopped on unavailable package-index access, and Win-builder results were still pending in the last recorded inspection. The final artifact gate remains open.

The same edit aligns G4's rights wording with G0: the ledger records Shinichi's maintainer warranty basis while keeping clear that the public BirdTree pages do not state a separate redistribution licence. The fail-closed release-ledger command returns `NOT READY`; product-contract, rights-and-consent, current-policy, and rendered-site evidence fields lack typed entries, and the artifact object is empty. `git diff --check` passed. Unlazy reports 8 of 11 release gates met, with G7–G9 open. The after-task structure check passed; its process exits 1 because the repo-wide scan also finds five unmet gates under `.unlazy/imputation-sim/gates/`, outside this CRAN slice. The naturalness check found zero findings. Independent ledger review found no mismatch in the correction and confirmed the eight-path G0 scope. No package or site input changed.

## Final site crawler and visual review, 2026-10-09

The site reviewer found that the original crawler missed responsive `srcset` and CSS assets. The first reruns surfaced separate build issues: the exact-HEAD script intentionally used committed files, so an uncommitted `--base-url` change was absent; after that change was committed, the crawler found a broken `/NA` script URL on the generated 404 page. Adding `.github/404.md` alone did not fix it because pkgdown rewrote the inline script's missing `src`; moving the script into `pkgdown/extra.js` fixed the rendered page. A later reviewer pass found that HTML symlinks could be read before root containment was checked. A red regression test reproduced that path, and the guard now rejects external-root HTML symlinks before reading page content while retaining in-root alias behavior.

At exact source commit `86efab58f1d36e1e133be8137fa5d44a545182cb`, the fresh build completed into `/private/tmp/pigauto-cran011-site-86efab5.gAlMYY`. Six crawler tests passed. The crawl checked 62 HTML pages and 3,584 HTML/CSS references; all 34 retired routes were absent from output, search, and sitemap, with zero errors and `SITE_CRAWL_OK`. The build wrapper exited successfully after `pkgdown::check_pkgdown()`; its retained build log does not include that command's own output. Build-log SHA-256: `0469f9e5f056183aec18ccaa85bff5fd86b2e699cb99571f355ea3531252e559`; crawl-log SHA-256: `666b188ae633c31388a96b76c6b704ae575a69281e6606c220539e00f8f159d4`.

Chrome visually reviewed the custom 404 page, Getting Started, the multiple-imputation article, and `trees300` at the immediately preceding `437e6e0` build. The four reviewed HTML files are byte-identical in the final `86efab5` output; hashes are recorded under G6. The final independent reviewer verified the 34 retired routes and four simulation-study variants are absent from output, search, and sitemap, confirmed the exact receipts and source commit, and repeated both external-root rejection and in-root symlink acceptance. CSS parsing remains regex-based; the static crawler does not cover JavaScript-generated links, deployed behavior, or resolution of every sitemap entry. G7 deployment verification and the artifact gates remain open. No merge, deployment, or CRAN submission occurred.
