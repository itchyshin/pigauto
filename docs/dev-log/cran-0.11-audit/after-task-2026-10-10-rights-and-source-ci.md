## 1. Goal

Record the fresh BirdTree rights-source recheck and PR #233 candidate-source CI without overstating CRAN readiness.

## 2. Implemented

Added a dated rights-source supplement and recorded the PR #233 workflow result in the release gate ledger. The result remains NOT READY, with exact-artifact platform evidence and independent closure open.

## 3a. Decisions and Rejected Alternatives

Kept the original rights report intact and added a dated supplement. The current megatrees MIT metadata is not treated as a licence for the underlying BirdTree data. Candidate-source CI is not treated as frozen-artifact evidence. Assumption: the previously recorded maintainer warranty remains the disposition; if Shinichi revises that warranty, the rights decision must be reopened.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/rights-recheck-2026-10-10.md`
- `.unlazy/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/release-ledger.json`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-10-rights-and-source-ci.md`

## 5. Checks Run

- `python3 ~/shinichi-brain/tools/slop_check.py "$PWD/docs/dev-log/cran-0.11-audit/provenance/rights-recheck-2026-10-10.md"`: passed, 269 words, 0 findings.
- `git diff --check`: passed.
- Reviewed current `rights-and-policy.md` and the release-gate ledger before adding the supplement.
- Current PR #233 workflow run #38067484909 was recorded from the prior Chrome review: three source-CI jobs succeeded and pkgdown was skipped by its configured pull-request guard.
- Unlazy reverify passed G1 (3 verifier tests), G2 (release-ledger JSON parse), G6 (16 verifier tests), and G7 (`CANDIDATE_SITE_OK`). It remains exit 1 because G3 and G4 are unmet.
- Chrome rechecked PR #228 at 68 commits, head `2921f94d1227a08ebcdbb960c386db338e914c05`, Draft/unmerged; run #38070572404 was skipped under the pull-request pkgdown guard. The JSON ledger now matches that snapshot and parses successfully.
- `python3 ~/shinichi-brain/tools/cran_release_gate.py docs/dev-log/cran-0.11-audit/release-ledger.json`: `READY FOR CLAIMED RUNG` at `source-clean`; next unproven rung is `tarball-clean`.
- `Rscript ~/shinichi-brain/tools/check-after-task.R ...`: report structure passed, but overall completion was withheld because the CRAN audit ledger has unmet gates and five unrelated `imputation-sim` ledgers remain open.

## 6. Tests of the Tests

No package code or user-facing behavior changed. G6's verifier includes 35 negative-control scenarios. The release validator was not rerun; G2 verifies only that the release-ledger JSON parses.

## 7a. Issue Ledger

- Preserved: independent BirdTree data redistribution terms are not published in the pages checked. The recorded maintainer warranty remains the basis; no independent legal determination is made.
- Open: inspect tree objects and notices in the eventual frozen artifact; G8 and G9 remain open.

## 8. Consistency Audit

Checked the new wording against the existing rights report, CRAN policy row, frozen-artifact caveat, and gate status. The evidence distinguishes BirdTree data from megatrees package metadata and source-CI from exact-artifact checks.

## 9. What Did Not Go Smoothly

The graft graph refresh failed with an operating-system permission error on its cache lock. The shared lane registry and approval records are outside the sandbox's writable paths; their scoped updates required approved escalated calls. The repository ignores `docs/` broadly, so new evidence files require explicit staging. The first gate rerun found two stale working-directory assumptions in the ledger; those commands were corrected and all four runnable local gates then passed.

## 10. Known Residuals

The tree redistribution basis is not independently established by the sources checked. PR #233's CI validates candidate source only. The existing tarball remains a predecessor artifact; no platform result is newly bound to its hash. G3 and G4 in the active Unlazy ledger remain open, so the full release-audit closeout is incomplete. The after-task verifier also reports five open `imputation-sim` ledgers outside this audit's ownership; they were left untouched. The original slice did not merge or deploy anything; since then PR #232 is merged and Pages has deployed `4f60a18`. No CRAN submission occurred.

## 11. Team Learning

Keep package-software licence metadata separate from upstream data rights, and bind each platform result to the artifact before moving the release rung. Memory receipt: the pigauto LOAD-FIRST manifest and exact-artifact rule were loaded and applied. Golden Set: not run because no code or user-facing behavior changed.

## 12. Cross-Product Coverage

Covers: rights-source wording in the release ledger and the candidate-source CI receipt.

Does NOT cover: frozen-artifact contents, Windows/platform results, deployed-site state, or independent exact-artifact verdicts.

## 13. Follow-up: deployed-site recheck after PR #232, 2026-10-10

Pat's independent live-site review was returned after the original report. PR #232 is merged as `4f60a18ae72e836b6c3e15b65abb1faa7409e731`; Pages run #38063988962 succeeded. Pat reports 62 rendered sitemap paths, 34 retired routes plus `/VALIDATION_LEDGER.html` showing the pkgdown not-found page, the corrected homepage callout, opt-in GPU wording on Getting Started, and no result for a live search of the retired ledger. The exact-source build/crawl counted 62 HTML pages, 3,584 references, and 34 retired pages absent from output, sitemap, and search. Direct live sitemap/search JSON access was blocked and Chrome did not expose numeric HTTP status codes. This closes G7 for the checked deployment only. Full current evidence and its limits are recorded in the new GATES entry above.

The source diff from frozen commit `d76804e` to `4f60a18` is limited, within the checked package/source paths, to `README.md` and `vignettes/getting-started.Rmd`; no `R/`, `tests/`, `DESCRIPTION`, `NAMESPACE`, or data changes occurred. The frozen archive was built from `d76804e` and is stale for current source. After updating the ledger's highest rung, `python3 ~/shinichi-brain/tools/cran_release_gate.py docs/dev-log/cran-0.11-audit/release-ledger.json` returned `READY FOR CLAIMED RUNG` at `source-clean`; the next unproven rung is `tarball-clean`. PR #233 remains open. A fresh archive, checksum-bound platform checks, and independent exact-artifact reviews remain open; the verdict remains NOT READY.
