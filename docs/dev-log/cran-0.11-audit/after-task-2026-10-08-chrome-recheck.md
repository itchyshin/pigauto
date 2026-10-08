## 1. Goal

Refresh the Chrome audit of deployed pigauto tree pages, compare it with merged main and the latest successful Pages build, and keep release gates tied to verified evidence.

## 2. Implemented

Added a dated Chrome site audit and appended current deployment evidence to G7. The record corrects an interim chat misread: the live `tree300` help remains stale.

## 3a. Decisions and Rejected Alternatives

Kept G4, G6, and G7 open. The successful pkgdown run binds the observed deployment to merged main but does not close the full route and sitemap checks or validate pending PR #231 changes. No source edits, merge, deployment, or publication action was taken.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/site-review-chrome-2026-10-08.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-chrome-recheck.md`

## 5. Checks Run

- Chrome opened five cache-busted public pages: `tree300`, `trees300`, `tree_full`, `read_tree`, and Getting Started.
- Chrome opened GitHub `main` source `R/data.R` and generated help `man/tree300.Rd`; both retain the older `tree300` wording.
- Chrome opened pkgdown Actions run #583. It completed successfully on `main` commit `0b0f71fee838c6ed51ef832ed819270eeafaf29b` and links to the public site.
- Chrome opened the official CRAN record, which lists 0.10.0 published 2026-07-30.
- PR #231 remains a draft, has no review, and reports no deployment for its branch.
- The shared lane lease was granted for only the three listed evidence files. `git diff --cached --check` passed after staging those paths; whitespace checks passed for both new files; the after-task structure check passed; and `slop_check.py` reported zero findings.
- `node ~/shinichi-brain/skills/unlazy/scripts/gate-check.mjs --status docs/dev-log/cran-0.11-audit/GATES.md` parsed 11 gates and reported six open: G0, G4, G6, G7, G8, and G9. This status-only command reads checkboxes and runs no gate checks.
- `closeout.py check` did not pass: its nested R validator reads the brain repository's unrelated acceptance ledgers and halts on their unmet gates. Running `check-after-task.R` directly from this pigauto worktree passed the report structure check; the pigauto release ledger remains separately open as recorded above.

## 6. Tests of the Tests

The source comparison prevents a live help page from being called corrected based only on a prior note or candidate branch. Each page was checked in the rendered browser, and the source link was opened on GitHub `main`.

## 7a. Issue Ledger

- Confirmed: all five checked reader pages retain at least one provenance, rights, or retrieval-guidance mismatch.
- Corrected: the interim assertion that live `tree300` was already corrected.
- Open: reconcile current source/manual/help with the provenance audit; user-directed retrieval and citation guidance remain needed on the novice path.
- Open: full live retained/retired route, search, and sitemap verification; visual review of the local candidate; data rights; final artifact and platform checks.

## 8. Consistency Audit

Compared the five live reader pages with `main` `R/data.R` and `man/tree300.Rd`, checked the public version badge against CRAN's current version, and checked the latest successful pkgdown workflow record. The live `tree300` page agrees with the older committed source; other tree pages remain inconsistent with candidate provenance guidance. The site badge is not CRAN publication evidence.

## 9. What Did Not Go Smoothly

I initially misread the live `tree300` wording as corrected. Rechecking the returned browser text showed it still describes a Hackett sample under the `megatrees` MIT licence. The corrected finding is now recorded. The shared lease registry was inaccessible from the ordinary sandbox; a narrow lease was granted through the approval path before repository writes. The global closeout wrapper reads the brain repo's open ledgers when called from this worktree, so I recorded its failure and ran the structure validator in the pigauto context directly.

## 10. Known Residuals

No product or source correction was made. This audit does not establish BirdTree redistribution permission, verify every live URL or sitemap entry, inspect the local rendered candidate, freeze or check a final tarball, or provide checksum-bound Windows evidence. G4, G6, G7, G8, and G9 remain open.

## 11. Team Learning

Compare current rendered text and its linked `main` source directly. A summary or previous browser observation can be stale or misread, and a site badge does not establish the CRAN version.

Memory receipt: loaded the pigauto route manifest, lane preflight, ultra-plan, and unlazy guidance; the release-ladder and ownership rules shaped this slice.

Golden Set: not applicable because this slice audited rendered documentation and deployment state.

## 12. Cross-Product Coverage

Covers: five live reader pages, the current `tree300` source/manual, CRAN's package record, and the latest successful pkgdown workflow.

Does NOT cover: the full website inventory, sitemap and search closure, local candidate visual review, source correction, redistribution rights, final artifact validation, platform results, or CRAN submission.
