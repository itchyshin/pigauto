# Open-issue verification and diagnostic closeout

## 1. Goal
Execute the approved issue plan on isolated source, verify existing fixes and keep missing-data, validation and release work separate.

## 2. Implemented
Added two end-to-end R regression tests for issues #258 and #259. Added a descriptive audit script and summary for the public seed-prediction attachment in #260. Recorded each issue's remaining work in README.md. No new package runtime defect was demonstrated, so no runtime behaviour changed.

## 3a. Decisions and Rejected Alternatives
Retained existing fixes rather than reimplementing them. Did not add seed averaging, prediction jitter, interval inflation or new defaults without supporting validation. Did not close issues awaiting reporters. Assumptions: the nine open issues form the scope; cached public predictions can check arithmetic but cannot validate current source. No raw private data or public CSV rows were committed.

## 4. Files Touched
All paths below are relative to this worktree:

- tests/testthat/test-open-issues-verification.R
- docs/dev-log/open-issues-20261108/README.md
- docs/dev-log/open-issues-20261108/summarize_seed_predictions.py
- docs/dev-log/open-issues-20261108/seed-summary.json
- docs/dev-log/open-issues-20261108/after-task.md
- docs/dev-log/open-issues-20261108/issue-replies-draft.md

The original dirty checkout was untouched. The user subsequently authorized committing, pushing and opening a PR, with draft replies kept separate from posting. Changes are on codex/open-issues-nov08. No merge or issue reply is authorized by that follow-up.

## 5. Checks Run
On base 30ab1e1b50b34a732485f6cccb195e9fa75299fc, R 4.6.0:

- Existing issue-input-clarity, issue-visuals, ordinal-levels, lambda-covariates and report-valid-denominators tests: 147 passing assertions, zero failures or warnings, one skipped 600-replicate recovery test under NOT_CRAN=false.
- New public-path tests: 12 passing assertions, zero failures, warnings or skips.
- Final affected evaluation-mode, check-pigauto, read-traits, ordinal-levels and public-path tests: 179 passing assertions, zero failures, warnings or skips. Counts overlap earlier runs and must not be added.
- Python summary reproduces 375 cells, 77 species, median seed range 0.85234 degrees C and 30 ranges above 3 degrees C.
- CRAN package page read on 8 November: version 0.11.0, so #261 retains its release dependency.

No exports or roxygen changed; no documentation regeneration needed. After the user authorized commit and PR creation, a staged-source export passed devtools::check with zero errors, warnings or notes in 14 minutes 30 seconds. That run included the package test suite and examples, including donttest examples; manual generation and vignette rebuilding were excluded. The initial dirty-directory check was stopped after it included untracked Graft caches; the clean export excluded those files. No site rebuild or release submission occurred.

The full testthat summary was 3,726 passing assertions, 34 skips and 184 warning conditions, with zero failures. These warning conditions are distinct from R CMD check's zero-warning status; this follow-up did not triage every full-suite warning. The new test file and final affected test run had zero warning conditions.

## 6. Tests of the Tests
The Python audit rejects a missing seed-cell, nonfinite prediction and inconsistent truth. A fresh read-only reviewer found no actionable issues in the R tests or Python summary and independently reproduced the attachment summary. A separate process temporarily removed independent-column selection: three assertions in the factor-coding regression failed while the retained-row test passed. The process restored the original function before exiting. No package file was mutated.

## 7a. Issue Ledger
#254: merged documentation represented and relevant regressions pass; reporter acceptance and slow-install cause remain open. #257: plot/report regressions pass; requested visual scope still awaits reporter feedback. #258: retained rows and warning verified. #259: public factor/numeric coding invariance verified in the reporter DGP. #260: cached arithmetic reproduced; new-source stability and remedies unvalidated. #250: tie reproducer still needed. #251: stored draws/model specifications still needed. #261: release outstanding. #135: broad comparative campaign separate. Issue states unchanged.

## 8. Consistency Audit
Checked current vignette sources for the CRAN evaluation workaround, posterior fields and repeated-observation/seed limitations. Checked report and ordinal paths alongside the two public regressions. No claim that every test of a function establishes scientific validity. The seed attachment comes from an older source than the tested package checkout.

## 9. What Did Not Go Smoothly
The local graph indexes the old dirty branch; source checks used origin/main and the isolated checkout. docs/ is ignored, so review artifacts require explicit staging. The first negative-control test harness tried to add helper bindings to a locked namespace; rerunning with a child test environment exercised the intended fault. A previous closeout wrapper resolved a different root; this report uses the absolute validator path and worktree.

## 10. Known Residuals
Private mammal inputs were not obtained. No new GNN fits, seed-remedy validation, large campaign, CRAN submission, public issue reply or issue closure occurred. The exact-draft prose assessment compiler was not used: NOTASSESSED. All scientific questions retain the limits in README.md.

## 11. Team Learning
Memory receipt: canonical contract already loaded once; route.py pigauto manifest, efficiency, ultra-plan, systematic-debugging, verification-before-completion, r-package-engineer, validation-harness and quantitative-analysis informed the scope. Prior-work retrieval found no matching brain note. Golden Set: not run; no implementation changed. Reused the existing fixes and requested missing evidence instead of duplicating work. A clean scoped lease protected these paths; no brain memory was edited. Reconciliation: planned runtime repairs were unnecessary because relevant existing fixes passed; the delivered additions strengthen regression coverage and clarify pending evidence.

## 12. Cross-Product Coverage
This slice covers public baseline-only factor coding, retained rows through memory/CSV input, existing ordinal/report/evaluation checks, the package suite as exercised by the recorded check and descriptive cached seed summaries. It does NOT cover GNN runtime stability validation, private mammal model agreement or interval calibration, a new all-trait scientific validation campaign, manual/vignette rebuilds, live-site rebuilds or a release. Data-dependent, validation-dependent and release-dependent work remains separately listed.

Review note: self-review of the whole report and README on 8 November 2026. Perceived AI-like style 3/10, moderate confidence, with mandated headings and a dense issue ledger making the report more formal than ordinary prose. The source SHA and separate track descriptions help the reader bound the findings; no wording repair needed. Science gate: pass for descriptive claims with explicit limits. Facts gate: pass against observed checks and the negative-control output. References gate: pass for issue and CRAN URLs in README. Formal exact-draft assessment remains NOTASSESSED.
