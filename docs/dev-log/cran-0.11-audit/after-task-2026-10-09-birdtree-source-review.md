# After-task: independent BirdTree source review

## 1. Goal
Close the two previously raised source-documentation findings against PR #231's current candidate without changing release gate status.

## 2. Implemented
Recorded the independent review and its focused test receipt in `GATES.md`; retained the raw test output with its SHA-256.

## 3a. Decisions and Rejected Alternatives
The candidate does not claim that member 69 is closest or best. Users control BirdTree download and cite the source. This source review does not prove deployed-site state or tree-data redistribution rights.

## 4. Files Touched
- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-birdtree-source-review.md`
- `docs/dev-log/cran-0.11-audit/provenance/readtree-source-review-dee8f028-2026-10-09.log`

## 5. Checks Run
At exact source head `dee8f02888fd2c4c22d426e595ea6451a28a280e`, `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla -e 'devtools::test(filter = "shipping-coverage|read-tree-multipho", stop_on_failure = TRUE)'` passed 52 testthat expectations across two test contexts, with 0 failures, 5 warnings, and 0 skips. Wall time was 8.1 seconds. After push, Chrome confirmed PR #228 includes commit `5ad0bf9`; Actions run [#37960057428](https://github.com/itchyshin/pigauto/actions/runs/37960057428) was skipped in 8 seconds under the PR pkgdown guard with no job or artifact. The after-task structure checker passed; the aggregate acceptance check remains red on five unrelated imputation-sim leaf gates (`leaf-campaign`, `leaf-env`, `leaf-prerun`, `leaf-results`, `leaf-runner`), which this receipt does not own. A whitespace check excluding the retained raw log passes; the byte-preserved test log has trailing spaces on live progress-bar lines, so a whole-diff whitespace check flags those log lines. The retained combined output is 5,406 bytes with SHA-256 `0ec6dc71c01533461dc43e4895e7d81f34d030a23d32cffb3f6f750cbceaa356`.

## 6. Tests of the Tests
The test filter exercises Newick/NEXUS reading and multiphylo handling; the reviewer confirmed these are focused behavior tests, not a proof of all tree workflows.

## 7a. Issue Ledger
Both prior source findings are resolved at the reviewed candidate head. G7 deployed-site verification, G8 frozen-artifact checks, and G9 independent final-artifact review remain open.

## 8. Consistency Audit
`inst/NOTICE`, `R/data.R`, `README.md`, `R/read_tree.R`, the Getting Started and tree-uncertainty articles, and generated help agree on example-tree status, user-run acquisition, formats, and citations. The historical Oct 8 finding is retained as a dated record.

## 9. What Did Not Go Smoothly
Chrome's GitHub tab timed out during a separate refresh attempt; the source reviewer completed the pinned-head review and test run from the exact source worktree. The raw test log contains padded live progress-bar output; it is retained byte-for-byte with its hash. No new live-site claim is made. PR #228 remains Draft and unmerged, and the site has not been deployed from this branch.

## 10. Known Residuals
The source PR remains unmerged. The deployed website still requires post-deployment route, search, and sitemap checks. The exact final tarball and platform checks do not yet exist.

## 11. Team Learning
Keep dated audit findings intact, then add a head-bound review receipt when a source change resolves them.

## 12. Cross-Product Coverage
This receipt covers source documentation and focused tests only. It does not cover deployment, the final artifact, CRAN submission, or clearance beyond the maintainer's recorded tree-use confirmation.
