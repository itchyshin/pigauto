# Approvals before a pigauto 0.11.0 submission

This record lists the approvals still needed before submission.

## Approval needed next: source merge and website deployment

Source: [PR #237](https://github.com/itchyshin/pigauto/pull/237), commit `259d7e98c86dd5aa7680367d9e621ad76037d323`. Evidence: [PR #236](https://github.com/itchyshin/pigauto/pull/236), kept unmerged.

The approved release plan reserves merges to Shinichi and leaves the evidence PR unmerged. Approve merging the separate tested source-fix PR into main when all applicable checks are settled green. PR #236 stays open and unmerged as the release-evidence PR. Merging PR #237 will trigger the automatic Pages deployment. A normal merge preserves the source and evidence history.

The added installation guidance distinguishes the CRAN release from the GitHub version. Primary CRAN metadata currently lists 0.10.0, while the candidate documentation describes 0.11.0 and functions introduced in that version. The existing c2 archive is predecessor evidence for the updated documentation.

After an approved merge, Codex can verify the live deployment, rebuild from clean merged source, freeze the final archive and rerun affected exact-artifact checks under the existing implementation/check authority. Estimated remaining time: 60-100 minutes, with deployment and platform checks run concurrently where independent. No new simulation campaign is proposed.

## Decision needed after final timing results

The previous Windows update checks completed with Status OK in 536.636 seconds on R-release and 629.689 seconds on R-devel. The devel test stage took 407 seconds and vignette rebuilds took 21 seconds. An earlier successful devel run took 656.047 seconds. CRAN will assess the final submission independently.

The final candidate will have its own timings and warning review. Shinichi must decide whether to proceed with those recorded update-check costs or request a scoped runtime reduction. A timing decision does not waive check failures, rights, identity or reviewer blockers.

## Separate final approval: CRAN upload and maintainer confirmation

Only after the deployed-site check, final frozen-artifact checks, fresh Grace/Rose/Pat reviews and submission-ready validator pass, Codex will present the exact archive path, SHA-256, size and final cran-comments text. Approve that specific upload and its maintainer confirmation then. Merge approval alone does not authorize it.

No new author or tree-rights approval is requested: Shinichi is the sole credited author/maintainer, and the recorded tree-use warranty with citation and credit remains the stated rights basis. No separate BirdTree licence is asserted. Public tagging and release announcements are outside this submission-preparation scope.
