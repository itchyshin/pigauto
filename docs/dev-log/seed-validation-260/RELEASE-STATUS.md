# Separate CRAN release track

This is a read-only receipt for the release owner, not a new release candidate or permission to bypass its pending gate.

The Pigauto chat owns the current 0.11.1 release. Its latest status is waiting for independent Win-builder results before final release review, the approved tag move and submission. No archive, tag, submission record or release-owner file was modified here.

## Exact artifact checked

Path: `/private/tmp/pigauto-release-0111-cross-mkl-evidence/pigauto_0.11.1.tar.gz`.

SHA-256 independently recomputed in this pass: `95998faf7cf7ee0ae5246407853aa153593626664e441bb931e7bb4e4c7ed903`.

The local full `--as-cran --run-donttest` log ends in Status: OK, including rebuilt vignettes and PDF/HTML manuals. The exact-archive GitHub Windows release and devel logs both end in Status: OK. The R-hub Ubuntu release and clang logs both end in Status: OK. All four external before/after SHA files match this artifact hash. Logs remain in the release owner's evidence directory; this receipt does not replace them.

The owner's `cran-comments-DRAFT.txt` still labels Windows and the R-hub checks PENDING. An updated wording draft is provided beside this receipt for reconciliation by that owner. It does not substitute GitHub Windows checks for the outstanding independent Win-builder results.

## Pending dependencies

- Independent Win-builder result review remains unproven in this pass. The upload log records HTTP 200 responses, not completed checks. Two targeted Gmail searches returned no matching new results; that is not proof no results exist in another account or source.
- The final rights/consent/timing/governance gate and required independent release panel are not repeated or claimed complete here.
- A concrete instruction/source inconsistency needs reconciliation in that product-contract review: AGENTS.md lines 44 and 345–346 describe covariates as ignored or GNN-only when the GNN is off. Current R/fit_pigauto.R lines 571–575 and the help explicitly permit covariates in the selected per-column numeric BM baseline. The real-data coding check in this follow-up uses that route. No shared instruction file was edited because that surface is outside this lane's lease and other lanes have live work.
- No submission receipt, maintainer confirmation or CRAN acceptance was established.
- Ayumi's latest issue #261 reply independently confirms that GitHub source 30ab1e1 fixes evaluation and the misspelt argument; she is retaining the issue to track CRAN publication. That is bounded reader verification, not release availability.

Highest supported evidence from this pass: the exact tarball has local and four external platform check receipts. The release lane's independent Win-builder and final submission gates remain open. Do not promote this to submission-ready or submitted.

Suggested replacement wording for the release owner to review: “Per-column numeric Brownian baselines can use supplied covariates. GNN-off disables neural covariate effects; it does not universally disable baseline covariates. Joint baseline paths omit them; see fit_baseline help for route-specific behaviour and warnings.” This is a draft correction to the factual description, not a proposed runtime change.
