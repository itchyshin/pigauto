# Final pigauto 0.11.0 artifact evidence

All four native checks and the final installed-model matrix have completed on
archive SHA-256
`14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb`.
The archive is 5,137,648 bytes with 250 entries, generated from clean merged
source `8e275ca738129c8d608f06ffc407b8141c1a2636`. The identical read-only
retained copy and its checked-path binding are in
`persistent-archive-retention.json`.

The current decision record is `closure-summary.json` together with the audit's
current `GATES.md` and `release-ledger.json`. Native macOS has one timing NOTE;
final independent timing disposition and the maintainer's submission decision
remain open. No CRAN submission or confirmation has occurred. Source PR #237 is
merged and deployed with explicit approval. Evidence PR #236 remains unmerged.

## Exact checks and their limits

Direct workflow 39082452558 checked the original gzip on Windows release,
Windows development and Linux release. All three report Status: OK with native
exit 0, all 29 declared dependencies available, NOT_CRAN unset, force-Suggests
enabled and full manuals. The retained platform bundles contain identity
records before and after execution, actual logs, environment/session details,
manual outputs and separate warning-capture data. Native totals were 645.813,
622.432 and 336.330 seconds respectively. Windows tests took 415/404 seconds;
vignette rebuilds took 22/22 seconds. Linux's brackets report rounded CPU/elapsed
pairs, including 206/206 for tests and 11/11 for vignettes.

The Mac native check exited0 with Status: 1 NOTE. Its strict driver returned 1 to
require explicit adjudication. The sole NOTE records compare_methods at
2.532 seconds user,0.992 system and 18.406 elapsed. Native total 3564.464 seconds
and driver total 3687.675 seconds are separate observations. Tests were logged
at 5 minutes CPU/31 minutes elapsed, with vignette rebuilds 17 seconds CPU/120 elapsed seconds.
No native run was repeated to select a faster result. The separate warning
capture completed in 687.196 driver seconds and 669.225 reporter seconds.

All four native suites report 3195 passing expectations,0 failures,161 warnings
and 86 skips. The four-platform warning audit passed using retained objects from
separate replays on the same checked installations. Messages, source locations,
callers and skip identities match the reviewed fixture contract; totals match
the native reporter. The original native reporter retained totals rather than
individual warning objects. The audit does not establish general coverage or
inferential validity. The largest group is 119 small calibration sample warnings;
other categories and their consequences are explicit in the actual result.

The first actual audit stopped before object review because its Python reader
expected JSON null while jsonlite encoded successful R NULL as an empty object.
Both the failure and successful fresh result are retained. The narrow reader
fix accepts only null or an empty object, requires the error field, and retains
all other environment/identity/provenance/count/problem checks. Independent
review exercised 142 cases with 136 refusals and checked 138 protected files.
Raw receipts and the comparison contract were unchanged.

## Installed workflows and ancillary proofs

The final optional-model matrix passed all 12 steps across neither package,
drmTMB only, gllvmTMB only and both. Both measured backend versions were 0.7.1.
Saved/reloaded objects, fresh unattached sessions, extraction/Rubin arithmetic,
provenance refusals and fit health are recorded. Fresh pigauto libraries used
recorded dependency symlinks. These are bounded three-model integration fixtures;
the separate n = 300, M = 20, ten-seed Gaussian recovery evidence remains in
`recovery/RESULTS.md` within its stated zero true phylogenetic covariance regime.
Neither backend is a pigauto installation dependency. Pooling supports termwise
fixed effects within the documented proper-imputation contract.

Missing-rgbif/terra fixtures passed 60 expectations with0 failures,0 warnings and
8 expected skips. They establish the recorded absence/cache/failure behaviour,
without proving live external-service availability. All 6 individual vignette
diagnostics passed with 74.8865 seconds render time,308.2825 summed R-process
seconds and 310.1696 seconds for the enclosing loop. Disabled teaching chunks
remained disabled. Native installed-size INFO was 5.9Mb remotely and 7.0Mb
locally; these are platform observations on unchanged archive bytes.

The URL checker retained one distinct DOI 403 in the Getting Started source and
HTML. Primary Crossref returned 200 and matched the Hackett 2008 Science citation.
The identifier is supported; automated DOI access remains denied. Review of
404 spelling flags in 1886 contexts across 43 archive-bound files identified no
ordinary-word typo among those flags. This does not certify unflagged prose or
all references. Zero direct reverse dependencies were found in the retained
primary CRAN-index refresh; the package has no owned compiled source.

## Current records and preserved snapshots

`warning-audit/actual-results-null-error-fix/` is the current successful warning
result. `warning-audit/actual-results-with-native-note/` is the preserved first
failed audit. `replay-serialization-fix/` records reproduction, source comparison,
controls and independent review. `external-platforms/` and
`checks-and-integration/` preserve the remote, service, spelling, vignette and
installed-model inputs with member hashes and standalone readback helpers.

`local-capture-snapshot/` preserves raw Mac capture proof and the packager's
earlier NOT_DONE statement about the first warning audit. That statement is a
historical snapshot superseded by the successful audit above. The native NOTE
and original strict-driver failure remain unchanged. Likewise,
`checks-synopsis-snapshot/` retains its capture-complete/consequences-pending
snapshot. Neither snapshot overrides the current decision records.

The earlier derived Linux timing parser and an ambiguous vignette process-total
label were corrected against raw evidence. A copied reproduction instruction
was corrected to name its actual helper and original scratch-input requirement;
prior versions and raw bundles remain preserved. The bundles support independent
byte/hash readback; regeneration of the original proof runs requires their
recorded inputs and software, including temporary paths where stated.

Actual deployment 8e275ca passed 44 retired and 62 retained routes, sitemap 62,
search 614 and bounded six-page own-browser review. Source CI's MPS-related step
does not prove GPU hardware execution; the separate local torch probe was CPU-only.
Rights and consent retain the recorded component ledger, sole-author consent
and maintainer tree-use warranty with credit. No broader BirdTree licence is
asserted. This evidence does NOT cover joint covariance/loading/variance-component
pooling, general interval coverage, untested missingness/families, live service
availability, GPU hardware validity or CRAN acceptance.

Assessment of this current wrapper remains pending exact-draft validation.
