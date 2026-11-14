The plot and simulation examples currently train the optional GNN even when its runtime is unavailable. They now demonstrate the package's baseline default. The comparison example retains its GNN comparison and first checks that torch can create a tensor; it reports a skip when Lantern is unavailable.

The installation guide also describes 0.11.0 functions while CRAN currently provides 0.10.0. README and Getting Started now check the installed version and explain when to use the GitHub version. Generated help and the vignette R script agree with their sources. No scientific default or runtime implementation changes.

This source-only followup extracts the example-portability changes already checked in PR236 and adds the installation guidance. PR236 remains the unmerged release-evidence PR.

Validation: R parsing and Rd checks pass; the source community-surface tests pass; the fresh local site build and crawl pass with 62 HTML pages, 3,587 local references and zero crawl errors. Independent source review is recorded in PR236's preparation receipt. The preceding c2 archive passed its exact-artifact checks, but the changed README and vignette require a new final archive after merge. Fresh CI on this PR is required before merge.

Merging will trigger the main-push website deployment. Deployment verification, clean merged-source freezing, final artifact checks, timing adjudication and the final release panel remain open. No CRAN submission is requested by this PR.
