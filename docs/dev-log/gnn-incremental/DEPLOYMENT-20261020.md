# GNN incremental-value push and deployment record

## Authorization and scope

On 2026-10-20, after the campaign's local closeout, the user explicitly authorized pushing the committed branch and deploying the pkgdown update. The repository is public. Its Pages workflow deploys from `main`; the reviewed branch will go through a pull request, and merge will wait for all required checks to settle successfully. This follow-up adds no simulation, public API, method, or default changes.

## Starting evidence

The Totoro campaign completed at 05:20 MDT with 2,500/2,500 datasets, 5,000/5,000 pigauto gate fits, and 30,000 successful method attempts. The strict audit and summary regeneration passed, with no failures, invalid outputs, skips, or nonconvergence. The workers exited. The scoped campaign ledger was independently reverified with all eight gates met, including a successful pkgdown build and structural check using a 600-second timeout.

The feature branch was fetched against current GitHub `main` at `b8153d8d44c40361ad87ee0770e4860e5a0e3860`. At this checkpoint the branch had 30 commits ahead and zero behind; no remote branch existed yet. Deployment gates D1-D4 remain open until the branch push, merge, Pages run, and live-page check are complete.

## Release outcome

Pending. The final record will name the pushed commit, pull request and merge commit, settled checks, Pages workflow result, and verified live article URL. No deployment will be reported from a local render alone.

## Residual limits

The results apply to the tested repeated-observation scenarios and current source. They do not establish general superiority over non-neural methods or performance across other traits, missingness mechanisms, or empirical datasets. The no-covariate-signal control showed a small loss from the neural contribution.
