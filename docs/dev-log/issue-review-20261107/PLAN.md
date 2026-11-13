# Open issue review, 7 November 2026

Goal: address actionable open reports, verify the changes, and reply to each
reporter for confirmation. Authors retain issue closure.

Owner: Codex, branch codex/address-open-issues-nov07, isolated from the dirty
handover checkout. Base origin/main 51e9a372. Lease codex:pigauto:nov07-issues.

1. Inventory: ten open issues and their comments read; existing PRs 252, 253,
   255 inspected. Preserve 252 and 253 contributions in this branch with author
   attribution. CI expansion in 255 is outside this change.
2. Repair factor-covariate lambda likelihood and prediction rank handling.
   Accept when equivalent full dummy and treatment designs agree in meaningful
   regression checks. Preserve full GNN covariate encoding.
3. Clarify released versus development evaluation, data-only warning,
   ordinal level input, inference documentation and README flow.
4. Add worked missingness views and useful baseline-only report charts.
5. Regenerate manuals, run focused tests then full suite and package check,
   render changed vignettes, inspect output and obtain independent review.
6. Publish reviewable changes and evidence-specific replies. No issue-closing
   keywords or direct closure. No CRAN submission without its release gate.

Expected implementation/check time: roughly 45 to 90 minutes, subject to the
first full-suite timing. Small baseline fits take seconds; no new GNN or
posterior campaign is required. Serious benchmarks are excluded pending a
bounded compute plan and applicable approval.

Prior-work receipt: GitHub issue search and all ten comment threads; local
origin/main log; PR metadata and file lists for 252, 253, 255; published GNN
campaign record docs/dev-log/gnn-incremental/RESULTS-20261020.md. The broad
external comparator campaign is not completed by this review. The private
mammal tree-1 ties require a shareable reproducer and remain unresolved.

Acceptance: tests and source checks support each claimed fix, documentation
states uncertainty limits, each issue receives a specific response, and all
issues remain open for their authors. Acknowledging an unresolved report is
not resolving it. No brain-vault changes are authorized.
