# GNN incremental-value push and deployment gates

Scope: following completion of the scientific campaign, the user explicitly authorized pushing the committed branch and deploying the pkgdown update on 2026-10-20. No further simulation, public API change, method change, default change, or publication is authorized by this follow-up.

- [ ] D1: The reviewed campaign branch is pushed without force to the intended remote branch, and its remote head equals the locally verified commit.
  REVIEWER: manual
  EVIDENCE: pending

- [ ] D2: The branch is merged to `main` through the repository's review path only after all checks for the exact PR head have settled successfully; no protection rule is bypassed.
  REVIEWER: manual
  EVIDENCE: pending

- [ ] D3: The pkgdown workflow for the merge commit completes successfully on `main`.
  REVIEWER: manual
  EVIDENCE: pending

- [ ] D4: The live validation-results article contains the measured GNN table and its scenario limits, and the deployed site corresponds to the verified `main` commit.
  REVIEWER: manual
  EVIDENCE: pending
