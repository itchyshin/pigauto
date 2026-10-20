# GNN incremental-value push and deployment gates

Scope: following completion of the scientific campaign, the user explicitly authorized pushing the committed branch and deploying the pkgdown update on 2026-10-20. No further simulation, public API change, method change, default change, or publication is authorized by this follow-up.

- [x] D1: The reviewed campaign branch is pushed without force to the intended remote branch, and its remote head equals the locally verified commit.
  REVIEWER: manual
  EVIDENCE: PR #240 remote head was `182717b19b253a87d96217570be4abbf75af4c80`; local reviewed head matched.

- [x] D2: The branch is merged to `main` through the repository's review path only after all checks for the exact PR head have settled successfully; no protection rule is bypassed.
  REVIEWER: manual
  EVIDENCE: PR #240 merged to `main` as `78b1306d18d2a27b6f1caef574c0782c1cc3b4dc` after Ubuntu R release, Ubuntu R devel, and macOS R release checks completed successfully.

- [x] D3: The pkgdown workflow for the merge commit completes successfully on `main`.
  REVIEWER: manual
  EVIDENCE: main `pkgdown.yaml` run #758 (`39474883346`) completed successfully; its `Deploy to GitHub Pages` step succeeded for merge commit `78b1306d18d2a27b6f1caef574c0782c1cc3b4dc`.

- [x] D4: The live validation-results article contains the measured GNN table and its scenario limits, and the deployed site corresponds to the verified `main` commit.
  REVIEWER: manual
  EVIDENCE: Live article https://itchyshin.github.io/pigauto/articles/validation-results.html inspected on 2026-10-20. The five-scenario table, measured values, and scenario limits were present.
