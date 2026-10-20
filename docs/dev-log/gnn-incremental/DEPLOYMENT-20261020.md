# GNN incremental-value push and deployment record

## Authorization and scope

On 2026-10-20, after the campaign's local closeout, the user explicitly authorized pushing the committed branch and deploying the pkgdown update. The repository is public. Its Pages workflow deploys from `main`; reviewed changes went through a pull request, and merge waited for all required checks to settle successfully. This follow-up adds no simulation or public API, method, or default changes. It records the user-authorized site deployment.

## Starting evidence

The Totoro campaign completed at 05:20 MDT with 2,500/2,500 datasets, 5,000/5,000 pigauto gate fits, and 30,000 successful method attempts. The strict audit and summary regeneration passed, with no failures, invalid outputs, skips, or nonconvergence. The workers exited. The scoped campaign ledger was independently reverified with all eight gates met, including a successful pkgdown build and structural check using a 600-second timeout.

The feature branch was pushed to `origin/codex/gnn-incremental-value` at commit `182717b19b253a87d96217570be4abbf75af4c80`. Pull request [#240](https://github.com/itchyshin/pigauto/pull/240) merged to `main` as `78b1306d18d2a27b6f1caef574c0782c1cc3b4dc` on 2026-10-20. Ubuntu R release, Ubuntu R devel, and macOS R release checks for the exact PR head completed successfully. The PR pkgdown job was skipped by the repository workflow design; the main push triggered the deployment workflow.

## Release outcome

Complete. The `pkgdown.yaml` workflow run [#758](https://github.com/itchyshin/pigauto/actions/runs/39474883346), triggered by the push of the merge commit to `main`, completed successfully. Its `pkgdown` job and `Deploy to GitHub Pages` step both succeeded. The live article at https://itchyshin.github.io/pigauto/articles/validation-results.html was opened on 2026-10-20. It showed the five-scenario table, the measured neural-on/off RMSE and relative reductions, the small null-control loss, and the explicit limits on general superiority. The deployed `main` SHA equals the merge commit above. No direct or manual deployment was used. A local render alone does not count as deployment evidence.

## Deployment gates

- D1 passed: remote PR head `182717b19b253a87d96217570be4abbf75af4c80` matched the reviewed commit.
- D2 passed: PR #240 merged through the review path after Ubuntu R release, Ubuntu R devel, and macOS R release checks settled successfully.
- D3 passed: main `pkgdown.yaml` run #758 and its GitHub Pages deployment step completed successfully for merge commit `78b1306d18d2a27b6f1caef574c0782c1cc3b4dc`.
- D4 passed: the live article was inspected and matched the tested table and stated limits.

## Residual limits

The results apply to the tested repeated-observation scenarios and current source. They do not establish general superiority over non-neural methods or performance across other traits, missingness mechanisms, or empirical datasets. The no-covariate-signal control showed a small loss from the neural contribution.
