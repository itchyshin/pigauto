# After-task — pigauto BACE wrap merge to main

Date: 2026-08-13
Lane: wrap merge G0 (Cursor)
PR: https://github.com/itchyshin/pigauto/pull/156 **MERGED**
Merge commit: `416561b` on `origin/main`

## 1. Goal

Land GitHub-dev BACE wrap (Option B-minus) on `main` after CRAN pigauto 0.10.0. Do not CRAN-submit. Do not merge P0 parent. Do not vendor-sync BACE.

## 2. Implemented

- NEWS lead-in on wrap worktree `.worktrees/bace-wrap-restore`: CRAN 0.10.0 shipped 2026-07-30; next CRAN cut must drop Suggests BACE or wait. Commit `a54e6a4`.
- Undrafted #156; retitled “Restore BACE wrap on GitHub-dev after CRAN 0.10.0”.
- Merged with merge commit `416561b`.
- Verified `origin/main`: Version `0.10.0.9000`, Suggests includes `BACE`, `R/fit_baseline_bace.R` present.

## 3. Evidence

- CRAN DESCRIPTION Date/Publication 2026-07-30 17:10:22 UTC.
- `gh pr view 156` → `state: MERGED`, `mergedAt: 2026-08-13T23:59:24Z`.
- Restore-SHA R-CMD-check (2026-08-09): ubuntu release/devel and macos release SUCCESS. NEWS commit merged while checks were UNSTABLE (re-running). Not claimed as a new `--as-cran` green.

## 4. Not done (parked)

P0 parent → `main` · CRAN-submit · vendor-sync · EM · Slack · DRM.jl · public pigauto-vs-BACE claim.

## 5. Files

Wrap branch: `NEWS.md` (`a54e6a4`). Docs on `handover/2026-08-09-cursor` (this after-task, board, 2026-08-13 handover). Dirty uinit / GNN on that checkout left unstaged.
