# Repository Hygiene Audit - 2026-05-11

Scope: repository source-of-truth cleanup for files that are tracked by git
while also being marked local/ignored by `.gitignore` or excluded from package
builds by `.Rbuildignore`. This pass does not rerun benchmarks and does not
change package behavior.

## Fixes Applied

### 1. Removed tracked ignored artifacts from git

`git ls-files -ci --exclude-standard` reported tracked files under paths that
the repo already treats as local-only:

- model checkpoints: `checkpoints/*.pt`, `checkpoints_bin/*.pt`
- generated benchmark result: `script/bench_joint_baseline.rds`
- Vulcan campaign bundles: `submit_v090_vulcan/`,
  `submit_v090_vulcan_gpu/`

These are now removed from the git index with `git rm --cached`, so they remain
available in the local working tree but will no longer be shipped in the main
branch history after this cleanup. A follow-up `git ls-files -ci
--exclude-standard` returned no files.

### 2. Removed duplicate ignore rule

`.gitignore` had two adjacent `script/data-cache/` entries. The duplicate was
removed while preserving the broader comment that benchmark data caches are
downloadable and non-portable.

### 3. Fixed a stale helper reference

`script/make_bench_avonet9993_bace_html.R` previously told users to run a
helper script inside `submit_v090_vulcan_gpu/`. Since that bundle is local /
branch-specific, the error message now says to copy the Vulcan GPU result RDS
into `script/` instead of naming an ignored helper path.

## Reviewed Without Patch

- Local `docs/` is ignored and untracked, but the working tree still contains
  generated stale pages such as `DOCS.html`, `CLAUDE.html`, TabPFN references,
  and legacy benchmark reports. Because `docs/` is not tracked, deleting or
  rebuilding it is a local developer action rather than a durable package patch.
- Tracked `pkgdown/assets/dev/*.html` and `script/*.html` remain in the repo by
  design for now. The methodology-claim pass in PR #72 tightened active
  rendered pages; hidden / pending assets still need a separate audit before
  re-entry into the navbar.

## Remaining Work

1. Decide whether to keep generated benchmark HTML artifacts tracked long-term
   or move toward generated-on-demand site assets.
2. Delete and rebuild local ignored `docs/` before any manual site review:
   `rm -rf docs && Rscript -e "pkgdown::build_site_local()"`.
3. Audit hidden methodology pages before exposing them:
   `bench_multi_proportion`, `bench_scaling_v090`, Phase 8 sweeps, and BACE
   head-to-head pages.
