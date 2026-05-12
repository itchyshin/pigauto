# Codex Cleanup Audit - pigauto

Date: 2026-05-11

Scope: low-risk cleanup planning for a GitHub/pkgdown release. This audit does
not propose algorithm changes, deep refactors, BACE changes, or new benchmark
claims. Findings below distinguish file evidence from command-derived repo
state.

## Executive Summary

pigauto has a coherent package core and a much cleaner public documentation
direction than it had before the pkgdown Phase 1 work. The main release risk is
not "the model is disorganized"; it is that source, generated documentation,
tracked artifacts, and user-facing claims are out of sync.

Most cleanup should be done as small documentation/repository hygiene PRs before
any code reorganization. The safest first sequence is:

1. Resolve the dirty worktree and generated roxygen drift.
2. Rebuild or remove stale generated pkgdown output, especially retired TabPFN
   and DOCS pages.
3. Decide which benchmark HTML artifacts are canonical sources and which are
   local build outputs.
4. Clean internal/public documentation boundaries.
5. Only then do behavior-preserving helper extraction in large files.

## Current State Observed

- Branch: `docs/pkgdown-reorg-phase1` tracking
  `origin/docs/pkgdown-reorg-phase1` (command: `git status --short --branch`).
- Dirty worktree: `DESCRIPTION` and `man/fit_pigauto.Rd` are modified;
  `AGENTS.md`, `script/bench_phase_g_double_prime_conditional.md`, and
  `useful/bench_summary_for_dan_2026-05-04.md` are untracked (command:
  `git status --short --branch`).
- `pkgdown::check_pkgdown()` returned `No problems found.` This is a useful
  reference-index check, not proof that ignored local `docs/` output is current.
- Source/test size from `wc -l`: `R/` has 15,956 lines; `tests/testthat/` has
  9,627 lines. The largest source files are `R/predict_pigauto.R` (1,193),
  `R/fit_pigauto.R` (1,041), `R/preprocess_traits.R` (882), `R/plot.R` (813),
  `R/evaluate.R` (785), and `R/active_impute.R` (743).

## Findings

| ID | Severity | Scope | Evidence | User impact | Fix size | Release-blocking? |
|---|---|---|---|---|---:|---|
| F1 | Required | Worktree hygiene | `git status --short --branch`; `git diff -- DESCRIPTION`; `git diff -- man/fit_pigauto.Rd` | A release PR can accidentally carry roxygen-version churn or malformed Rd text. | Small | Yes |
| F2 | Required | Tracked ignored artifacts | `.gitignore:32-34`, `.gitignore:65-73`, `.Rbuildignore:7-9`, command `git ls-files -ci --exclude-standard` | Future contributors cannot tell which large/local artifacts are source versus leftovers. | Small/Medium | Yes for release branch hygiene |
| F3 | Required | Local generated pkgdown output | `.gitignore:42-43`, `docs/reference/setup_tabpfn.html:2-6`, `docs/reference/fit_baseline_tabpfn.html:2-6`, `docs/sitemap.xml:3-4`, `docs/sitemap.xml:35`, `docs/sitemap.xml:60` | Local docs advertise retired TabPFN/DOCS pages and can be mistaken for the release site. | Small | Yes if `docs/` is deployed or reviewed locally |
| F4 | Required | Methodology navbar/assets policy | `_pkgdown.yml:69-110`, command `comm -23 pkgdown/assets/dev docs/dev` | Evidence pages exist in one tree but are hidden/pending in another, so method readers see an arbitrary subset. | Medium | No, but fix before public evidence push |
| F5 | Required | User-facing docs drift | `DESCRIPTION:3`, `README.md:311-317`, `README.md:350`, `README.md:354-355`, `R/fit_baseline.R:340-347` | README can tell users an implemented Phase F path is still queued and cites an older package version. | Small | Yes |
| F6 | Required | Public/internal reference boundary | `NAMESPACE:21-51`, `R/bm_internal.R:23-29`, `R/bm_internal.R:36-56`, `R/mask_missing.R:190-197`, `man/bm_impute_col.Rd:1-33`, `man/expand_trait_idx_to_latent.Rd:1-22` | Internal helpers appear as installed help topics/local reference pages even though they are not exported. | Small | No, but high polish value |
| F7 | Suggestion | Code organization | `R/predict_pigauto.R:143-245`, `R/predict_pigauto.R:641-740`, `R/predict_pigauto.R:746-879`, `R/predict_pigauto.R:884-1045`, `R/fit_baseline.R:80-86`, `R/fit_baseline.R:179-228`, `R/fit_baseline.R:419-426` | Large files slow review and make safe changes harder, but behavior is test-covered enough to postpone refactor. | Medium | No |
| F8 | Suggestion | Validation/test depth | `tests/testthat/test-shipping-coverage.R:1-9`, `tests/testthat/test-shipping-coverage.R:29-59`, `tests/testthat/test-shipping-coverage.R:130-170` | Some release-hardening tests are intentionally smoke-level; extraction/refactor work needs targeted behavior tests first. | Medium | No |

## Details and Recommendations

### F1. Resolve dirty worktree before cleanup PRs

Evidence:

- `git status --short --branch` reports modified `DESCRIPTION` and
  `man/fit_pigauto.Rd`, plus three untracked files.
- `DESCRIPTION` currently has `Config/roxygen2/version: 8.0.0` at
  `DESCRIPTION:57`, but the diff shows this replaced `RoxygenNote: 7.3.2`.
- The `man/fit_pigauto.Rd` diff changes the continuation line for the
  `safety_floor` equation from an indented line to a line beginning with `+`.

Recommendation:

- Do not bundle these changes into the audit cleanup unless intentional.
- For the Rd change, regenerate with the agreed roxygen2 version rather than
  hand-editing generated files.
- For the roxygen metadata change, choose one policy: keep the current local
  roxygen2 output and commit it explicitly, or restore the previous field before
  release.

### F2. Stop tracking files that the repo now treats as ignored/local

Evidence:

- `.gitignore` says model checkpoints are ignored at `.gitignore:32-34`.
- `.gitignore` says `script/data-cache/` is ignored at `.gitignore:65-68`.
- `.gitignore` says Vulcan bundles are ignored at `.gitignore:70-73`.
- `.Rbuildignore` excludes `avonet`, `checkpoints`, `checkpoints_bin`, `dev`,
  `script`, `submit_v090_vulcan`, and `submit_v090_vulcan_gpu` at
  `.Rbuildignore:6-11` and `.Rbuildignore:26-28`.
- `git ls-files -ci --exclude-standard` still reports tracked ignored files,
  including checkpoint `.pt` files, `script/bench_joint_baseline.rds`, and
  Vulcan scripts.

Recommendation:

- Create a repository-hygiene PR that does only `git rm --cached` for artifacts
  that are not canonical source.
- Keep true package data in `data/` and `inst/extdata/`; those are already
  exceptions in `.gitignore:38-40`.
- If any ignored tracked artifacts are intentionally canonical, document that in
  a short `README` under the relevant directory and tighten the ignore pattern.

### F3. Rebuild or remove stale generated `docs/`

Evidence:

- `_pkgdown.yml` builds into `docs` (`_pkgdown.yml:3`).
- `.gitignore` ignores `docs/` (`.gitignore:42-43`), and `git ls-files docs`
  returns zero tracked files.
- Local generated pages still contain retired TabPFN references:
  `docs/reference/setup_tabpfn.html:2-6`,
  `docs/reference/fit_baseline_tabpfn.html:2-6`, and
  `docs/dev/bench_tabpfn.html:29-37`.
- Local `docs/sitemap.xml` still lists retired pages:
  `CLAUDE.html` and `DOCS.html` at `docs/sitemap.xml:3-4`,
  `fit_baseline_tabpfn.html` at `docs/sitemap.xml:35`, and
  `setup_tabpfn.html` at `docs/sitemap.xml:60`.
- NEWS says those TabPFN pages were removed at `NEWS.md:1926-1931`, and the
  new Phase 1 NEWS says `DOCS.md` was retired at `NEWS.md:24-26`.

Recommendation:

- Before release, delete local `docs/` and run the pkgdown build from source,
  or rely only on the GitHub Pages workflow.
- Add a small developer note or Make target: `rm -rf docs && Rscript -e
  "pkgdown::build_site_local()"`.
- Do not review stale local HTML as evidence unless it has just been rebuilt.

### F4. Decide the methodology-page source of truth

Evidence:

- The active Methodology dropdown lists 7 per-trait pages and 6 cross-dataset
  pages (`_pkgdown.yml:69-97`).
- `_pkgdown.yml:103-110` explicitly says ten older bench HTMLs are pending
  regeneration and are not in the navbar.
- `pkgdown/assets/dev/` has 34 tracked HTML files, including pages absent from
  local `docs/dev/`.
- `comm -23 pkgdown/assets/dev docs/dev` reports missing-from-local-docs pages
  such as `bench_multi_proportion`, `bench_scaling_v090`,
  `bench_bace_avonet_head_to_head`, `bench_pantheria_bace_head_to_head`,
  `phase8_summary`, and `pantheria_summary`.
- The workflow watches `pkgdown/**` at `.github/workflows/pkgdown.yaml:23-25`,
  so these tracked assets are treated as site-affecting source.

Recommendation:

- Pick one policy before release:
  1. Rebuild the missing `docs/dev` pages and re-add the methodology links.
  2. Keep them hidden, but move non-current assets out of the tracked
     `pkgdown/assets/dev/` source tree.
- Do not partially expose benchmark evidence. If a page is not current enough
  for the navbar, label it archived or remove it from the build source.

### F5. Update user-facing drift in README/DESCRIPTION

Evidence:

- `DESCRIPTION:3` says the package version is `0.9.1.9014`.
- README citation still says version `0.9.1.9009` at `README.md:354-355`.
- README says Phase F is queued for v0.9.2 at `README.md:311-317`.
- Source contains the Phase F LP path in `R/fit_baseline.R:340-347`.
- README still links "Architecture notes" to `CLAUDE.md` at `README.md:350`,
  while `.Rbuildignore:23` excludes `CLAUDE.md` from the built package. That is
  fine for GitHub readers, but not for installed-package readers.
- `DESCRIPTION:10-17` describes continuous, counts, binary, ordered, and
  unordered categories, while `_pkgdown.yml:21-29` advertises eight types plus
  uncertainty mechanisms.

Recommendation:

- Update README citation version as part of release prep.
- Replace internal "Phase B3/Phase F queued" language with a current,
  user-facing caveat. If the caveat still holds, state the exact current
  condition without phase labels.
- Either keep the `CLAUDE.md` link clearly contributor-only or point users to
  live pkgdown architecture pages instead.
- Align DESCRIPTION and pkgdown home language so the package description
  mentions the same type support at the same level of precision.

### F6. Clean the installed help/reference boundary

Evidence:

- `NAMESPACE:21-51` exports 31 functions.
- `R/bm_internal.R` marks `phylo_cor_matrix()` and `bm_impute_col()` as
  internal via `@keywords internal` at `R/bm_internal.R:23-29` and
  `R/bm_internal.R:36-56`.
- `R/mask_missing.R:190-197` similarly documents
  `expand_trait_idx_to_latent()` as internal.
- Generated Rd files still exist for these internal topics:
  `man/bm_impute_col.Rd:1-33`,
  `man/expand_trait_idx_to_latent.Rd:1-22`, and
  `man/phylo_cor_matrix.Rd:1-19`.

Recommendation:

- For internal-only helpers, add `@noRd`, run `devtools::document()`, and
  confirm the internal Rd files are removed.
- If a helper is intended for advanced users, export it deliberately, add it to
  `_pkgdown.yml`, and write examples. Do not leave it halfway visible.

### F7. Postpone large-file refactors until after release hygiene

Evidence:

- `R/predict_pigauto.R` contains model reconstruction and prediction setup
  (`R/predict_pigauto.R:143-245`), decoding (`R/predict_pigauto.R:641-740`),
  pooling (`R/predict_pigauto.R:746-879`), and SE computation
  (`R/predict_pigauto.R:884-1045`).
- `R/fit_baseline.R` owns baseline API setup (`R/fit_baseline.R:80-86`),
  Level-C dispatch (`R/fit_baseline.R:179-228`), and OVR categorical dispatch
  (`R/fit_baseline.R:419-426`).

Recommendation:

- Do not refactor these before docs/site hygiene unless a behavior bug requires
  it.
- When ready, extract by responsibility without changing signatures:
  `R/predict_decode.R`, `R/predict_se.R`, and possibly
  `R/baseline_dispatch.R`.
- Make each extraction a separate PR with no semantic edits and targeted tests
  around before/after equality on small fixtures.

### F8. Treat smoke tests as smoke tests

Evidence:

- `tests/testthat/test-shipping-coverage.R:1-9` says the file covers smoke
  tests for exports, non-default option paths, and edge cases.
- Examples include shape/type smoke checks for `confusion_matrix()` and
  `calibration_df()` at `tests/testthat/test-shipping-coverage.R:29-59`.
- Non-default path checks run short end-to-end fits at
  `tests/testthat/test-shipping-coverage.R:130-170`.

Recommendation:

- Keep the smoke tests. They are useful.
- Before any code movement in `predict_pigauto.R` or `fit_baseline.R`, add
  targeted tests for the exact behavior being moved: decoded class/probability
  equality, SE equality, gate override validation, and fallback path selection.

## Quick Wins Patch Queue

### PR 1: Release hygiene and stale generated docs

Files/paths:

- Resolve or revert current `DESCRIPTION` and `man/fit_pigauto.Rd` drift.
- Remove stale ignored local `docs/` before any local site verification.
- Add a short docs-build note to contributor guidance if desired.

Verification:

- `git status --short` contains only intended changes.
- `Rscript -e "pkgdown::check_pkgdown()"` returns no problems.
- Optional: `rm -rf docs && Rscript -e "pkgdown::build_site_local()"`, then
  verify no TabPFN/DOCS pages appear in `docs/sitemap.xml`.

### PR 2: Methodology asset policy

Files/paths:

- `_pkgdown.yml`
- `pkgdown/assets/dev/`
- relevant `script/make_*_html.R` drivers only if regenerating pages

Verification:

- Every Methodology navbar link resolves locally after a clean pkgdown build.
- `git grep -n "bench_tabpfn\\|fit_baseline_tabpfn\\|setup_tabpfn"` returns only
  historical NEWS entries or explicitly archived pages.

### PR 3: User-facing docs alignment

Files/paths:

- `README.md`
- `DESCRIPTION`
- `_pkgdown.yml`
- possibly `vignettes/common-pitfalls.Rmd`

Verification:

- README citation version matches `DESCRIPTION`.
- README no longer says Phase F is queued if the code path is already present.
- `pkgdown::check_pkgdown()` passes.

### PR 4: Internal/public documentation boundary

Files/paths:

- `R/bm_internal.R`
- `R/mask_missing.R`
- generated `man/*.Rd`
- `_pkgdown.yml` only if any helper is intentionally promoted

Verification:

- `devtools::document()` removes unwanted internal Rd topics.
- `pkgdown::check_pkgdown()` passes.
- No exported function is removed in this low-risk pass.

### PR 5: Tiny behavior-preserving extractions

Files/paths:

- `R/predict_pigauto.R` into decode/SE helpers first.
- `R/fit_baseline.R` only after prediction extraction is complete.

Verification:

- Add focused tests before extraction.
- Run targeted test files for predict, fit-predict, mixed-types, and
  multi-impute.
- Run full `devtools::test()` before merging.

## Do Later

- Deep baseline architecture redesign.
- New performance claims or benchmark headline changes.
- CRAN-specific dependency/runtime policy.
- BACE modifications.
- Any public API removal or deprecation.

## Verification Performed During This Audit

Commands run:

- `git status --short --branch`
- `git diff --stat`
- `git diff -- DESCRIPTION`
- `git diff -- man/fit_pigauto.Rd`
- `Rscript -e "pkgdown::check_pkgdown()"` -> no problems found
- Static file inspections with `nl -ba`, `rg`, `find`, `git ls-files`,
  `git ls-files -ci --exclude-standard`, and file line counts

Commands intentionally not run:

- `pkgdown::build_site_local()`: skipped because it writes `docs/`, and the
  audit only needed static evidence plus `pkgdown::check_pkgdown()`.
- `devtools::test()` and `devtools::check()`: skipped because no package code
  was changed and the requested deliverable was an audit report.
- Benchmarks: skipped because this cleanup audit does not change performance
  claims.

## Bottom Line

The first cleanup pass should not touch algorithms. The fastest path to a
cleaner GitHub/pkgdown release is to make the repository tell one story:
source docs, generated docs, tracked assets, NEWS, README, and DESCRIPTION all
need to agree about what pigauto currently is.
