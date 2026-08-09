# GOAL — pigauto BACE wrap, Option B-minus (IMMUTABLE)

## Mission
`fit_baseline_bace()` gains an opt-in proper-MI path: `final_imp = FALSE` (default) +
`n_final = 15L`. When `final_imp = TRUE`, call `BACE::bace_final_imp()` on the `bace_imp`
object the wrapper already builds and derive `mu` / `se` from `all_datasets`.

## Headline
Make the Study-B-aligned BACE object reachable from pigauto **without changing the default
path by a single bit**.

## Invariants
- Default path (`final_imp = FALSE`) must be byte-identical to pre-change. Verify against
  `git show HEAD:R/fit_baseline_bace.R`, same seed — not by inspection.
- `n_final = 15L` is Shinichi's explicit AskQuestion pick (the `/goal` paste also used 15).
  Matches Dan's Study B. Do **not** flip to 50 without a new ask.
- Scope: `R/fit_baseline_bace.R` + `man/` + `tests/` + `NEWS.md` + docs-only
  `R/ovr_categorical.R` footnote. Nothing else.
- NEWS must **not** claim this fixes "imputed-as-observed". That defect lived in BACE's
  `bace_final_imp()`, which pigauto never called — no pigauto output was ever affected.
- Slack remains parked per Shinichi's `/goal` paste ("No Slack to Dan"). Do not draft or send.

## Forbidden
Editing either BACE tree · vendor-sync `pigauto/BACE` · touching PR #155 /
`fix/p0-review-blockers` · restoring EM (`max_iter > 0`) · rebase or merge to `main` ·
`git add -A` · staging uinit files or `gnn_attribution` artefacts · `ScheduleWakeup` /
archiving this project.

## Authoritative WHAT
`LOOP/ultra-plan.md` (frozen copy) · source of truth
`docs/dev-log/handover/2026-08-09-bace-wrap-g0-ultra-plan.md`.

## Definition of done
S0–X2 complete; default path proven bit-identical; focused tests green with output pasted;
S3 re-bench numbers exist (numbers only — a public pigauto-vs-BACE claim is a Shinichi gate).
