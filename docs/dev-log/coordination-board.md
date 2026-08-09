# Coordination Board — pigauto

## Active Lane Split

| Lane | Owner | Branch | Handover / truth | Status |
|---|---|---|---|---|
| **P0 review blockers + B1–B3** | prior Cursor `/goal` | `fix/p0-review-blockers` @ `5800b01` | `docs/dev-log/handover/2026-08-08-p0-rose.md` · PR [#155](https://github.com/itchyshin/pigauto/pull/155) → `fix/ci-install-libtorch` | CI in progress; merge to parent only; **not** to `main` |
| **BACE wrap / re-bench (pigauto remit)** | fresh Cursor | `handover/2026-08-09-cursor` | `docs/dev-log/handover/2026-08-09-cursor-handover.md` | START HERE for wrap G0. Do not edit `BACE/` |

Rehydrate must read **both** rows. A single AGENTS.md snapshot pointer would orphan a sibling — there is no Live Phase Snapshot in AGENTS.md; **this table is the split**.

## Current Rule

- One lane edits pigauto `R/` at a time. Wrap lane does not push onto `fix/p0-review-blockers` while #155 is open.
- Do not modify standalone BACE or in-tree `pigauto/BACE`.
- No rebase onto `main` without a new G0. No EM restore (`max_iter>0`) without a new G0.

## Status

- 2026-08-09: P0 landed on origin; #155 merge-from-parent done; R-CMD-check running. BACE wrap lane opened via Cursor handover. Dirty uinit / DRAC scripts / README banners stay unstaged.
