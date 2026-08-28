# pigauto 0.11 pkgdown visual gate

Date: 2026-08-28. Worktree `/Users/z3437171/local-scratch/lanes/pigauto-011-ship` @ `6fddd79` plus scoped CSS/yml fixes. PR [#174](https://github.com/itchyshin/pigauto/pull/174) remains **draft**.

## Verdict

**PASS (after fail-to-ship CSS/yml fixes).** Florence/Tufte gate, not redesign. #174 may be considered for undraft on visual grounds. This file is an internal verdict; it is not a public “Florence-clean” claim. CI green was never this gate.

## Destination check

A written Florence+Tufte verdict exists for the 0.11 user-facing pages; fail-to-ship defects were fixed in `pkgdown/extra.css` and `_pkgdown.yml`; #175 was not touched; #174 was not undrafted or merged.

## S1 — Local site build

`pkgdown::build_site(devel = TRUE, preview = FALSE)` in this worktree. `devel = TRUE` so the render is the #174 source, not the installed library.

- Exit 0. Warnings only: two `@examplesIf interactive()` are FALSE (not visual).
- `docs/` had to be emptied first: pkgdown refused a non-pkgdown destination that already held tracked `docs/dev-log/` files. Tracked ledgers were restored after the build.
- `pkgdown::check_pkgdown()` after yml edit: **No problems found.**

## S0 — Page inventory

See `docs/dev-log/2026-08-27-pkgdown-visual-inventory.md`. In-gate: home, five articles, reference landing, navbar. Methodology = layout + claim-banner only. Tombstones = withdrawn check only. Get started / pitfalls plot chunks remain `eval=FALSE` (G0 default).

## S2 — Emmy (structure vs `_pkgdown.yml`)

| Check | Result |
| --- | --- |
| Navbar left: Get started, Articles (5), Methodology (live benches), Reference, Changelog | PASS |
| Tree-sensitivity bench **not** in navbar; article label “not inferentially validated” | PASS |
| Reference groups match yml (Start here → experimental MI → stochastic diagnostics → …) | PASS |
| Tombstones still say “Historical development page withdrawn” (sampled `tests_overview`, `bench_bien`, `bench_avonet9993_bace`) | PASS |
| `dev/bench_tree_uncertainty.html` is a withdrawn historical page, not a live inferential bench | PASS |
| Root `AGENTS.md` / `CLAUDE.md` / `VALIDATION_LEDGER.md` / `goodagents.md` built as orphan HTML | polish-later (not in navbar) |

## S3 — Florence / Tufte see-pass

Seen: 1400×900 screenshots of home, getting-started, common-pitfalls, mixed-types, tree-uncertainty, gnn-architecture, reference, `dev/bench_continuous.html`, `dev/tests_overview.html`; plus home after the CSS fix. Layout widths measured via headless DOM (`innerWidth` 500; no element `scrollWidth` > viewport).

| Check | Bucket | Notes |
| --- | --- | --- |
| Contrast, body text / navbar white on `#2c5e4f` | PASS | Brand green + Inter |
| Hierarchy (H1 / H2 / code / TOC) | PASS | |
| Code blocks styled, not raw | PASS | JetBrains Mono |
| Mixed-types type table: light rules, no grid junk | PASS | |
| Article hex logos overlap or clip title (desktop) | PASS | |
| Get started “Visualise” is unevaluated code | polish-later | G0: not an automatic fail |
| Methodology live benches: own CSS, claim banner present, “not a general dominance claim” | PASS (layout/claim only) | Not a figure rewrite |
| Tombstone pages look withdrawn, not live benches | PASS | |
| Literal `[!WARNING]` on home (GFM alert not converted by pkgdown) | **fail-to-ship → fixed** | now a peach callout; token stripped |
| Forced 260px home logo (`extra.css` pre-fix) | **fail-to-ship → fixed** | `max-width` only; stack on small screens |
| Long `<pre>` stretching Bootstrap flex columns | **fail-to-ship → fixed** | `min-width: 0` on row/cols |
| Empty `alt=""` on pkgdown logos | polish-later | pkgdown default |
| Search placeholder “Search for” | polish-later | |
| Live `pkgdown.yaml` badge can read failing | polish-later | workflow skipped on PRs by design; not a page-layout defect |
| 390×844 headless PNGs looked clipped | **not a page bug** | Chrome laid out at 500px CSS and cropped the bitmap to 390 |

## S4 — Rose claim fence

Articles and home were grepped for `#175`, sentinel, attribution, outperform / beats / state-of-the-art, certified coverage. **No #175 evidence bleed.**

- Home + sidebar: experimental; point estimates supported; intervals not certified.
- Tree article: descriptive sensitivity only; do not `pool_mi()`; not analysis-aware MI.
- GNN article 33.7% RMSE lift is **caveated** as pre-in-house-solver and not a current shipped-path number.
- Methodology continuous bench banner: scenario-specific deltas, not general dominance.
- Reference “Experimental analysis-aware MI” vs “Stochastic prediction diagnostics” split is intact.

No new public capability claim is made by this pass.

## S5 — Fixes applied (visual-audit lease only)

1. `pkgdown/extra.css` — header flex, logo max-width, mobile stack, column `min-width: 0`, blockquote callouts, pre overflow.
2. `_pkgdown.yml` `template.includes.after_body` — strip leftover GFM `[!WARNING]` / `[!NOTE]` / … tokens and classify the blockquote.

Vignette / README source left to Codex #174 (do not take that lease). Suggested ticket: replace README `[!WARNING]` with a normal paragraph so the JS strip is unnecessary.

## S6 — Re-verify

- `pkgdown::build_home()` after yml change; `extra.css` copied into `docs/`.
- Home re-seen: experimental callout is a real warning box; no `[!WARNING]` token.
- `pkgdown::check_pkgdown()`: no problems.
- `gh pr view 174`: `draft=true state=OPEN`.

## S7 — Close

Plan-vs-actual: `docs/dev-log/plan-actual/2026-08-27-pkgdown-visual-gate.md`.

**OPEN GATE (human):** undraft / merge of #174. This lane stops.
