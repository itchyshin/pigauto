# pigauto 0.11 pkgdown-011-repair (approved contract)

This LOOP file is the binding WHAT for the repair run. It **replaces** the earlier visual-gate-only ultra-plan (`pkgdown_visual_gate_866c3f72`) as the execution contract. Do not re-plan.

## Destination

Named visual, contrast, table, and claim defects are gone on the rebuilt 0.11 pkgdown site. A written close report exists. PR #174 remains draft.

## Approved defaults

1. Drop Step 7 GNN-beats-BM table in `vignettes/getting-started.Rmd`; interpret via closed gate on high-signal morphometrics.
2. Tombstone `pkgdown/assets/validation_suite.html` in place (withdrawn banner + `noindex`).
3. Body links in `pkgdown/extra.css` use brand green `#2c5e4f` (claimed 7.4:1 on white).
4. README and home sidebar use Florence cream/amber warning banners, not raw unstyled `[!WARNING]` text.
5. `vignettes/gnn-architecture.Rmd`: drop BACE 42%→72% attribution; stale `phylopars()` wording matches the in-house solver; §5 UQ table becomes scannable blocks.
6. `AGENTS.md` / `CLAUDE.md` out of scope.

## Sequence

S0 inventory → parallel S1–S4 → S5 rebuild + `check_pkgdown()` → S6 see-pass → S7 close.

## Out of scope

#175 merge · manuscript · evidence compute · Dropbox dirty checkout · undraft/merge/push #174 · `R/` `tests/` `NEWS.md` `DESCRIPTION` `man/` · new `_brand.yml` · knitting `eval=FALSE` plots.
