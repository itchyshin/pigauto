# Melissa LIGHT — plan-vs-actual: #157 fix + P0 land (2026-08-15)

Plan: `~/.claude/plans/read-agents-md-and-docs-dev-log-handover-replicated-whisper.md`
(approved via ExitPlanMode) + the `/goal` block. Actual: PR #158 merged, `main` = `3677a85`.

| Axis | Verdict | Notes |
|---|---|---|
| Scope | **adaptive ×1** | Second floor layer in `fit_pigauto()` added mid-execution when the experiment re-run showed the calibrate_gates-only floor could not see the post-refine surface. Within the plan's stated intent ("guards assert on the delivered surface"); recorded in after-task §6.1. |
| Evidence / verification | **adaptive ×2, no drift** | All pre-registered checks ran: suite ×2, regression tests, --as-cran, PR CI. The experiment-script acceptance was met after correcting the script's own mis-specified surface (§6.2) — a measurement fix, not a bar change. --as-cran ran twice (library drift, §6.3). |
| Model routing | **as planned** | 2 children this checkpoint (Explore recon @ Sonnet-class; doc-bundle @ Sonnet), 0 Opus children, inline implementation by the orchestrator. Within the 6/1 budget. Luna/Haiku: none — both slices needed judgment (recorded in the retrofit block). |
| Safety gates | **no drift** | Sweep receipt present (retrofit block). Land gated on the full pre-registered evidence chain; merged only after all CI legs passed. Stop rules never triggered. |
| Public claims | **no drift** | PR body claims match measured evidence verbatim (ratios, check counts). No pigauto-vs-BACE claim. NEWS wording scoped to what shipped. |
| Handoff state | **adaptive ×1** | Docs-checkout commits routed via `GIT_INDEX_FILE` around a stale Cursor-gitWorker `index.lock` (§6.5). All close-out docs landed and pushed. |

**Drift: none.** Five adaptive deviations, all justified and recorded in the after-task §6.
Nothing silently dropped: every DEFER item remains fenced and listed in after-task §5.
