# Plan-vs-actual — pigauto BACE wrap CLOSEOUT (Melissa LIGHT)

Date 2026-08-09 · Reconciler: Melissa lens (closeout arc)
Plans: `LOOP/ultra-plan.md` (frozen) ·
`docs/dev-log/handover/2026-08-09-bace-wrap-g0-ultra-plan.md` · locked `/goal`
closeout paste.
Prior X1: `docs/dev-log/plan-actual/2026-08-09-bace-wrap-reconcile.md` (HEAD of
implementation slice, handover `3bfd740` era).

Materiality: scope · evidence/verification · model routing · safety gates ·
public claims · handoff state. Cosmetic wording skipped.

Rule: ≥2 NOT-DONE verdicts withhold a whole-plan "done" claim (D-43).

## Git reality (receipts)

| Object | SHA / URL |
|---|---|
| `origin/main` | `bf46991` |
| Wrap API | `b57da54` (handover) |
| Loud `bace_final_imp` failures | `3bfd740` (handover) |
| Restore on main's parent | `b180555` `feat/bace-wrap-restore` |
| Restore worktree | `/tmp/pigauto-bace-wrap-restore` clean, tracking origin |
| Handover parent before this closeout commit | `88357ef` `handover/2026-08-09-cursor` |
| Draft PR | https://github.com/itchyshin/pigauto/pull/156 (`isDraft: true`, base `main`) |
| Standalone BACE | `ce8bc87d5b1dd4059bacfba23ab42c3f4dfe6080` |

Focused tests this closeout (paste, not memory):

```
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-fit-baseline-bace-final-imp.R")'
SKIP: 'test-fit-baseline-bace-final-imp.R:127:3' ----------
Reason: installed BACE does export bace_final_imp(); nothing to assert
[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]
```

## Row-by-row (S0–X2 + restore + closeout)

| Row | Planned | Actual | Verdict | Tag |
|---|---|---|---|---|
| SCAFFOLD | LOOP kit | LOOP present; P0 kit archived | **DONE** | — |
| S0 | Install standalone BACE `@ce8bc87` | SHA matches; in-tree clone untouched | **DONE** | — |
| R1 | Recon `bace_final_imp` shape | Written recon; no BACE edits | **DONE** | — |
| S1 | `final_imp` + `n_final`, default byte-identical | `b57da54` + restore fold at `b180555` | **DONE** | — |
| S2 | OVR footnote docs-only | On handover and restore | **DONE** | — |
| M2 | Focused test file | 18 pass / 1 skip, re-run this closeout | **DONE** | — |
| M1a | Paths-scoped drift, no rebase | Ran; discovered `b615579` deleted wrap on main | **DONE** | adaptive (finding changed landing) |
| M1b | Focused tests + T4 | Focused file green; T4 restored on restore branch | **DONE** | — |
| M1c | `devtools::check()` clean on wrap surface | Wrap Rd `\usage` OK. Full check **2 WARN / 3 NOTE**, pre-existing. Not 0/0/0 | **NOT DONE** (parked) | adaptive — pre-existing, not wrap-caused; claiming check-green would be drift |
| M1d | Default path unchanged | Bit-identity receipt on handover S1; formals test still green | **DONE** | — |
| S3 | Re-bench wrap both paths vs `BACE::bace()` | Wrap-config vs wrap-config coverage 0.672 → 0.940; **not** vs `BACE::bace()` | **DONE (regime-fenced)** | adaptive — harness scope; public-claim axis held |
| S3b | Unplanned robustness | `3bfd740` folded into restore | **DONE** | adaptive |
| N1 | NEWS, no overclaim | Restore NEWS has the non-fix sentence | **DONE** | — |
| X1 | Melissa reconcile | First file exists | **DONE** | — |
| X2 | After-task + brain-write **proposal** | Proposal file exists; vault untouched | **DONE (proposal only)** | — |
| KEEP | Keep wrap; no merge | Locked; landing executed; still no merge to main | **DONE as locked** | adaptive — landing path was not in Phase 2 table |
| RESTORE | New branch from `origin/main` | `feat/bace-wrap-restore` @ `b180555` | **DONE** | adaptive (Shinichi A-after-#155) |
| CLOSEOUT | After-task + Melissa LIGHT + draft DO-NOT-MERGE PR; stop | This file + after-task + PR #156 draft | **DONE** (closeout only) | adaptive — restore checkpoint said "no PR this turn"; closeout `/goal` required draft |
| wrap→main | Implicit in "landing"; DEFER fence said no merge to main | **Not merged.** Draft PR is wait-for-CRAN | **NOT DONE** (parked / deferred) | adaptive — Shinichi + CRAN surface |

## Material deviations (six axes)

1. **Scope — S3 target.** Planned wrap vs `BACE::bace()`. Actual wrap-default vs
   wrap-`final_imp`. Tag: **adaptive**. Owner: Ada. Public-claim fence held.
2. **Scope — S3 harness path.** Planned reuse `script/bench_avonet_bace.R`.
   Actual purpose-built harness under gitignored `docs/dev-log/`. Tag:
   **adaptive**. Owner: Ada. Residual: not reproducible from `script/`.
3. **Evidence — M1c.** Planned check-clean. Actual wrap Rd OK, full check still
   2 WARN / 3 NOTE. Tag: **adaptive** (pre-existing). Owner: Rose. Must not be
   narrated as check-green.
4. **Safety gates — wrap→main.** Plan DEFER + later KEEP said do not merge.
   Restore checkpoint preferred "no PR." Closeout `/goal` required draft
   DO-NOT-MERGE PR. Tag: **adaptive**. Owner: Ada. Merge still forbidden.
5. **Public claims.** No pigauto-vs-BACE sentence shipped. NEWS + PR body
   repeat the fence. Tag: none (plan held).
6. **Handoff state.** Restore executed on isolated worktree; handover left
   dirty (uinit + gnn). Closeout docs land on handover. Tag: **adaptive**.
7. **`n_final` provenance.** Frozen ultra-plan G0 table: Ada's IF-YOU-DO-NOT-MIND
   default. Later AskQuestion: Shinichi picked 15. `LOOP/GOAL.md` corrected;
   X2 brain proposal still says Ada. Tag: **unclear** on the proposal file
   (stale), **adaptive** on GOAL/checkpoint. Owner: Rose (vault write still
   D-37 gated).
8. **Slack.** Ultra-plan later G0: "draft for him to send." `/goal` paste:
   "No Slack to Dan." Actual: parked, no draft sent this closeout. Tag:
   **adaptive** (user override). Owner: Ada.

9. **Handoff — `docs/` gitignore.** `.gitignore` has a blanket `docs/` (pkgdown
   site). X1 reconcile + X2 proposal lived on disk but were **never tracked**.
   Closeout force-adds after-task + both reconciles + X2 proposal by explicit
   path (same pattern as `coordination-board.md`). Tag: **drift** (X1/X2 claimed
   DONE without a git object) → **adaptive** this closeout (force-add). Owner:
   Rose / Ada.

No other silent drops. Parked items remain named.

## D-43 verdict

Two load-bearing NOT-DONE remain (plus the now-corrected X1/X2 untracked-docs drift):

1. M1c full check not 0/0/0.
2. wrap→main deferred (PR #156 is draft / DO NOT MERGE).

Therefore: **do not claim the whole ultra-plan done.** Implementation + restore
+ closeout of the GitHub-dev wrap lane are complete. CRAN landing and
check-green are not.

## DECISION RECEIPT

```
DECISION RECEIPT
  Questions asked
    — Phase 0.4: wrap API chain-only vs opt-in final_imp (G0).
    — AskQuestion: n_final 15 vs 50; Slack draft vs park; KEEP wrap vs drop;
      landing timing vs #155.
    — Closeout `/goal` (this arc): not a new ask; execute closeout then STOP.
  Answers received
    — "it should be final yes — OK opt-in cool /ultra-plan it"
    — B-minus; n_final = 15 explicit; KEEP wrap; landing A after #155;
      GitHub-dev until CRAN; `/goal` "No Slack to Dan".
  Defaults accepted
    — Closeout paste: draft DO-NOT-MERGE PR (restore checkpoint had preferred
      branch-only). Ran the paste.
    — check-after-task 12-section protocol (section 12 required) even though
      the paste said "11 sections" — followed the protocol file on disk.
  Adaptive decisions
    — S3 wrap-vs-wrap not wrap-vs-BACE::bace; S3b loud errors; restore from
      origin/main rather than rebase handover; draft PR #156; n_final
      provenance correction in GOAL.md; no vault write (D-37).
  Unresolved → Rose / Shinichi
    — When may #156 leave draft (v0.10.0 shipped and/or BACE on CRAN)?
    — Promote S3 harness to script/ (needs scope OK)?
    — Approve X2 brain DECISIONS.md entry (and fix n_final attribution first)?
    — M1c pre-existing WARN/NOTE (other lanes; not wrap closeout).
```
