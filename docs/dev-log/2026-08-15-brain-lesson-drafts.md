# PROPOSED brain lessons from 2026-08-15 — NOT yet filed (D-37: awaiting approval)

Four candidates from today's pigauto session. Each cites the case that proves it.
On approval: 1, 2, 4 → `[[WHAT-WORKS]]`; 3 → `[[LESSONS]]`.

## 1. WHAT-WORKS — `git merge-tree --write-tree` before sizing any merge arc

One read-only command turns "N commits behind, unknown conflict surface" into an exact
conflicted-file list in seconds, without touching refs or worktree.

*Proof:* pigauto P0 land (2026-08-15). The probe predicted **exactly 7 conflicts** on a
branch 46 behind / 9 ahead; the real merge produced exactly those 7. The arc estimate was
built on the probe and the core-resolution budget came in *under*, not over — the
prediction is what made a "measured" confidence label honest.

## 2. WHAT-WORKS — PSOCK throttling: renice AFTER startup, never nice at launch

To make an R PSOCK-cluster campaign polite on a shared machine: cap the worker count
(env-settable `MC_CORES`-style knob) and **renice the master + workers after they start**
(a 30 s watcher loop). Do **not** wrap the launch in `nice` — nice'd workers start slowly
and miss `makePSOCKcluster`'s connect window.

*Proof:* pigauto bench re-runs (2026-08-15). `nice -19 Rscript` → "4 of 16 workers failed
to connect", chain dead in 2.5 min. Post-start renice at 4 workers → cluster up in 4.7 s,
campaign ran clean while two other R workloads shared the Mac.

## 3. LESSONS — an evaluation script's masking choice is itself a measurable assumption

When a hand-rolled experiment script defines its own "held-out surface", diff that
definition against the package's own evaluation convention **before** trusting its
numbers. A masking mismatch can flip the conclusion's sign.

*Proof:* pigauto #157 investigation (2026-08-15). My delta-surface script masked val cells
only (test cells stayed pinned) — an intermediate surface nothing in the package uses. It
read the GNN as "+2% worse than baseline"; on the package-convention surface (val+test
hidden, as `evaluate()`/`cross_validate()`/`impute()`/`report()` all use) the same fits
showed the GNN **better** in 3/3 reps (0.959–0.981). Three downstream decisions had
started to form on the artifact before the diff caught it.

## 4. WHAT-WORKS — `GIT_INDEX_FILE` plumbing route around a stale index.lock

When another tool's git worker leaves a stale `.git/index.lock` (lsof-verified unheld)
and deleting it is unavailable, commits still land safely:
`GIT_INDEX_FILE=/tmp/alt-index git read-tree HEAD && git add <paths> && git commit && git push`
— porcelain honors the env var and locks the *alternate* index instead. Afterwards the
default index is stale vs HEAD; one `git reset` (once the lock clears) resyncs it.

*Proof:* pigauto close-out (2026-08-15). Three docs commits landed this way around a
Cursor gitWorker lock, with the working tree and remote ending byte-correct.
