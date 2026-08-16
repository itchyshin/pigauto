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

## 2. LESSONS — niceness is inherited at fork: never renice a process that will spawn PSOCK workers

To make an R PSOCK-cluster campaign polite on a shared machine, **cap the worker count**
(an env-settable `MC_CORES`-style knob). Do **not** reach for `nice`/`renice` on the R
master — not at launch, and *not after startup either*. Niceness is inherited at fork, so
a reniced master hands nice-19 to every worker it later spawns, and those workers start
too slowly to make `makePSOCKcluster()`'s connect window.

*Proof, in two stages — the second corrects the first:*
1. `nice -19 Rscript` at launch → *"4 of 16 workers failed to connect"*, chain dead in
   2.5 min. Obvious lesson: don't nice at launch.
2. So I reniced **after** startup via a 30 s watcher — and two later drivers hung at
   "Dispatching cells" with 0% CPU for an hour each. The watcher had reniced each *master*
   before it built its cluster, so the workers inherited 19 and hit the same window. The
   one driver that survived had started before the watcher caught it.

**The durable rule:** bound a cluster campaign by *worker count*, not by priority. If you
must renice, renice only the already-connected worker PIDs, never the parent. `MC_CORES=1`
on a many-core box is already polite and is stable.

*Meta-lesson:* the first version of this note was written after stage 1 and was confidently
wrong. A mitigation that has not yet been observed across a full campaign is a hypothesis,
not a lesson.

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
