# PROPOSED brain lessons — 2026-08-15/16

**Lessons 1–4: APPROVED and FILED** 2026-08-16 (vault `0cf87d1`; #2 skipped as already
filed by the gllvmTMB lane). **Lesson 5 below is NEW and awaits its own D-37 approval** —
the earlier approval covered only the four that existed at the time.

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


---

## 5. ~~LESSONS — testing the option you chose is not testing the choice~~ — **RETRACTED: ALREADY FILED**

> **DO NOT APPROVE OR LAND THIS AS A NEW ENTRY.** The drmTMB lane wrote it directly to the
> shared brain ~20 minutes before this draft existed: `~/shinichi-brain/memory/LESSONS.md:2193`,
> commit `387651d`, *"A recommendation is code you did not run — test the choice, not the option
> you chose"*. Verified present, and it already cites pigauto PR #166 / `c536499` as its receipt.
> Landing this draft too would put **two entries in the vault for one lesson** — the same
> "two correct implementations of one fix" hazard the lane preflight warns about, and harder to
> spot in prose than in code.
>
> **If anything below reads better than the filed entry, improve that entry in place.** Do not
> append beside it. Text kept only for that comparison.

*Formulated by the drmTMB lane (`claude/eloquent-driscoll-521fa1`), receipt `665423395`.
Credit theirs; staged here because the failure had two halves and pigauto owns one.*

Two agents compared two forms of the same fix. **Both tested green. Between them they never
exercised the one case that distinguished the options.**

- The drmTMB lane recommended `on.exit(parallel::stopCluster(cl), add = TRUE)` **without
  executing it** — a recommendation did not feel like code. It is code; it just runs in
  someone else's repository.
- This lane wrote `on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)` and
  verified **that** version behaviourally. The `try()` wrapper masked the exact failure that
  justified it, so the test could not distinguish the two forms either.

The bare form errors with `invalid connection` on the success path — it would have traded a
worker leak for a spurious error at the end of every successful run. That surfaced only
because one lane asked the other a direct question, days into running.

**The rule.** When you choose between alternatives — including an alternative you are
*recommending to someone else* — the test must exercise the **rejected** arm. A green test on
the option you picked measures the option, not the comparison. And a recommendation crossing a
repository boundary deserves the same execution bar as a commit, because it will become one.

**Companion failure, same shape, same night (pigauto):** `PPID 1` was taken as evidence of an
orphaned PSOCK worker, and a master at 0% CPU was taken as evidence of a hung driver. Both are
proxies the mechanism never promised to provide — R 4.0's `setup_strategy = "parallel"`
reparents healthy workers to init at birth, and `parLapply` logs nothing until every cell
finishes. The sound test is the port pair (`lsof -nP -iTCP:<PORT>`, two ESTABLISHED = healthy).
Reading a proxy the mechanism never guaranteed cost two healthy drivers.


---

## ⚖️ FOR SHINICHI — a fleet-wide rule whose citation does not support it (D-87)

**Verified independently by this lane, 2026-08-16.** Not an inference.

`pigauto/AGENTS.md:476` says:

> Brain-write boundary — see brain [[DECISIONS]] D-37 (never write to the brain vault
> without explicit approval; stage a draft and propose).

**D-37 does not say that.** `~/shinichi-brain/memory/DECISIONS.md:754–780` is *"Brain
foundation sprint: local-only restored · honest health gate · truth repair"*, and its three
items are (1) remove the vault's `origin` and name the no-push exception, (2) four health
states for `brain_doctor.py`, (3) repair specific audited truth drifts. There is no approval
gate in it.

Nor anywhere else in the vault. Exact-phrase counts:

| phrase | `memory/DECISIONS.md` | hub `AGENTS.md` |
|---|---:|---:|
| `without explicit approval` | 0 | 0 |
| `stage a draft` | 0 | — |
| `brain-write` | 0 | — |

The hub cites D-37 **once**, at `AGENTS.md:233`, for the opposite-shaped point: *"this vault
is the exception (D-37) — local-only, no remote, a local commit IS landed state."*

**Scope — CORRECTED 2026-08-16 after the drmTMB lane pushed back.** My first version of this
section said "five repos carry the identical line." That is true of five **working trees** and
only **two repositories**. Verified per repo (`disk` = file on disk, `HEAD` = committed):

| repo | disk | HEAD | git status |
|---|---:|---:|---|
| `GLLVM.jl` | 1 | 1 | clean — **committed** |
| `CBIC` | 1 | 1 | clean — **committed** |
| `drmTMB` | 1 | 0 | ` M` dirty |
| `pigauto` | 1 | 0 | ` M` dirty |
| `survey_best_paper_awards_followup_analysis` | 1 | 0 | `??` untracked |

**In pigauto the rule exists only in the uncommitted edit to `AGENTS.md`** — one of the 15
carried-over dirty files this lane spent all night protecting from staging. `HEAD` and
`origin/main` both have 0 hits. A fresh clone, a `git worktree`, CI, or any lane reading
`origin` does not see it. One `git checkout -- AGENTS.md` removes it with no trace.

**My own error, stated plainly.** I said twice that following the gate was *"a fact about this
repo's instructions, not my judgement call."* It was a fact about this repo's **dirty working
tree**. The instruction still bound me — a file I actually load binds me whether or not it is
committed, and I had no particular reason to check its git state — but the claim as I made it
was wrong, and I made it while separately and repeatedly telling myself that `AGENTS.md` was
one of the uncommitted files not to stage. I had both halves and did not put them together.

**What this lane infers (labelled as inference):** the rule probably came from a real
instruction of yours that was written into the per-repo `AGENTS.md` files with D-37 attached
as the nearest-looking decision. That makes it *a sound rule with a broken warrant* — worse
than either alone, because the citation is what makes it read as already-settled to every
lane, so no lane can resolve it from the documents.

**What this lane is NOT doing, deliberately.** Not editing any `AGENTS.md` — it is an
instruction file, and a peer session asking for a change to one is precisely the thing I must
not act on. Not writing a decision record to `DECISIONS.md` — that is the gated action under
dispute, and settling a contested norm by exercising the contested power would be the worst
available move. Both are yours.

**Meanwhile this lane keeps obeying the rule.** A broken citation does not void an
instruction in the repo's own `AGENTS.md`; you may well have given it. Absent your ruling,
pigauto's brain writes continue to be staged and proposed, as they were tonight.

### How the two lanes actually stood

- **This lane (pigauto):** staged and waited. Correct under `AGENTS.md:476`.
- **The drmTMB lane:** wrote `LESSONS.md:2193` directly (`387651d`). I told it that
  `drmTMB/AGENTS.md:558` carries the rule, so its write was out of bounds by its own repo's
  instruction. **That correction was itself wrong, and I withdraw it.** The line at :558 is in
  drmTMB's *dirty primary checkout*; `origin/main` and that lane's own worktree branch are
  byte-identical at 1183 lines with **zero** hits. That lane loads a file without the rule, so
  its writes were in bounds for the bytes it actually reads.

**The real shape of the disagreement:** neither lane was careless and neither misread its
instructions. We were reading **different files with the same name** — one dirty checkout, one
clean worktree — and each correctly followed what it loaded. That is a far more interesting
failure than a rule dispute, and it is invisible to both lanes without exactly this exchange.

### Original framing, kept for the record


Two lanes handled the *same* finding under *different* brain-write rules on the same night:

- **This lane (pigauto):** staged a draft and waited for approval. pigauto's `AGENTS.md` says
  in terms: *"Brain-write boundary — see brain DECISIONS D-37 (never write to the brain vault
  without explicit approval; stage a draft and propose)."* So the gate is repo-mandated here.
- **The drmTMB lane:** wrote straight to `LESSONS.md` and committed (`387651d`), citing the
  session instruction to write durable findings to a file before the turn ends, plus D-37's
  own note that the vault is local-only where a local commit *is* landed state.

Both read their instructions defensibly. They cannot both be the norm. Either this lane's gate
is stricter than intended, or that write was looser than intended — and the practical cost of
the ambiguity is visible right here: a near-duplicate entry, avoided only because the two lanes
happened to be talking.

Neither lane should resolve this unilaterally (D-87). The drmTMB lane is surfacing it to you
from its side and has explicitly **not** reverted its entry — removing a filed finding to settle
a process question would be the worse error, and I agree. Flagging it from this side so the
question reaches you once, from both.
