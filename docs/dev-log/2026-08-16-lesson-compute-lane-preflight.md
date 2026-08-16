# Lesson: lane preflight checks the REPO, not the shared COMPUTE HOST

2026-08-16 · pigauto Claude lane · **DRAFT for the brain — staged, not written**
(brain-write boundary: proposing an edit to `~/shinichi-brain`, not making it).

## What happened

Ultra-plan Phase 0.2 ran `tools/lane_preflight.sh` and reported correctly: 3 lanes live in
the *repository*, foreign main-direct lane active, don't touch its files. I took a clean
lane, launched two Totoro campaigns at 100 + 40 workers, and recorded the D-143 arithmetic
in the design docs: 140 ≤ 150. Compliant, on the numbers I had.

Ten minutes later the machine's load average was 235 and `ps` showed the **snakagaw account
at 192 cores** — over D-143's cap. The extra was not mine:

| processes | script | whose |
|---:|---|---|
| 87 | `inst/sim/b1-calibration/run-b1-shard.R` | another repo's lane (gllvmTMB/drmTMB family) |
| 48 | `cov119-driver.R` | another lane, launched ~3 s before I looked |
| ~102 | `lambda_cell.R` | mine |

Two other lanes were running campaigns on the same host, under the same account, and
**nothing in the plan-time checks looks for that**. `lane_preflight.sh` answers "is another
lane editing this repo?" — a different question from "is another lane already consuming the
shared compute budget I am about to spend?"

D-143 is an *account-level* cap on a *shared machine*. A per-lane worker count cannot
satisfy it, because each lane can be individually reasonable and collectively over.

## Why the existing checks miss it

- **Phase 0.2 (Shannon)** is repo-scoped: branches, worktrees, coordination board, ledger IDs.
- **D-139** gates *my* run's duration and correctness (estimate, pre-run, approval). It asks
  "how long will this take?", never "who else is on the box?"
- **D-143** states the cap but names no checker — so it is enforced by whoever happens to
  remember, using only the numbers in front of them.

Per D-90, a binding rule wants a named checker. This one has none.

## Proposed addition to ultra-plan (Phase 0.2b — compute preflight)

Before launching any campaign on Totoro or DRAC, and cheap enough to be unconditional:

```sh
# Totoro: what is the ACCOUNT already using, and who else is on it?
ssh totoro 'ps -eo user,pcpu --no-headers | awk "{a[\$1]+=\$2} END {for (u in a) if (a[u]>50) print u, int(a[u]/100)\" cores\"}" | sort -k2 -rn
  ps -u $USER -o args --sort=-pcpu | grep -E "^/usr/lib/R|Rscript|julia" | sed "s/ --.*//" | sort | uniq -c | sort -rn | head'
```

Then: **my budget = 150 − (what the account is already using)**, not 150. If that leaves too
few workers to finish in reasonable time, that is a **coordination question for Shinichi
(D-87)** — do not silently take the cores, and do not silently abandon the campaign.

Rules this implies:
1. Record the *measured* pre-launch account usage in the design doc beside the D-139
   estimate. "140 ≤ 150" is only true if the account was at ~0.
2. Re-measure after launch. A campaign that was compliant at launch can be pushed over by a
   lane that starts afterwards — as happened here.
3. Never kill another lane's jobs. Throttle your own and surface the overlap (D-87: ownership
   is Shinichi's call).

## Two operational gotchas worth keeping

- **`pkill -f <pattern>` over SSH kills your own session** when the pattern appears in your
  own remote command line (`ssh host 'pkill -f run_campaign.sh'` → the wrapper matches
  itself → exit 255, nothing killed). Use the bracket trick: `pkill -f "run_campaign_arc[2]"`.
- **Per-cell resumable campaigns make throttling cheap and safe.** Killing the dispatcher
  costs at most the in-flight jobs; the skip-existing check means a relaunch at lower
  parallelism loses nothing. That property is what turned a D-143 breach into a 2-minute fix,
  and it is worth building into every campaign runner for that reason alone.

## Postscript (the honest ending)

By the time the throttle landed, the λ campaign had already finished — 1,920/1,920 in ~19
minutes. The breach was real and lasted roughly ten minutes; the remedy arrived after the
fact. That is precisely why this belongs at *plan time* rather than in a monitoring habit:
a short campaign can be over before a reactive check fires.
