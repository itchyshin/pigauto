GOAL: see GOAL.md (+ arcs.md).   STATE: campaign running on nibi + fir; rorqual/narval validating (2026-09-24 ~17:50 nibi).
ARCS DONE (verified): pre-run (240/240), freq timing (36/36), results page v1.
ARC IN PROGRESS: A1 (nibi 22634021), A2 (nibi 22634022/22634023), A3 (nibi 22633907/909, root pigauto_rubin_freq),
  A4 freq campaign on fir (61365818/19/20; n1000 seeds 1-200, n100/300 seeds 51-200; root ~/pigauto_rubin; env ~/pigauto_sim/env.sh).
  rorqual/narval: R lib copied from fir (BACE digest identical), env ~/projects/def-snakagaw/snakagaw/pigauto_rubin/env.sh,
  validation jobs rorqual 21771783, narval 3935076 -> if G-S4a PASS, they take part of A5.
STATUS: bash script/rubin_status.sh (Mac). Totoro: this user 155 cores busy in other lanes (>150, D-143) -> not used; recheck.
NEXT: A3 done -> pull nibi pigauto_rubin_freq/results/freq -> compare.R -> interim page; A1 done -> A5 split by seed across
  nibi / fir / rorqual / narval (per-host seed ranges, never overlapping; pool per host subdir, dedupe by filename).
OPEN GATES: A5 overrun gate (A1 > 30% over 6.9 h/fit). Never contact Dan; never edit R/ or BACE/.
TRUTH LIVES IN: branch arc/rubin-freq-bace; results on each host under pigauto_rubin/results/{bace,freq}; page https://claude.ai/artifact/WMHxQ8idafD5sLAVEJtg9Q
RESUME: Read .unlazy/rubin-freq-bace/GOAL.md, arcs.md, checkpoint.md; run script/rubin_status.sh; continue from NEXT.
HOST SPLIT (2026-09-24 19:5x nibi): BACE n300 seeds 1-102 on nibi (22634023, 204 tasks after scancel of the 198 PENDING
  seed 103-200 elements), seeds 103-200 on fir (61366948, 294 tasks, BLOCK 2). BACE n100 all seeds on nibi (22634022).
  freq: n100/300 seeds 1-50 on nibi root pigauto_rubin_freq (22633907/909), the rest on fir (61365818/19/20).
  A5 (n1000 BACE) will split seeds 2-100 across nibi/fir/rorqual/narval by disjoint ranges.
HOST SPLIT v2 (20:0x nibi; nibi BACE arrays still all PENDING after 2 h, fir/rorqual start in minutes):
  BACE n100: seeds 1-100 nibi (22634022, 120 tasks), 101-200 rorqual (21772063, 120 tasks).
  BACE n300: seeds 1-51 nibi (22634023, 102 tasks), 52-102 rorqual (21772064, 102 tasks), 103-200 fir (61366948, 294 tasks).
  rorqual validation 21771783 G-S4a PASS. narval validation FAILED: illegal instruction in ape (fir-built library not portable
  to narval CPUs) -> narval not used (rebuild would be needed). Totoro still over this user's 150 cap.
