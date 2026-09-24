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
