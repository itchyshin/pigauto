# Mondrian real-data confirmation: run log

One line per event, oldest first. `ESTIMATE` lines are written before the run they
describe; `MEASURED` lines after. Times are wall clock on the named host.

ESTIMATE pantheria pre-run: one split fit, n = 4,027, 500 epochs, 1 thread, Totoro. No measured basis at 500 epochs (08-18 timing run epochs unrecorded). Stop and re-report if it exceeds 3 h.
BLOCKED tamia 2026-09-23: ControlMaster socket ~/.ssh/cm-snakagaw@tamia.alliancecan.ca:22 absent; BatchMode login denied (keyboard-interactive). FishBase pre-run waits for Shinichi to open one Tamia session. Flagged once (D-64).
LANE NOTE 2026-09-23: a lambda-implementation lane runs in another Claude account (baseline files). Overlap only possible at R/fit_pigauto.R if the default flips; diff their branch first.
ESTIMATE mi-sim pre-run: 5 reps, n = 1000, 500 epochs, m = 20, 2 fits per rep (split, mondrian) + reference, 1 thread each, 5 in parallel, Totoro, source bfecd84. Basis: 08-16 mech cells 10-20 min per fit at n = 1000. Expect 20-40 min wall; stop and re-report past 60 min.
FAILED-SMOKE mi-sim pre-run 05:54: R --vanilla dropped ~/R/lib (torch); all 5 reps died at load; logs kept in results/mi-sim-prerun/failed-libpath. Fix: R_LIBS=lib_<sha>:$HOME/R/lib. Relaunched 05:54:32; 5 procs at ~97% CPU at 45 s. The same R_LIBS fix applies to 10_launch_totoro.sh.
KOHAKU 2026-09-23: Shinichi asked to move FishBase GPU to kohaku. Checked: no R/Rscript (needs root via John, per tools/kohaku-setup.md); only per-user Julia in ~/kohaku_work/gpu_env; no sudo; neighbour llama-server holds 50.5 GB of 97.9 GB VRAM. Not installed anything. FishBase stays on hold pending Shinichi's choice (kohaku user-space R attempt vs Tamia vs ask John).
FISHBASE TARGET 2026-09-23: Shinichi chose Tamia over kohaku. Pre-run on Tamia once the instrumented build exists; full campaign still waits for Q1.
ESTIMATE pantheria pre-run confirmed as the campaign cell pantheria-mcar-m20260818, mondrian method only first (split resumes after). Launch 2026-09-23.
