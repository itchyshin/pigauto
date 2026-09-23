# MI under phylogenetic GLS: run log

ESTIMATE sweep (2026-09-23): 16 regimes x 120 reps = 1,920 cells; ~918 s per cell from the local smoke (GNN methods scaled 150 -> 500 epochs); ~490 core-hours; ~5 h wall at ~100 concurrent fir tasks. Approved by Shinichi 2026-09-23 ("Full sweep on fir").
LAUNCH fir: source 67a8df8 on /scratch/snakagaw/pigauto-mi-gls; install job then a 2-task smoke (task 1 = regime 1 rep 1, cheapest; task 1681 = regime 15 rep 1, n = 1000 MAR_phylo both-missing) before the full array.
MEASURED fir smoke (500 epochs, m = 20, 1 core): regime 1 rep 1 in 4 min 27 s (MaxRSS 1.8 GB); regime 15 rep 1 (n = 1000, MAR_phylo, both missing) in 8 min 47 s (MaxRSS 6.5 GB). Both receipts written; draw_cond and oracle PGLS slopes 0.643 and 0.660 in regime 1.
ESTIMATE revised: ~6.5 min per cell average -> ~210 core-hours (was 490 from the local 150-epoch smoke). Full array launched for tasks 2-1680 and 1682-1920.
