# Two-worker gllvmTMB feasibility check, 2026-10-06

This bounded pre-run used two independent registered seeds, 2026100602 and
2026100603, concurrently on Totoro. Each seed used n = 300, four genuine
Gaussian outcomes, 60 masked cells per outcome, M = 20 posterior imputations,
and four serial sampler chains. The two processes each used one CPU core with
`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`.
Observed resident memory was about 440 MB per process near completion.

Both supervisors exited zero. Seed 2026100602 took 1,307.395 seconds and seed
2026100603 took 1,301.696 seconds; the concurrent supervisor took 1,309 seconds
(21.82 minutes). Each seed passed posterior convergence, 20 genuine downstream
fits, independent coefficient and standard-error oracles, independent Rubin
arithmetic, and pooling after saved-fit reload. The logs and pooled tables are
beside this note. `remote-sha256.txt` records the four copied CSV/log digests
and two retained raw-RDS digests; all four local copies matched their remote
digests. The raw RDS files remain at
`/home/snakagaw/pigauto_cran_011_audit/recovery-parallel/` on Totoro.

The remote driver, adapter source, and DESCRIPTION SHA-256 digests matched this
branch's `script/cran-0.11-recovery/run.R`, `R/pool_mi.R`, and `DESCRIPTION`.
The installed source is the audit's 0.11.0 adapter snapshot, not a frozen
post-merge release tarball.

The first registered gllvmTMB seed had taken 1,275.364 seconds alone. The two
concurrent seeds therefore support a two-worker planning estimate of five
waves for the original nine remaining seeds: about 109 minutes of wall time.
Seven seeds remain after this pre-run, requiring four more waves, about 87
minutes at the measured pace. Automatic sampler extension or later resource
contention could lengthen that estimate. This pre-run is a mechanism and
resource check, not a ten-seed bias or coverage result. No further gllvmTMB
seed was launched here.
