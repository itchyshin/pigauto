# Two-worker gllvmTMB recovery, 2026-10-06

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

The first registered gllvmTMB seed had taken 1,275.364 seconds alone. These
two concurrent seeds supported a four-wave estimate of about 87 minutes for
the remaining seven seeds. This pre-run was a mechanism and resource check.

The remaining seeds, 2026100604 to 2026100610, completed in four supervised
waves: seeds 04/05 in 1,332 seconds, 06/07 in 1,298 seconds, 08/09 in 1,321
seconds, and seed 10 in 1,285 seconds. All seven supervisors exited zero. The
continuation took 5,236 seconds (87.27 minutes), within its three-hour bound.
All ten registered seeds passed posterior convergence, 20 genuine downstream
fits per seed, independent coefficient/SE and Rubin checks, and saved-fit
reload. The locked `gllvm-summary.csv` passed the prespecified ±0.15 mean-bias
margin for all eight fixed effects. Coverage was 9/10 for two terms and 10/10
for six terms; ten seeds do not certify nominal 95% coverage.

`full-remote-sha256.txt` records ten retained raw-RDS digests and 21 copied
CSV/log digests. All 21 local copies matched. The raw RDS files remain on
Totoro. The summary is a bounded result for the audit source and is not
exact-tarball release evidence.
