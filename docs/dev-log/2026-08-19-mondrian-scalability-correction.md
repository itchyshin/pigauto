# Stage C: FishBase scalability correction

## What happened

The full FishBase timing receipt used the canonical 10,484-tip, six-trait
input, one fixed observed-only 20% mask, `seed = 20260818`, ten epochs, and
one compute thread.  It began the `split` fit but did not finish or write a
`split.rds` receipt within 24 hours and 19 minutes.  The process was CPU-bound
(rather than waiting for I/O or memory), peaked at about 16 GB RSS, and was
terminated by the approved operator request.  The only retained output is the
complete `mask_receipt.rds`.

This is an operational scalability result, not a failed scientific result and
not evidence about interval calibration.  In particular, no method comparison
can be reconstructed from the stopped fit.

## Correction now implemented

The masked-confirmation runner now:

1. validates and reuses a pre-existing `mask_receipt.rds` for the same input,
   tree labels, data set, and seed;
2. accepts a comma-separated method subset (`split`, `mondrian`); and
3. never overwrites an existing completed method receipt.

A local 30-tip smoke ran `split`, resumed with `mondrian`, then invoked
`split` again.  It verified that the fixed mask was reused and both receipts
survived.  This protects completed *method* receipts.  It does not checkpoint
inside `impute()`, so an interrupted individual fit still needs a fresh fit.

## Revised compute gate

Before any renewed full FishBase campaign, run a separately retained,
one-method CPU-versus-Tamia-GPU scaling ladder on deterministic FishBase
subsets. Tamia's H100 policy allocates a whole four-GPU node, so the immediate
bounded feasibility receipt is four independent `split` fits (300, 600, 1,000,
and 1,500 tips) dispatched one per GPU, with a hard one-hour Slurm limit. Its
purpose is to measure device use and identify whether the graph/GNN portion
materially changes wall time. It is not a calibration run.

Only after that receipt reports elapsed time, peak memory, and a credible
CPU/GPU comparison may a revised multi-mask design be costed.  Production
would use resumable per-method/per-mask tasks on persistent DRAC project
storage with explicit time and memory limits.  It remains separately approval
gated and cannot make a Mondrian calibration claim merely because the
feasibility ladder completes.

## GPU ladder execution record

The first submitted ladder attempt (Tamia job `419940`) acquired a full H100
node but failed after three seconds, before any R task started or any
per-size receipt was written.  The outer allocation used typed
`h100:4` resources whereas each nested `srun` requested an untyped GPU;
Tamia Slurm rejected that inconsistent GRES request.  `nvidia-smi` confirms
that all four H100s were available in the allocation, but this failed attempt
is not a timing, CUDA, or model result.

The task request is now explicitly typed as `--gres=gpu:h100:1`.  Local shell
syntax checking and `sbatch --test-only` both passed.  The same bounded
four-task, one-hour ladder was resubmitted as Tamia job `419946` after the
approved repair.

Job `419946` completed all four `split` tasks successfully in 54 seconds;
each method receipt has `status = "ok"` and each task confirmed
`cuda_available=TRUE`.  Its runner elapsed times were 5.83, 4.84, 5.88, and
8.17 seconds for 300, 600, 1,000, and 1,500 tips respectively.  Slurm
accounted one H100 per task, but the node-level GPU accounting trace places
all four CUDA processes on physical GPU 0.  This is valid evidence that the
CUDA route and the four inputs run, but it is **not** valid four-GPU scaling
evidence and cannot be compared with the CPU timing receipt.

The next same-scope repair adds `--gpu-bind=single:1` and records each task's
`CUDA_VISIBLE_DEVICES` value.  A successful rerun must show four distinct
visible-device assignments before its elapsed times are interpreted as a
four-GPU feasibility ladder.  It remains a predeclared operational receipt,
not a calibration result.
