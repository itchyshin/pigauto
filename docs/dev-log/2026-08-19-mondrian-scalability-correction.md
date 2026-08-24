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

The task request was then made explicitly typed as `--gres=gpu:h100:1`.
Local shell syntax checking and `sbatch --test-only` both passed.  Tamia job
`419946` was the resulting same-scope rerun.

Job `419946` completed all four `split` tasks successfully in 54 seconds;
each method receipt has `status = "ok"` and each task confirmed
`cuda_available=TRUE`.  Its runner elapsed times were 5.83, 4.84, 5.88, and
8.17 seconds for 300, 600, 1,000, and 1,500 tips respectively.  Slurm
accounted one H100 per task, but the node-level GPU accounting trace places
all four CUDA processes on physical GPU 0.  This is valid evidence that the
CUDA route and the four inputs run, but it is **not** valid four-GPU scaling
evidence and cannot be compared with the CPU timing receipt.

Tamia job `419948` added `--gpu-bind=single:1`, but repeated the same
structural mistake in a subtler form: it started four *separate* one-rank
`srun` steps.  Each step has local rank zero, so Slurm validly exposed a
single logical device numbered `0` to each task but bound all four processes
to physical GPU 0.  It completed in 54 seconds, all four task receipts have
`status = "ok"`, and CUDA was available, but node-level accounting records
all four processes against PCI bus `00000000:4E:00.0`.  Thus job `419948` is
also only a CUDA/input smoke—not a valid four-GPU scaling ladder.

The repair launched **one four-rank `srun` step** with
`--gpus-per-task=h100:1 --gpu-bind=single:1`; rank selects the nested input
size.  The next attempt will record Slurm rank plus GPU UUID and PCI bus ID.
Because Slurm may renumber a task's isolated device as
`CUDA_VISIBLE_DEVICES=0`, acceptance is **four distinct UUID/PCI-bus values**,
not four distinct local logical-device numbers.

Tamia job `424950` completed this corrected ladder in 24 seconds. All four
retained `split.rds` receipts report `status = "ok"` and `cuda_available=TRUE`;
their measured runner elapsed times were 7.542, 7.118, 7.030, and 9.650
seconds for 300, 600, 1,000, and 1,500 tips, respectively. The four tasks
recorded four distinct physical GPU UUID/PCI pairs (PCI buses `4E:00.0`,
`CB:00.0`, `5F:00.0`, and `DB:00.0`). It is therefore valid four-device
**operational** scaling evidence for these bounded split inputs. It remains
neither a CPU-versus-GPU comparison nor evidence about interval calibration,
and it does not authorise a full multi-mask campaign.
