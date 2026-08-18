# Stage A: active-recovery implementation and smoke receipt

This is an implementation receipt, not a performance result.

On source SHA `f925e17`, a local 30-tip, five-epoch structural smoke completed
for both continuous and binary outcomes. For each family it fit the initial
data, selected active/random/uncertainty candidates from the same candidate
set, restored exactly one true value, refit each policy with the same fitting
seed, and produced finite predictions on the protected fixed test set. Each
three-policy smoke took about 11 seconds on the local machine.

The 30-tip smoke emits the expected small-validation warning (two validation
cells); its purpose is treatment integrity and execution viability only. It
does not estimate effect size, interval calibration, failure rate, or support
the active-imputation headline. The next permitted empirical step is the
pre-registered 100-replicate-per-cell pilot on Totoro, followed by a measured
wall-time/variance estimate and explicit approval before any full campaign.
