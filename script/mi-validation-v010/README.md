# pigauto v0.10.0 multiple-imputation validation

This harness tests downstream fixed-effect inference under known data-generating
processes. It does **not** validate random-effect variances, correlations, BLUPs,
latent loadings, MNAR missingness, missing discrete covariates, or tree
uncertainty.

## Design

The deterministic manifest uses base seed `20260710` and crosses:

| DGP | Sample size | Analysis | Truth |
|---|---:|---|---|
| Gaussian | 250 species | `lm(y ~ x + z)` | beta-x = 0.60, beta-z = -0.40 |
| Binary | 400 species | `glm(y ~ x + z, binomial())` | beta-x = 0.80, beta-z = -0.50 |
| Mixed | 100 species x 4 observations | `lmer(y ~ x + z + (1 | species))` | beta-x = 0.60, beta-z = -0.40, tau-squared = 0.64, sigma-squared = 1 |

Each DGP has a phylogeny-dominant regime and an auxiliary nonlinear regime.
Exactly 30% of `x` is masked using probabilities determined by fully observed
`y` and `z` (MAR). The outcome is included in pigauto's trait block and `z` is
passed as a fully observed covariate.

Each replicate trains pigauto once with `M = 50` by default. The fitted object
then supplies both paired draw sets:

- `mc_dropout`: Brownian posterior draws blended with stochastic GNN forward
  passes;
- `conformal`: conformal-width Normal draws regenerated from the same fit,
  calibrated gates, mask, and point prediction.
- `pmm` (redesign candidate): stochastic Brownian/MC-dropout predictions
  followed by predictive mean matching to one of the five nearest observed
  donors. This is evaluated in the harness before any public API promotion.
- `oracle_conditional` (positive control): draws from the known DGP-specific
  conditional distribution of missing `x`, including the true random intercept
  in the mixed DGP.
- `standard_smc` (positive control): `smcfcs` for `lm` and logit `glm`, and
  `jomo::jomo.smc(model = "lmer")` for the Gaussian mixed model. The fully
  observed `z^2` term is included as an auxiliary variable.

This pairing isolates the draw rule from model-training variability and avoids
retraining for each candidate. Complete-data oracle and complete-case analyses
are retained as comparators. For the mixed DGP, variance-component estimates and
singularity/boundary rates are diagnostics only; they are never Rubin-pooled.

## Commands

Run commands from the repository root. Generated manifests, task results, and
summaries are `.rds` files under `script/mi-validation-v010/results/` and are
ignored by git.

```sh
# Pure DGP/manifest smoke; works without torch or libtorch.
Rscript script/mi-validation-v010/0_prepare.R --profile=smoke
Rscript script/mi-validation-v010/1_run.R --profile=smoke --all-pending --dry-run
Rscript script/mi-validation-v010/2_summarise.R --profile=smoke

# One real, two-draw/two-epoch smoke when libtorch is available.
Rscript script/mi-validation-v010/1_run.R --profile=smoke --task=1 --force

# Economical positive-control gate: 50 replicates per cell, no pigauto training.
Rscript script/mi-validation-v010/0_prepare.R --profile=pilot --reps=50 \
  --output=script/mi-validation-v010/results/controls-pilot
seq 1 300 | OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 xargs -P 96 -I{} \
  Rscript script/mi-validation-v010/1_run.R --profile=pilot --reps=50 \
    --output=script/mi-validation-v010/results/controls-pilot \
    --controls-only --task={}
Rscript script/mi-validation-v010/2_summarise.R --profile=pilot --reps=50 \
  --output=script/mi-validation-v010/results/controls-pilot

# Ten-replicate pilot: 60 tasks, M=50, 500 epochs.
Rscript script/mi-validation-v010/0_prepare.R --profile=pilot
seq 1 60 | OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 xargs -P 24 -I{} \
  Rscript script/mi-validation-v010/1_run.R --profile=pilot --task={}
Rscript script/mi-validation-v010/2_summarise.R --profile=pilot

# Full campaign: 3,000 trained models, 500 epochs, and paired method results.
Rscript script/mi-validation-v010/0_prepare.R --profile=full
seq 1 3000 | OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 xargs -P 96 -I{} \
  Rscript script/mi-validation-v010/1_run.R --profile=full --task={}
Rscript script/mi-validation-v010/2_summarise.R --profile=full
```

Totoro is capped at 96 workers. Use a DRAC array only if the ten-replicate
pilot projects more than 12 hours or Totoro is unavailable. On DRAC, set
`SLURM_ARRAY_TASK_ID`; `1_run.R` reads it automatically. Never run the campaign
on a login node.

The CLI options have matching environment variables:

| CLI | Environment | Default by profile |
|---|---|---|
| `--profile=` | `PIGAUTO_MI_PROFILE` | `full` |
| `--reps=` | `PIGAUTO_MI_REPS` | 1 / 10 / 500 |
| `--m=` | `PIGAUTO_MI_M` | 2 / 50 / 50 |
| `--epochs=` | `PIGAUTO_MI_EPOCHS` | 2 / 500 / 500 |
| `--task=` | `PIGAUTO_MI_TASK_ID` | no implicit task |
| `--output=` | `PIGAUTO_MI_OUTPUT` | profile-specific results directory |

`--controls-only` (or `PIGAUTO_MI_CONTROLS_ONLY=true`) runs only the oracle and
standard substantive-compatible controls. It deliberately skips torch and
pigauto training, making the gate-attainability check much cheaper than another
candidate campaign. It requires the optional harness dependencies `smcfcs` and
`jomo`; neither is added to the package DESCRIPTION.

`--force` overwrites an existing task atomically. Otherwise completed task files
are skipped, so interrupted campaigns are resumable. `--all-pending` is intended
only for bounded sequential local runs.

The pilot is also the epoch-budget gate. Stay at 500 epochs for the full campaign
unless more than 5% of successful pilot fits both (a) reach epoch 500 without
early stopping and (b) improve validation loss by more than 1% across their last
three evaluation intervals (150 epochs at the default `eval_every = 50`). The
summariser reports this as `training$escalate_above_500`. If triggered, run a
second ten-replicate pilot with the proposed `--epochs` override before changing
the full campaign; do not jump directly to 2,000 epochs.

Project wall time from the measured pilot rather than an assumed per-fit cost:
`median(elapsed_seconds) * 3000 / workers`, with 20% scheduling/I/O headroom.
Use 24 workers for the pilot and at most 96 on Totoro. Move to DRAC only if that
projection exceeds 12 hours or Totoro is unavailable.

## Pre-specified fixed-effect gates

Every `x` and `z` cell for every DGP, regime, and method must have:

- at least 500 processed replicates and at least 95% successful pooled analyses;
- at least 95% successful downstream fits within every counted pooled analysis;
- absolute bias no larger than 0.10 empirical standard deviations;
- mean pooled SE / empirical SD between 0.90 and 1.10;
- 95% interval coverage between 92.5% and 97.5%;
- coverage Monte Carlo SE no larger than one percentage point;
- at least 99% finite estimates, SEs, intervals, and valid Rubin df/FMI/RIV.

The summariser refuses to issue a release pass from smoke or pilot evidence. It
reports `pilot_criteria_pass` without the 500-replicate completeness and 0.01
coverage-MCSE requirements. The oracle and standard SMC controls must first pass
that directional 50-replicate screen across all cells. Only then is a public
analysis-aware backend worth implementing and validating at 500 replicates per
cell. Current conformal, MC-dropout, and PMM methods remain experimental unless
they independently pass the full gate.

Variance-component behavior is reported descriptively. A method is flagged when
it adds more than 0.20 absolute relative bias over the oracle or increases the
boundary rate by more than five percentage points. A flag does not create a
variance-component pooling claim.
