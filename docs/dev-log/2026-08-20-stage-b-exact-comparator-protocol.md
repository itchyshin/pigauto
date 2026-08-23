# Stage B: post-merge exact-conditional comparator protocol

## Status and purpose

This is a pre-run protocol.  PR #173 is merged to `main`; its
`predict_method = "exact"` route is therefore eligible for comparison, but
no Stage-B outcome has been run or claimed.  The study asks a bounded
competitiveness question: in the labelled regimes below, where does opt-in
exact conditional prediction sit relative to the package default and external
methods?  It cannot justify a default change, a parity claim, or a general
superiority claim.

## Locked data regimes and masks

Use the bundled labelled `avonet300` data with `tree300`, retaining only cells
that are originally observed before masking.  The two labelled regimes are
reported separately:

| ID | Regime | Traits scored | Why separate |
|---|---|---|---|
| B-continuous | Four continuous AVONET traits: `Mass`, `Beak.Length_Culmen`, `Tarsus.Length`, and `Wing.Length` | Every method with a continuous output | It is the fair common-method comparison. |
| B-mixed | AVONET's four continuous plus `Trophic.Level`, `Primary.Lifestyle`, and `Migration` | Each method only on traits it can validly return | It preserves pigauto's mixed-type scope rather than letting continuous-only methods erase it. |

For each regime, construct exactly five shared, independently seeded 30% MCAR
masks over originally observed cells.  The mask seeds are `20260901` through
`20260905`; within a seed, every method receives exactly the same masked input
and every scored method sees the same truth cells.  Persist the full mask,
original-observed indicator, trait map, source SHA, package versions, fit seed,
and wall-clock receipt before examining aggregate metrics.

## Methods

1. `pigauto::impute(predict_method = "per_column")`: production default.
2. `pigauto::impute(predict_method = "exact")`: opt-in exact conditional
   route; it must not change the default's configuration, data split, epochs,
   or stochastic seed.
3. `Rphylopars::phylopars(model = "BM")`: continuous traits only.
4. `phylolm::phylolm(model = "lambda")` intercept-only conditional predictor:
   continuous traits only.
5. `missForest::missForest()`: mixed-type, non-phylogenetic comparator.
6. `BACE::bace()`: phylogenetic chained-equation comparator where it can fit.

No unavailable method is silently substituted.  A missing dependency, invalid
trait type, convergence error, malformed output, or timeout is retained as a
method--trait--mask failure with the verbatim error and elapsed time.  A method
is marked **not applicable** only where its declared data type makes the
comparison invalid; that is distinct from a failed applicable fit.

## Estimands and reporting

For continuous traits report one row per method, trait, and mask: normalized
RMSE, raw-scale RMSE, Pearson correlation, runtime, and status.  For ordinal
and categorical traits report Brier score where probabilities are available,
plus accuracy and runtime; do not fabricate a continuous RMSE.  Aggregate only
within a trait and regime, displaying all five mask values, mean, SD, and
Monte-Carlo standard error (`SD / sqrt(5)`).  Show failures and unavailable
cells beside, not beneath, the performance tables.

The primary exact-solver contrast is paired, per trait and mask:

`loss(exact) - loss(default)`.

Negative loss contrasts favour the exact route, but are descriptive at five
masks; they are not a default-change threshold.  The external-method tables
describe a competitiveness boundary only.  In particular, neutral or unstable
`Mass` results remain an unresolved case and cannot be pooled away.

## Execution and verification gates

1. Before the campaign, exercise each method on one shared-mask smoke and
   verify that default and exact have identical input/mask receipts and
   differ only in `predict_method`.
2. Run the five masks as separate, resumable receipts; do not reuse a receipt
   unless source SHA, configuration, package versions, and mask hash match.
3. Preserve raw per-mask values, warnings, errors, runtimes, and the rendered
   summary.  A regenerated report must replace the stale pre-exact report;
   the latter is not evidence for this protocol.
4. A reviewer who did not build the summary checks the receipts, paired
   contrasts, failure table, and final wording before any public claim.

## Implementation status

`script/bench_external_continuous_core.R` implements the B-continuous
component; `script/bench_external_mixed_core.R` implements B-mixed.  Both use
the locked five masks and immutable per-mask RDS receipts.  The mixed runner
records continuous and discrete metrics separately and marks continuous-only
methods as **not applicable**, rather than failed, for the three discrete
traits.  Neither script is authorised to run the claim-bearing comparison
until a fresh timing smoke, cost estimate, and compute approval are recorded.

## Explicit non-claims

Stage B cannot establish equality, general parity, universal superiority,
mixed-type dominance by a continuous-only competitor, an exact-route default
change, interval calibration, or a CRAN readiness claim.  It also makes no
claim about the source of any difference (e.g., solver versus GNN) beyond the
locked `predict_method` contrast.
