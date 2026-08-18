# Stage B: post-merge exact-solver comparator protocol

Blocked until PR #173 is merged and its CI is green.  Pre-register two labelled
datasets/regimes, five shared masks each, and preserve per-mask values, elapsed
time, unavailable fits, and errors.  Compare pigauto default, opt-in
`predict_method = "exact"`, Rphylopars, phylolm, BACE, and missForest.

Continuous traits and mixed-type traits are reported separately.  The outcome
is a competitiveness boundary only: it cannot establish parity, trigger a
default change, or erase pigauto's mixed-type scope.  If Mass remains neutral,
it is retained as an unresolved case.
