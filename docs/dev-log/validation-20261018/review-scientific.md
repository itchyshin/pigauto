# Independent scientific review of the adapters

Reviewer: separate Sol agent, high effort; reviewer did not write runner adapters. Read-only review; no fits run.

## Findings and fixes

Analysis initially duplicated formula variables as auxiliary variables. It now uses the supported narrow API and checks actual downstream fit failures. Complete-data references are retained. Finite-sample t intervals are used for lm and GLS.

Mixed repeated-observation scenarios now apply their named lambda before generating traits and mask within-species observations separately by trait. Family tests use the existing three-trait ZI and compositional generators, with gate, count, magnitude, Aitchison, CLR and simplex metrics. GNN-off is explicit for these family checks.

Non-finite masked predictions cannot disappear through na.rm scoring. Prediction coverage is emitted even when absent, as an explicit failure. Posterior inference requires sampler convergence and all requested downstream fits; actual DGP and sampler seeds are retained beside the manifest identity seed.

## Verdict

READY for the bounded pilot on final read. Production scientific claims remain pending the measured campaign. This does NOT cover scientific calibration or successful completion of the full programme.
