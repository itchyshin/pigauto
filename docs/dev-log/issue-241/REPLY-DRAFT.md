Thanks, Santiago, for the careful report and the full reproducer. I reproduced all three examples locally. You are right: the vignette workflow scored production predictions at cells whose true values were still available, so those perfect scores were misleading. Sorry we missed this when separating production imputation from evaluation.

I have prepared a fix that adds `predict(fit, evaluation = TRUE)`, using the held-out baseline and masking validation/test cells from the prediction inputs. The scorer now rejects production prediction objects when they are used for held-out scoring. Both vignette workflows have been corrected, including type-aware splitting and scoring in getting-started.

The corrected mixed-types example gives test RMSE 0.9285 for mass, binary accuracy 0.50, and categorical accuracy 0.3333, matching your masked results. Regression tests also cover the GNN-on path with an active GNN contribution. The fix is local so far; it is not yet in a released package or on the website. For the current release, `evaluate(fit)` remains the supported workaround.

Thank you also for explaining how AI was used and what you checked yourself. The reproducible evidence made this straightforward to verify.
