Hi @Santiago-0rtega,

Thanks for the careful report and full reproducer. I reproduced all three examples, and you were right: the vignette workflow was scoring predictions at cells whose true values were still available. Those perfect scores were misleading. Sorry we missed this.

I've prepared and tested a fix that adds `predict(fit, evaluation = TRUE)`, blocks unsafe scoring of production prediction objects, and corrects both vignettes. The corrected results match yours. I also checked the GNN-on prediction path with an active GNN contribution, and the full test suite now passes with no failures (3,672 passed, 184 warnings, eight skips).

The fix is local for now and has not been merged, so it is not yet on CRAN or the website. Until it is available, `evaluate(fit)` remains the workaround.

Thanks also for being clear about the AI assistance and what you verified yourself. Your reproducible examples were very helpful.

Best,
Shinichi
