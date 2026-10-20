pkgname <- "pigauto"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "pigauto-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('pigauto')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("build_phylo_graph")
### * build_phylo_graph

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: build_phylo_graph
### Title: Build a phylogenetic graph representation from a tree
### Aliases: build_phylo_graph

### ** Examples

set.seed(1)
tree <- ape::rtree(30)
g <- build_phylo_graph(tree, k_eigen = 4)
dim(g$adj)    # 30 x 30
dim(g$coords) # 30 x 4
dim(g$D)      # 30 x 30



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("build_phylo_graph", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calibration_df")
### * calibration_df

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calibration_df
### Title: Compute calibration data for probability predictions
### Aliases: calibration_df

### ** Examples

calibration_df(
  truth = c(0, 0, 1, 1),
  prob = c(0.1, 0.3, 0.7, 0.9),
  n_bins = 2
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calibration_df", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("check_pigauto")
### * check_pigauto

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: check_pigauto
### Title: Check whether traits and a tree are ready for pigauto
### Aliases: check_pigauto

### ** Examples

## No test: 
tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
traits <- data.frame(x = c(1, NA, 3), row.names = tree$tip.label)
check_pigauto(traits, tree)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("check_pigauto", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compare_methods")
### * compare_methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compare_methods
### Title: Compare BM baseline and pigauto methods across replicates
### Aliases: compare_methods

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
torch_ready <- requireNamespace("torch", quietly = TRUE) &&
  isTRUE(tryCatch(torch::torch_is_installed(), error = function(e) FALSE)) &&
  isTRUE(tryCatch({
    probe <- torch::torch_tensor(1)
    as.numeric(probe$item()) == 1
  }, error = function(e) FALSE))
if (torch_ready) {
  cmp <- compare_methods(preprocess_traits(traits, tree), tree, seeds = 1L,
                         epochs = 5L, verbose = FALSE, gnn = TRUE)
  aggregate(value ~ method + trait + metric, data = cmp, FUN = mean)
} else {
  message("Skipping GNN comparison: install torch and its Lantern runtime.")
}
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compare_methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("completed_data")
### * completed_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: completed_data
### Title: Extract completed trait data
### Aliases: completed_data

### ** Examples

## Not run: completed_data(result)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("completed_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("confusion_matrix")
### * confusion_matrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: confusion_matrix
### Title: Compute a confusion matrix for categorical or binary predictions
### Aliases: confusion_matrix

### ** Examples

confusion_matrix(
  truth = c("a", "a", "b", "b"),
  predicted = c("a", "b", "b", "b")
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("confusion_matrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cross_validate")
### * cross_validate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cross_validate
### Title: k-fold cross-validation for pigauto trait imputation
### Aliases: cross_validate

### ** Examples

## No test: 
data(avonet300, tree300, package = "pigauto")
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(12L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
cv <- cross_validate(preprocess_traits(traits, tree), tree, k = 3L,
                     seeds = 1L, epochs = 1L, verbose = FALSE,
                     k_eigen = 4L, hidden_dim = 8L,
                     use_transformer_blocks = FALSE)
print(cv)
summary(cv)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cross_validate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("evaluate")
### * evaluate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: evaluate
### Title: Evaluate a fitted pigauto model on its test set
### Aliases: evaluate

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
pd <- preprocess_traits(traits, tree)
splits <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map)
fit <- fit_pigauto(pd, tree, splits = splits, epochs = 5L,
                   verbose = FALSE)
eval_df <- evaluate(fit)
eval_df[eval_df$metric == "rmse", ]
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("evaluate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("evaluate_imputation")
### * evaluate_imputation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: evaluate_imputation
### Title: Evaluate imputation performance against known values
### Aliases: evaluate_imputation

### ** Examples

## No test: 
truth <- matrix(c(1, 2, 3, 4), nrow = 2L)
splits <- list(test_idx = 1L, val_idx = 2L, n = 2L, p = 2L)
evaluate_imputation(truth, truth, splits)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("evaluate_imputation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_baseline")
### * fit_baseline

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_baseline
### Title: Fit the phylogenetic baseline
### Aliases: fit_baseline

### ** Examples

## No test: 
data(avonet300, tree300, package = "pigauto")
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
pd     <- preprocess_traits(traits, tree)
splits <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map)
bl     <- fit_baseline(pd, tree, splits)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_baseline", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_pigauto")
### * fit_pigauto

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_pigauto
### Title: Fit a pigauto model for trait imputation
### Aliases: fit_pigauto

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
data <- preprocess_traits(traits, tree)
splits <- make_missing_splits(data$X_scaled, trait_map = data$trait_map,
                              seed = 1)
fit <- fit_pigauto(data, tree, splits = splits, epochs = 5L,
                   verbose = FALSE, seed = 1)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_pigauto", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("impute")
### * impute

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: impute
### Title: Impute missing phylogenetic traits (convenience wrapper)
### Aliases: impute

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
df <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                c("Mass", "Wing.Length"), drop = FALSE]
rownames(df) <- tree$tip.label
df$Mass[seq_len(3L)] <- NA_real_
result <- impute(df, tree, epochs = 5L, verbose = FALSE)
head(result$completed)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("impute", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("load_pigauto")
### * load_pigauto

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: load_pigauto
### Title: Load a saved pigauto model
### Aliases: load_pigauto

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
traits$Mass[seq_len(3L)] <- NA_real_
fit <- impute(traits, tree, epochs = 5L, verbose = FALSE)$fit
path <- tempfile(fileext = ".pigauto")
save_pigauto(fit, path)
pred <- predict(load_pigauto(path))
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("load_pigauto", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("make_missing_splits")
### * make_missing_splits

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: make_missing_splits
### Title: Split cells into train/val/test for imputation evaluation
### Aliases: make_missing_splits

### ** Examples

X <- matrix(rnorm(100), nrow = 20)
splits <- make_missing_splits(X, missing_frac = 0.25, seed = 1)
length(splits$val_idx)

# MAR: missingness depends on another trait
splits_mar <- make_missing_splits(X, mechanism = "MAR_trait",
  mechanism_args = list(driver_col = 1, beta = 2))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("make_missing_splits", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mask_missing")
### * mask_missing

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mask_missing
### Title: Create an observed/missing mask matrix
### Aliases: mask_missing

### ** Examples

X <- matrix(c(1, NA, 3, 4), nrow = 2)
mask_missing(X)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mask_missing", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("multi_impute")
### * multi_impute

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: multi_impute
### Title: Generate experimental stochastic completion datasets
### Aliases: multi_impute

### ** Examples

## No test: 
library(pigauto)
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
df <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                c("Mass", "Wing.Length"), drop = FALSE]
rownames(df) <- tree$tip.label
df$Mass[seq_len(3L)] <- NA_real_
mi <- multi_impute(df, tree, m = 2L, draws_method = "conformal",
                   epochs = 5L, verbose = FALSE)
print(mi)
lapply(mi$datasets, head)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("multi_impute", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("multi_impute_analysis")
### * multi_impute_analysis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: multi_impute_analysis
### Title: Analysis-aware multiple imputation for narrow regression models
### Aliases: multi_impute_analysis

### ** Examples

set.seed(11)
dat <- data.frame(
  y = stats::rnorm(80),
  x = stats::rnorm(80),
  z = stats::rnorm(80)
)
dat$x[seq(4, 80, by = 5)] <- NA_real_
mi <- multi_impute_analysis(
  dat, y ~ x + z, missing = "x", model = "lm", m = 5, seed = 12
)
fits <- with_imputations(
  mi, function(d) stats::lm(y ~ x + z, data = d), .progress = FALSE
)
pool_mi(fits)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("multi_impute_analysis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("multi_impute_trees")
### * multi_impute_trees

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: multi_impute_trees
### Title: Posterior-tree prediction sensitivity
### Aliases: multi_impute_trees

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
trees <- structure(list(tree, tree), class = "multiPhylo")
df <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                c("Mass", "Wing.Length"), drop = FALSE]
rownames(df) <- tree$tip.label
df$Mass[seq_len(3L)] <- NA_real_
mi <- multi_impute_trees(df, trees, m_per_tree = 1L, epochs = 5L,
                         verbose = FALSE)
print(mi)
mass_by_tree <- vapply(mi$datasets, function(dat) dat$Mass,
                       numeric(nrow(df)))
apply(mass_by_tree, 1L, stats::sd)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("multi_impute_trees", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pigauto_report")
### * pigauto_report

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pigauto_report
### Title: Generate an HTML imputation report from a pigauto fit
### Aliases: pigauto_report

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
traits$Mass[seq_len(3L)] <- NA_real_
result <- impute(traits, tree, epochs = 5L, verbose = FALSE)
pigauto_report(result, output_path = tempfile(fileext = ".html"),
               open = FALSE)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pigauto_report", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.pigauto_benchmark")
### * plot.pigauto_benchmark

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.pigauto_benchmark
### Title: Plot a pigauto benchmark
### Aliases: plot.pigauto_benchmark

### ** Examples

benchmark <- structure(list(summary = data.frame(
  scenario = rep("BM", 2), trait = rep("trait1", 2),
  method = c("baseline", "pigauto"), metric = rep("rmse", 2),
  mean = c(1, 0.8)
)), class = "pigauto_benchmark")
plot(benchmark)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.pigauto_benchmark", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.pigauto_fit")
### * plot.pigauto_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.pigauto_fit
### Title: Plot diagnostics for a fitted pigauto model
### Aliases: plot.pigauto_fit

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
data <- preprocess_traits(traits, tree)
splits <- make_missing_splits(data$X_scaled, trait_map = data$trait_map)
# The package default disables the optional GNN, so this example works
# without torch or its Lantern runtime.
fit <- fit_pigauto(data, tree, splits = splits, epochs = 5L,
                   eval_every = 1L, verbose = FALSE, gnn = FALSE)
plot(fit)
plot(fit, type = "gates")
plot(fit, type = "conformal")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.pigauto_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.pigauto_pred")
### * plot.pigauto_pred

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.pigauto_pred
### Title: Plot predictions from a pigauto model
### Aliases: plot.pigauto_pred

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
traits$Mass[seq_len(3L)] <- NA_real_
pred <- impute(traits, tree, epochs = 5L, verbose = FALSE)$prediction
plot(pred)
plot(pred, type = "scatter", data = traits)
plot(pred, type = "intervals", trait = "Mass")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.pigauto_pred", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_comparison")
### * plot_comparison

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_comparison
### Title: Forest-plot style comparison of benchmark results
### Aliases: plot_comparison

### ** Examples

## No test: 
results <- data.frame(
  trait = rep("Mass", 2L), type = rep("continuous", 2L),
  metric = rep("rmse", 2L), method = c("baseline", "pigauto"),
  value = c(1, 0.8)
)
plot_comparison(results)
plot_comparison(results, metric = "rmse")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_comparison", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_history_gg")
### * plot_history_gg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_history_gg
### Title: Plot training history (ggplot2, deprecated)
### Aliases: plot_history_gg

### ** Examples

fit <- structure(
  list(history = data.frame(epoch = 1:3, val_loss = c(1, 0.7, 0.5))),
  class = "pigauto_fit"
)
plot_history_gg(fit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_history_gg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_uncertainty")
### * plot_uncertainty

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_uncertainty
### Title: Plot uncertainty ribbons for imputed trait values
### Aliases: plot_uncertainty

### ** Examples

pred <- list(
  imputed = data.frame(mass = c(10, 12, 15)),
  se = matrix(c(1, 1.5, 2), ncol = 1,
              dimnames = list(NULL, "mass"))
)
plot_uncertainty(pred, trait_name = "mass")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_uncertainty", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pool_mi")
### * pool_mi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pool_mi
### Title: Pool downstream model fits across multiple imputations (Rubin's
###   rules)
### Aliases: pool_mi

### ** Examples

## No test: 
dat <- data.frame(y = stats::rnorm(30L), x = stats::rnorm(30L),
                  z = stats::rnorm(30L))
dat$x[seq(3L, 30L, by = 5L)] <- NA_real_
mi <- multi_impute_analysis(
  data = dat, formula = y ~ x + z, missing = "x",
  model = "lm", m = 2L
)
fits <- with_imputations(mi, function(d) lm(y ~ x + z, data = d))
pool_mi(fits)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pool_mi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.pigauto_fit")
### * predict.pigauto_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.pigauto_fit
### Title: Impute missing traits using a fitted pigauto model
### Aliases: predict.pigauto_fit

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
data <- preprocess_traits(traits, tree)
splits <- make_missing_splits(data$X_scaled, trait_map = data$trait_map)
fit <- fit_pigauto(data, tree, splits = splits, epochs = 5L,
                   verbose = FALSE)
pred <- predict(fit, return_se = TRUE)
pred$imputed
pred$se
pred$probabilities
pred2 <- predict(fit, n_imputations = 2L)
pred2$imputed_datasets
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.pigauto_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("preprocess_traits")
### * preprocess_traits

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: preprocess_traits
### Title: Preprocess trait data: align to tree, encode into latent space
### Aliases: preprocess_traits

### ** Examples

# Single-obs per species (backward compatible)
data(avonet300, tree300, package = "pigauto")
traits <- avonet300
rownames(traits) <- traits$Species_Key
traits$Species_Key <- NULL
pd <- preprocess_traits(traits, tree300)
dim(pd$X_scaled)   # 300 x p_latent

# Multi-obs per species (via species_col)
pd2 <- preprocess_traits(avonet300, tree300, species_col = "Species_Key")
pd2$n_obs      # number of observations
pd2$n_species  # number of unique species




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("preprocess_traits", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pull_gbif_centroids")
### * pull_gbif_centroids

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pull_gbif_centroids
### Title: Fetch species range-centroid covariates from GBIF
### Aliases: pull_gbif_centroids

### ** Examples

## Don't show: 
if (interactive()) withAutoprint({ # examplesIf
## End(Don't show)
sp <- c("Quercus alba", "Pinus taeda", "Acer saccharum")
cov <- pull_gbif_centroids(sp, cache_dir = "script/data-cache/gbif")
cov_num <- cov[, c("centroid_lat", "centroid_lon"), drop = FALSE]
## Don't show: 
}) # examplesIf
## End(Don't show)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pull_gbif_centroids", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pull_worldclim_per_species")
### * pull_worldclim_per_species

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pull_worldclim_per_species
### Title: Fetch per-species bioclim covariates from WorldClim v2.1
### Aliases: pull_worldclim_per_species

### ** Examples

## Don't show: 
if (interactive()) withAutoprint({ # examplesIf
## End(Don't show)
sp <- c("Quercus alba", "Pinus taeda", "Acer saccharum")
gbif_df <- pull_gbif_centroids(sp,
  cache_dir = "script/data-cache/gbif")
wc_df   <- pull_worldclim_per_species(
  species             = sp,
  gbif_cache_dir      = "script/data-cache/gbif",
  worldclim_cache_dir = "script/data-cache/worldclim")
# Combine for impute()
cov <- cbind(gbif_df[, c("centroid_lat", "centroid_lon")],
             wc_df[, grep("^bio", colnames(wc_df))])
## Don't show: 
}) # examplesIf
## End(Don't show)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pull_worldclim_per_species", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_traits")
### * read_traits

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_traits
### Title: Read trait data from a CSV file or data frame
### Aliases: read_traits

### ** Examples

df <- data.frame(species = c("Sp_a", "Sp_b"), mass = c(10, 20))
traits <- read_traits(df, species_col = "species")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_traits", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_tree")
### * read_tree

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_tree
### Title: Read a phylogenetic tree from a file
### Aliases: read_tree

### ** Examples

## No test: 
path <- tempfile(fileext = ".tre")
ape::write.tree(ape::rtree(10L), path)
tree <- read_tree(path)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_tree", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("save_pigauto")
### * save_pigauto

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: save_pigauto
### Title: Save a fitted pigauto model
### Aliases: save_pigauto

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
traits$Mass[seq_len(3L)] <- NA_real_
fit <- impute(traits, tree, epochs = 5L, verbose = FALSE)$fit
path <- tempfile(fileext = ".pigauto")
save_pigauto(fit, path)
fit2 <- load_pigauto(path)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("save_pigauto", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_benchmark")
### * simulate_benchmark

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_benchmark
### Title: Run a simulation benchmark for pigauto
### Aliases: simulate_benchmark

### ** Examples

## No test: 
bench <- simulate_benchmark(n_species = 20L, n_traits = 2L,
                            scenarios = "BM", epochs = 5L, n_reps = 1L,
                            verbose = FALSE, gnn = FALSE)
bench$summary
plot(bench)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_benchmark", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_non_bm")
### * simulate_non_bm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_non_bm
### Title: Simulate non-BM trait data for benchmarking
### Aliases: simulate_non_bm

### ** Examples

tree <- ape::rtree(12)
traits <- simulate_non_bm(tree, n_traits = 2, scenario = "OU", seed = 1)
head(traits)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_non_bm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("suggest_next_observation")
### * suggest_next_observation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: suggest_next_observation
### Title: Suggest model-based candidate observations
### Aliases: suggest_next_observation

### ** Examples

## No test: 
data(avonet300, tree300, package = "pigauto")
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
traits$Mass[seq_len(3L)] <- NA_real_
res <- impute(traits, tree, epochs = 5L, verbose = FALSE)
suggest_next_observation(res, top_n = 5L)
suggest_next_observation(res, top_n = 5L, by = "species")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("suggest_next_observation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.pigauto_fit")
### * summary.pigauto_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.pigauto_fit
### Title: Summary method for pigauto_fit objects
### Aliases: summary.pigauto_fit

### ** Examples

## No test: 
data(avonet300, tree300)
tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
                     c("Mass", "Wing.Length"), drop = FALSE]
rownames(traits) <- tree$tip.label
pd <- preprocess_traits(traits, tree)
splits <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map)
fit <- fit_pigauto(pd, tree, splits = splits, epochs = 5L,
                   verbose = FALSE)
summary(fit, data = pd)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.pigauto_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("with_imputations")
### * with_imputations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: with_imputations
### Title: Fit a downstream model on every imputed dataset
### Aliases: with_imputations

### ** Examples

## No test: 
dat <- data.frame(y = stats::rnorm(30L), x = stats::rnorm(30L),
                  z = stats::rnorm(30L))
dat$x[seq(3L, 30L, by = 5L)] <- NA_real_
mi <- multi_impute_analysis(
  data = dat, formula = y ~ x + z, missing = "x",
  model = "lm", m = 2L
)
fits <- with_imputations(mi, function(d) stats::lm(y ~ x + z, data = d))
pool_mi(fits)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("with_imputations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
