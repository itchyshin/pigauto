args <- commandArgs(TRUE)
out <- args[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
base <- paste0("https://raw.githubusercontent.com/daniel1noble/BACE/",
               "b4b3cde02c4b171947042100e26b13b9f7ac425d/vignettes/amphibian_ctmax/data/")
for (f in c("ctmax_data.rds", "ctmax_tree.rds")) {
  path <- file.path(out, f)
  if (!file.exists(path)) download.file(paste0(base, f), path, mode = "wb")
}
d <- readRDS(file.path(out, "ctmax_data.rds"))
d <- d[d$row_type == "experimental", ]
rownames(d) <- NULL
d$medium[is.na(d$medium)] <- "body_water"
d$medium <- droplevels(d$medium)
tree <- ape::keep.tip(readRDS(file.path(out, "ctmax_tree.rds")), unique(d$species))
stopifnot(nrow(d) == 2661L, length(unique(d$species)) == 524L,
          !anyDuplicated(d$row_n), all(is.finite(d$CTmax)))
covs <- d[, c("acclimation_temp", "medium", "endpoint", "acclimated", "life_stage", "ecotype")]
stopifnot(!anyNA(covs))
# Stratified row folds: repeated observations are distributed across folds.
# Species with fewer than five observations can occur in fewer folds.
set.seed(2601109L)
fold <- integer(nrow(d))
for (idx in split(seq_len(nrow(d)), d$species)) {
  assignment <- rep(sample.int(5L), length.out = length(idx))
  fold[idx] <- assignment[sample.int(length(idx))]
}
stopifnot(setequal(unique(fold), 1:5))
saveRDS(list(data = d, tree = tree, covariates = covs, fold = fold,
             split_seed = 2601109L, fit_seeds = 51001:51006,
             input_md5 = tools::md5sum(file.path(out, c("ctmax_data.rds", "ctmax_tree.rds")))),
        file.path(out, "prepared.rds"))
write.csv(data.frame(row_n=d$row_n, species=d$species, fold=fold),
          file.path(out, "folds.csv"), row.names = FALSE)
cat("PREPARED", nrow(d), "rows", length(tree$tip.label), "species\n")
