# Generate a synthetic fixed-lambda regression oracle using historical pigauto.
# Run from the repository root with:
#   Rscript data-raw/generate-lambda-reference.R
# The required historical commit must be available in this repository.
# This replaces a historical fixture whose input provenance was not recorded.
# The seed, generated inputs, and historical solver output are retained in the
# fixture. Exact floating-point output may vary with R/BLAS/package versions.
stopifnot(length(commandArgs(trailingOnly = TRUE)) == 0L)
oracle_commit <- "ab02e316ec72eabb57dad5607d9f43fb0a6fbcc1"
archive <- tempfile(fileext = ".tar")
oracle_dir <- tempfile("pigauto-oracle-")
dir.create(oracle_dir)
archive_status <- system2("git", c("archive", "--format=tar", oracle_commit),
                          stdout = archive)
if (!identical(archive_status, 0L)) stop("Failed to archive the pinned oracle commit")
extract_status <- system2("tar", c("-xf", shQuote(archive), "-C", shQuote(oracle_dir)))
if (!identical(extract_status, 0L)) stop("Failed to extract the pinned oracle commit")
pkgload::load_all(oracle_dir, quiet = TRUE)
set.seed(20261108L)
n <- 60L
tree <- ape::rtree(n)
continuous <- matrix(rnorm(n * 3L), n, 3L)
colnames(continuous) <- c("c1", "c2", "c3")
df <- data.frame(continuous,
                 b1 = factor(sample(c("no", "yes"), n, replace = TRUE)),
                 k1 = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
                 row.names = tree$tip.label)
for (j in seq_along(df)) df[sample.int(n, 12L), j] <- NA
L <- scale(as.matrix(df[, 1:3]))
rownames(L) <- tree$tip.label
solver <- pigauto:::fit_mvn_bm_inhouse(L, tree = tree)
pd <- pigauto::preprocess_traits(df, tree)
splits <- pigauto::make_missing_splits(pd$X_scaled, trait_map = pd$trait_map,
                                      seed = 20261109L)
baseline <- pigauto::fit_baseline(pd, tree, splits = splits)
ref <- list(commit = "ab02e31", seed = 20261108L, split_seed = 20261109L,
            generator = "data-raw/generate-lambda-reference.R",
            R_version = R.version.string,
            ape_version = as.character(utils::packageVersion("ape")),
            tree = tree, L = L,
            solver = list(anc_recon = solver$anc_recon,
                          anc_var = solver$anc_var,
                          phylocov = solver$pars$phylocov),
            df = df, splits = splits,
            baseline = list(mu = baseline$mu, se = baseline$se,
                            path = baseline$path))
saveRDS(ref, "tests/testthat/fixtures/lambda_fixed1_reference_ab02e31.rds",
        version = 3)
unlink(c(archive, oracle_dir), recursive = TRUE)
