## data-raw/make_avonet_full.R
## Creates avonet_full (data.frame) and tree_full (phylo) package datasets.
## Run from the package root: Rscript data-raw/make_avonet_full.R
##
## This is the full-scale counterpart to avonet300 / tree300: all 9,993 bird
## species for which AVONET3 and the BirdTree Stage2 Hackett MCC phylogeny
## agree. The schema is identical to avonet300 (same trait columns, same
## Species_Key column, same factor encodings) so any code that runs on
## avonet300 runs on avonet_full with no modification.
##
## Unlike the 300-species subset, native AVONET missingness is PRESERVED --
## users see the real-world missingness pattern as it comes from AVONET3.

library(ape)
library(here)

# ---- Load raw data -----------------------------------------------------------
tree_path <- here("avonet", "Stage2_Hackett_MCC_no_neg.tre")
csv_path  <- here("avonet", "AVONET3_BirdTree.csv")

if (!file.exists(tree_path)) stop("Tree file not found: ", tree_path)
if (!file.exists(csv_path))  stop("CSV file not found: ", csv_path)

tree   <- ape::read.tree(tree_path)
avonet <- read.csv(csv_path, stringsAsFactors = FALSE)
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

# ---- Select traits -----------------------------------------------------------
# Same schema as avonet300: 4 continuous morphometric + 2 categorical + 1 ordinal.
cont_cols  <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
trait_cols <- c(cont_cols, "Trophic.Level", "Primary.Lifestyle", "Migration")

# ---- Clean and prepare -------------------------------------------------------
# Require all continuous traits complete; allow NAs in ecological traits.
# This matches the avonet300 build so that any NA cells in avonet_full come
# from the ecological (categorical/ordinal) columns only.
df <- avonet[stats::complete.cases(avonet[, cont_cols]), ]
df <- df[, c("Species_Key", trait_cols)]

df$Trophic.Level <- trimws(df$Trophic.Level)
df$Trophic.Level[df$Trophic.Level == ""] <- NA
df$Trophic.Level <- factor(df$Trophic.Level)

df$Primary.Lifestyle <- trimws(df$Primary.Lifestyle)
df$Primary.Lifestyle[df$Primary.Lifestyle == ""] <- NA
df$Primary.Lifestyle <- factor(df$Primary.Lifestyle)

df$Migration <- ordered(df$Migration, levels = c(1, 2, 3),
                        labels = c("Resident", "Partial", "Full"))

# ---- Align to tree -----------------------------------------------------------
common <- intersect(tree$tip.label, df$Species_Key)
df     <- df[df$Species_Key %in% common, ]

cat("Species with complete continuous data in both tree and CSV:",
    length(common), "\n")

# ---- Build objects (no subsampling) -----------------------------------------
tree_full   <- ape::keep.tip(tree, common)
avonet_full <- df[match(tree_full$tip.label, df$Species_Key), ]
rownames(avonet_full) <- NULL

cat("avonet_full:", nrow(avonet_full), "rows x", ncol(avonet_full), "cols\n")
cat("tree_full:  ", length(tree_full$tip.label), "tips\n")
cat("\nTrait types:\n")
for (col in names(avonet_full)) {
  cat(sprintf("  %-20s %s\n", col, paste(class(avonet_full[[col]]), collapse = "/")))
}
cat("\nMissing values:\n")
for (col in names(avonet_full)) {
  n_na <- sum(is.na(avonet_full[[col]]))
  if (n_na > 0) {
    cat(sprintf("  %-20s %d NAs (%.2f%%)\n",
                col, n_na, 100 * n_na / nrow(avonet_full)))
  }
}

stopifnot(all(tree_full$tip.label == avonet_full$Species_Key))

# ---- Save -------------------------------------------------------------------
# xz compression (level 9) gives the smallest on-disk footprint for character
# + factor vectors and for phylo tip-label vectors. Combined size is ~360 KB,
# well under the CRAN 5 MB tarball limit.
usethis::use_data(avonet_full, tree_full,
                  overwrite = TRUE, compress = "xz")

cat("Saved data/avonet_full.rda and data/tree_full.rda\n")
cat("File sizes:\n")
for (f in c("data/avonet_full.rda", "data/tree_full.rda")) {
  if (file.exists(f)) {
    cat(sprintf("  %-25s %7.1f KB\n", f, file.size(f) / 1024))
  }
}
