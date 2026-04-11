## data-raw/make_delhey5809.R
## Creates delhey5809 (data.frame) and tree_delhey (phylo) package datasets.
## Run from the package root: Rscript data-raw/make_delhey5809.R
##
## Source: Delhey et al. (2019) "A comparative analysis of the colour of
## bird plumage in relation to environmental variables."
## 5,809 passerine bird species with plumage lightness + climate covariates.

library(ape)
library(megatrees)
library(here)

# ---- Load raw data -----------------------------------------------------------
csv_path <- here("useful", "Delhey2019_plumage_lightness.csv")
if (!file.exists(csv_path)) stop("CSV not found: ", csv_path)

d <- read.csv(csv_path, stringsAsFactors = FALSE)
cat("Raw data:", nrow(d), "rows x", ncol(d), "cols\n")

# Clean up: drop the row-number column 'X'
d$X <- NULL

# Rename for clarity
names(d)[names(d) == "TipLabel"] <- "Species_Key"

# ---- Build MCC tree ----------------------------------------------------------
# Use the MCC tree from the bundled AVONET Stage2 Hackett tree
tree_path <- here("avonet", "Stage2_Hackett_MCC_no_neg.tre")
if (!file.exists(tree_path)) stop("Tree not found: ", tree_path)

tree_full <- ape::read.tree(tree_path)
cat("Full tree tips:", length(tree_full$tip.label), "\n")

# Check overlap
n_match <- sum(d$Species_Key %in% tree_full$tip.label)
cat("Species matching MCC tree:", n_match, "/", nrow(d), "\n")

# Prune tree to Delhey species
tree_delhey <- ape::keep.tip(tree_full, d$Species_Key)
cat("Pruned tree:", length(tree_delhey$tip.label), "tips\n")

# Reorder data to match tree tip order
d <- d[match(tree_delhey$tip.label, d$Species_Key), ]
rownames(d) <- NULL

stopifnot(all(d$Species_Key == tree_delhey$tip.label))

# ---- Build the dataset -------------------------------------------------------
# Trait columns: lightness_male, lightness_female (continuous, to be imputed)
# Covariate columns: environmental variables (fully observed conditioners)
# Also keep family for potential random-effect use

delhey5809 <- d

cat("\nFinal dataset:\n")
cat("  Rows:", nrow(delhey5809), "\n")
cat("  Cols:", ncol(delhey5809), "\n")
for (col in names(delhey5809)) {
  cl <- class(delhey5809[[col]])
  n_na <- sum(is.na(delhey5809[[col]]))
  cat(sprintf("  %-45s %s   NAs: %d\n", col, cl, n_na))
}

# ---- Save --------------------------------------------------------------------
save(delhey5809, file = here("data", "delhey5809.rda"), compress = "xz")
save(tree_delhey, file = here("data", "tree_delhey.rda"), compress = "xz")

f1 <- file.size(here("data", "delhey5809.rda"))
f2 <- file.size(here("data", "tree_delhey.rda"))
cat(sprintf("\nSaved delhey5809.rda (%.0f KB) and tree_delhey.rda (%.0f KB)\n",
            f1 / 1024, f2 / 1024))
