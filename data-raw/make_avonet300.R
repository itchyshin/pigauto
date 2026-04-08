## data-raw/make_avonet300.R
## Creates avonet300 (data.frame) and tree300 (phylo) package datasets.
## Run from the package root: Rscript data-raw/make_avonet300.R
##
## The dataset includes 4 continuous morphometric traits plus 4 ecological
## traits (categorical, ordinal, binary) to demonstrate mixed-type imputation.

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
# Continuous morphometric traits
cont_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

# Ecological traits (categorical, ordinal)
# Trophic.Level: 4 categories (Carnivore, Herbivore, Omnivore, Scavenger)
# Primary.Lifestyle: 5 categories (Insessorial, Terrestrial, Generalist, Aerial, Aquatic)
# Habitat: many categories — collapse rare ones for robustness
# Migration: 3 ordinal levels (1=resident, 2=partial migrant, 3=full migrant)

trait_cols <- c(cont_cols, "Trophic.Level", "Primary.Lifestyle", "Migration")

# ---- Clean and prepare -------------------------------------------------------
# Require all continuous traits complete; allow NAs in ecological traits
df <- avonet[stats::complete.cases(avonet[, cont_cols]), ]
df <- df[, c("Species_Key", trait_cols)]

# Clean Trophic.Level: trim whitespace, convert to factor
df$Trophic.Level <- trimws(df$Trophic.Level)
df$Trophic.Level[df$Trophic.Level == ""] <- NA
df$Trophic.Level <- factor(df$Trophic.Level)

# Clean Primary.Lifestyle: trim whitespace, convert to factor
df$Primary.Lifestyle <- trimws(df$Primary.Lifestyle)
df$Primary.Lifestyle[df$Primary.Lifestyle == ""] <- NA
df$Primary.Lifestyle <- factor(df$Primary.Lifestyle)

# Migration: ordinal (1=resident, 2=partial, 3=full migrant)
df$Migration <- ordered(df$Migration, levels = c(1, 2, 3),
                        labels = c("Resident", "Partial", "Full"))

# ---- Align to tree -----------------------------------------------------------
common <- intersect(tree$tip.label, df$Species_Key)
df     <- df[df$Species_Key %in% common, ]

cat("Species with complete continuous data in both tree and CSV:", length(common), "\n")

# ---- Sample 300 species ------------------------------------------------------
set.seed(42)
keep_300 <- sample(common, 300)

# ---- Build objects -----------------------------------------------------------
avonet300 <- df[df$Species_Key %in% keep_300, ]
rownames(avonet300) <- NULL

# Re-order to match the pruned tree tip order
tree300   <- ape::keep.tip(tree, keep_300)
avonet300 <- avonet300[match(tree300$tip.label, avonet300$Species_Key), ]
rownames(avonet300) <- NULL

cat("avonet300:", nrow(avonet300), "rows x", ncol(avonet300), "cols\n")
cat("tree300:  ", length(tree300$tip.label), "tips\n")
cat("\nTrait types:\n")
for (col in names(avonet300)) {
  cat(sprintf("  %-20s %s\n", col, paste(class(avonet300[[col]]), collapse = "/")))
}
cat("\nMissing values:\n")
for (col in names(avonet300)) {
  n_na <- sum(is.na(avonet300[[col]]))
  if (n_na > 0) cat(sprintf("  %-20s %d NAs (%.1f%%)\n", col, n_na, 100 * n_na / 300))
}

stopifnot(all(tree300$tip.label == avonet300$Species_Key))
stopifnot(nrow(avonet300) == 300)

# ---- Save -------------------------------------------------------------------
usethis::use_data(avonet300, tree300, overwrite = TRUE)
cat("Saved data/avonet300.rda and data/tree300.rda\n")
