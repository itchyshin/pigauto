## data-raw/make_avonet300.R
## Creates avonet300 (data.frame) and tree300 (phylo) package datasets.
## Run from the package root: Rscript data-raw/make_avonet300.R

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

# ---- Select traits and align -------------------------------------------------
trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
df <- stats::na.omit(avonet[, c("Species_Key", trait_cols)])

# Keep only species present in both tree and data, with all traits complete
common <- intersect(tree$tip.label, df$Species_Key)
df     <- df[df$Species_Key %in% common, ]

cat("Species with complete data in both tree and CSV:", length(common), "\n")

# ---- Sample 300 species ------------------------------------------------------
set.seed(42)
keep_300 <- sample(common, 300)

# ---- Build objects -----------------------------------------------------------
avonet300 <- df[df$Species_Key %in% keep_300, ]
rownames(avonet300) <- NULL  # clean row names for a nice data.frame

# Re-order to match the pruned tree tip order
tree300   <- ape::keep.tip(tree, keep_300)
avonet300 <- avonet300[match(tree300$tip.label, avonet300$Species_Key), ]
rownames(avonet300) <- NULL

cat("avonet300:", nrow(avonet300), "rows x", ncol(avonet300), "cols\n")
cat("tree300:  ", length(tree300$tip.label), "tips\n")

stopifnot(all(tree300$tip.label == avonet300$Species_Key))
stopifnot(nrow(avonet300) == 300)

# ---- Save -------------------------------------------------------------------
usethis::use_data(avonet300, tree300, overwrite = TRUE)
cat("Saved data/avonet300.rda and data/tree300.rda\n")
