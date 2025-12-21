# ==============================================================================
# AVONET: Mixed-type imputation engines (species-wise holdout)
# - Targets: 4 continuous traits + 1 categorical variable (cat_var)
# - Engine A: Rphylopars(BM) on continuous + dummy-coded categoricals -> softmax probs
# - Engine B: rMIDAS multiple imputations -> categorical probs by frequency
# ==============================================================================
suppressPackageStartupMessages({
  library(ape)
  library(Rphylopars)
  library(dplyr)
  library(Matrix)
  library(here)
})

# ------------------------- USER SETTINGS --------------------------------------
trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

cat_var <- "Diet"      # <--- CHANGE THIS if needed (one categorical variable)

holdout_frac <- 0.20   # species-wise holdout fraction
seed_split   <- 555

# Predictor selection caps (keep phylopars tractable)
max_num_pred <- 20      # numeric predictors (besides the 4 targets)
max_cat_pred <- 6       # other categorical predictors (low-cardinality)
max_levels_cat_pred <- 10

# rMIDAS settings
M_imputations <- 30     # number of completed datasets for categorical probabilities
use_rMIDAS <- TRUE      # set FALSE if you only want Engine A

# -----------------------------------------------------------------------------

# ------------------------- Helpers -------------------------------------------
softmax_rows <- function(mat) {
  # stable softmax by row
  mat <- as.matrix(mat)
  mat <- mat - apply(mat, 1, max, na.rm = TRUE)
  ex <- exp(mat)
  ex / rowSums(ex)
}

rmse <- function(truth, pred) sqrt(mean((truth - pred)^2))
acc  <- function(truth, pred) mean(truth == pred)

as_unknown_factor <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  factor(x)
}

onehot_from_factor <- function(f, prefix) {
  f <- as_unknown_factor(f)
  mm <- model.matrix(~ f - 1)
  colnames(mm) <- paste0(prefix, "=", sub("^f", "", colnames(mm)))
  mm
}

zscore <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
  if (is.finite(sdv) && sdv > 0) (x - mu)/sdv else (x - mu)
}

# -----------------------------------------------------------------------------

# ------------------------- Load AVONET + Tree --------------------------------
tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))

avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

stopifnot(all(trait_cols %in% names(avonet)))

if (!(cat_var %in% names(avonet))) {
  cat("\ncat_var not found. Here are some candidate categorical columns:\n")
  cand <- names(avonet)[sapply(avonet, function(x) !is.numeric(x))]
  cand <- setdiff(cand, c("Species3","Species_Key"))
  print(head(cand, 60))
  stop("\nSet cat_var to one of the columns above (or another categorical column).")
}

# Keep species with complete truth for continuous targets AND non-missing truth for cat_var
df_truth <- avonet %>%
  transmute(
    Species_Key,
    across(all_of(trait_cols), ~ .),
    .cat = .data[[cat_var]]
  ) %>%
  filter(if_all(all_of(trait_cols), ~ !is.na(.))) %>%
  filter(!is.na(.cat) & .cat != "")

common <- intersect(tree$tip.label, df_truth$Species_Key)
tree_pruned <- keep.tip(tree, common)
df_truth <- df_truth[match(tree_pruned$tip.label, df_truth$Species_Key), ]

n <- nrow(df_truth)
cat(sprintf("Aligned dataset: %d species | targets: %d continuous + 1 categorical (%s)\n",
            n, length(trait_cols), cat_var))

# ------------------------- Build predictor pool -------------------------------
# We use "all other variables" but cap them to keep things feasible.
# - numeric predictors: choose those with lowest missingness among aligned species
# - categorical predictors: choose low-cardinality columns with lowest missingness

df_all <- avonet[match(df_truth$Species_Key, avonet$Species_Key), , drop=FALSE]

# numeric predictor candidates
num_cols <- names(df_all)[sapply(df_all, is.numeric)]
num_cols <- setdiff(num_cols, c(trait_cols))   # don't duplicate targets

# rank numeric cols by missingness (lower is better)
num_miss <- sapply(df_all[, num_cols, drop=FALSE], function(x) mean(is.na(x)))
num_cols <- num_cols[order(num_miss)]
num_cols <- head(num_cols, max_num_pred)

# categorical predictor candidates
cat_cols <- names(df_all)[!sapply(df_all, is.numeric)]
cat_cols <- setdiff(cat_cols, c("Species3","Species_Key"))
cat_cols <- setdiff(cat_cols, cat_var)

# keep low-cardinality + not-too-missing
cat_ok <- c()
for (cc in cat_cols) {
  x <- df_all[[cc]]
  x <- as.character(x)
  x[x == ""] <- NA
  miss <- mean(is.na(x))
  levs <- length(unique(x[!is.na(x)]))
  if (miss < 0.5 && levs >= 2 && levs <= max_levels_cat_pred) cat_ok <- c(cat_ok, cc)
}
# rank by missingness
cat_miss <- sapply(df_all[, cat_ok, drop=FALSE], function(x) mean(is.na(x) | x==""))
cat_ok <- cat_ok[order(cat_miss)]
cat_ok <- head(cat_ok, max_cat_pred)

cat(sprintf("Using predictors: %d numeric + %d categorical (capped)\n",
            length(num_cols), length(cat_ok)))
if (length(cat_ok) > 0) print(cat_ok)

# ------------------------- Prepare target matrices ----------------------------
# Continuous targets: standardized log scale (as before)
X_raw <- as.matrix(df_truth[, trait_cols])
X_log <- log(X_raw)
X_truth <- scale(X_log)

# Categorical truth
y_truth <- as_unknown_factor(df_truth$.cat)
y_levels <- levels(y_truth)

# ------------------------- Species-wise holdout split -------------------------
set.seed(seed_split)
holdout_species <- sample.int(n, size = floor(holdout_frac * n), replace = FALSE)
is_holdout <- rep(FALSE, n); is_holdout[holdout_species] <- TRUE

cat(sprintf("Holdout species: %d (%.1f%%)\n", sum(is_holdout), 100*mean(is_holdout)))

# ------------------------- Build predictor matrices/dataframes ----------------
# numeric predictors (z-scored, no log transform assumed)
Z_num <- NULL
if (length(num_cols) > 0) {
  Z_num <- sapply(num_cols, function(cc) zscore(df_all[[cc]]))
  Z_num <- as.matrix(Z_num)
  colnames(Z_num) <- paste0("NUM__", num_cols)
}

# categorical predictors one-hot (with Unknown)
Z_cat <- NULL
if (length(cat_ok) > 0) {
  oh_list <- lapply(cat_ok, function(cc) onehot_from_factor(df_all[[cc]], prefix=paste0("CAT__", cc)))
  Z_cat <- do.call(cbind, oh_list)
}

# one-hot for target categorical (for Engine A)
Ycat_oh_truth <- onehot_from_factor(df_truth$.cat, prefix=paste0("TGT__", cat_var))
tgt_oh_cols <- colnames(Ycat_oh_truth)

# ------------------------- ENGINE A: Rphylopars BM ----------------------------
# Build a big multivariate trait matrix:
# [continuous targets] + [numeric predictors] + [categorical predictors one-hot] + [target cat one-hot]
# Then set targets missing for holdout species, keep predictors observed.
cat("\n--- Engine A: Rphylopars(BM) mixed (continuous + dummy categorical) ---\n")

T_all <- cbind(
  X_truth,
  if (!is.null(Z_num)) Z_num else NULL,
  if (!is.null(Z_cat)) Z_cat else NULL,
  Ycat_oh_truth
)

T_in <- T_all

# Mask targets for holdout species:
# - continuous targets columns 1:4
T_in[is_holdout, 1:ncol(X_truth)] <- NA
# - target categorical dummies at the end
T_in[is_holdout, tgt_oh_cols] <- NA

df_ph <- data.frame(species = df_truth$Species_Key, T_in)
fitA <- phylopars(trait_data = df_ph, tree = tree_pruned, model = "BM", pheno_error = FALSE)

T_hat <- fitA$anc_recon[df_truth$Species_Key, colnames(T_all), drop=FALSE]

# Continuous predictions (on standardized log scale)
X_hat_A <- as.matrix(T_hat[, colnames(X_truth), drop=FALSE])
rmse_A <- rmse(X_truth[is_holdout, ], X_hat_A[is_holdout, ])

# Categorical probabilities via softmax over dummy scores
S_hat <- as.matrix(T_hat[, tgt_oh_cols, drop=FALSE])
P_hat_A <- softmax_rows(S_hat)
colnames(P_hat_A) <- sub(paste0("^TGT__", cat_var, "="), "", tgt_oh_cols)

# Align probability columns to y_levels (ensure all levels present)
for (lv in setdiff(y_levels, colnames(P_hat_A))) {
  P_hat_A <- cbind(P_hat_A, setNames(rep(0, n), lv))
}
P_hat_A <- P_hat_A[, y_levels, drop=FALSE]
P_hat_A <- P_hat_A / rowSums(P_hat_A)  # renormalize

y_hat_A <- factor(y_levels[max.col(P_hat_A, ties.method="first")], levels=y_levels)
acc_A <- acc(y_truth[is_holdout], y_hat_A[is_holdout])

cat(sprintf("Engine A | continuous RMSE (holdout species): %.4f\n", rmse_A))
cat(sprintf("Engine A | categorical accuracy (holdout species): %.4f\n", acc_A))

# Show a few probabilistic imputations (holdout)
cat("\nEngine A: example categorical probabilities (first 6 holdout species):\n")
print(round(P_hat_A[which(is_holdout)[1:min(6, sum(is_holdout))], , drop=FALSE], 3))

# ------------------------- ENGINE B: rMIDAS (probabilistic via MI) ------------
# We estimate category probabilities by repeated completions (M imputations).
# This works even if rMIDAS doesn't directly expose soft probabilities.
if (use_rMIDAS) {
  if (!requireNamespace("rMIDAS", quietly = TRUE)) {
    stop("rMIDAS not installed. Run: remotes::install_github('MIDASverse/rMIDAS')")
  }
  library(rMIDAS)
  
  cat("\n--- Engine B: rMIDAS mixed-type (multiple imputations) ---\n")
  
  # Build mixed dataframe:
  df_mixed <- data.frame(
    Species_Key = df_truth$Species_Key,
    X_truth
  )
  colnames(df_mixed)[-1] <- trait_cols
  
  # add predictors
  if (!is.null(Z_num)) df_mixed <- cbind(df_mixed, as.data.frame(Z_num))
  if (length(cat_ok) > 0) {
    for (cc in cat_ok) df_mixed[[paste0("CAT__", cc)]] <- as_unknown_factor(df_all[[cc]])
  }
  
  # target categorical
  df_mixed[[paste0("TGT__", cat_var)]] <- y_truth
  
  # species-wise missingness on targets (continuous + target cat)
  for (cc in trait_cols) df_mixed[is_holdout, cc] <- NA
  df_mixed[is_holdout, paste0("TGT__", cat_var)] <- NA
  
  # Fit MIDAS model
  # NOTE: rMIDAS API can vary; this is written to be reasonably robust.
  midas_fit <- rMIDAS::midas(data = df_mixed, seed = 1)
  midas_fit <- rMIDAS::train(midas_fit)
  
  # Generate M completed datasets:
  completes <- vector("list", M_imputations)
  
  # Try native multi-imputation first; if not supported, fall back to repeated complete() calls
  got_multi <- FALSE
  try({
    tmp <- rMIDAS::complete(midas_fit, m = M_imputations)
    if (is.list(tmp) && length(tmp) == M_imputations) {
      completes <- tmp
      got_multi <- TRUE
    }
  }, silent = TRUE)
  
  if (!got_multi) {
    for (m in 1:M_imputations) {
      set.seed(1000 + m)
      completes[[m]] <- rMIDAS::complete(midas_fit)
    }
  }
  
  # Continuous: average across imputations
  X_hat_list <- lapply(completes, function(dd) as.matrix(dd[, trait_cols, drop=FALSE]))
  X_hat_B <- Reduce(`+`, X_hat_list) / M_imputations
  rmse_B <- rmse(X_truth[is_holdout, ], X_hat_B[is_holdout, ])
  
  # Categorical probabilities: frequency across imputations
  y_imp <- sapply(completes, function(dd) as.character(dd[[paste0("TGT__", cat_var)]]))
  # y_imp: n x M matrix (as characters)
  probs_B <- matrix(0, nrow=n, ncol=length(y_levels))
  colnames(probs_B) <- y_levels
  for (i in 1:n) {
    tab <- table(factor(y_imp[i, ], levels = y_levels))
    probs_B[i, ] <- as.numeric(tab) / sum(tab)
  }
  y_hat_B <- factor(y_levels[max.col(probs_B, ties.method="first")], levels=y_levels)
  acc_B <- acc(y_truth[is_holdout], y_hat_B[is_holdout])
  
  cat(sprintf("Engine B | continuous RMSE (holdout species): %.4f\n", rmse_B))
  cat(sprintf("Engine B | categorical accuracy (holdout species): %.4f\n", acc_B))
  
  cat("\nEngine B: example categorical probabilities (first 6 holdout species):\n")
  print(round(probs_B[which(is_holdout)[1:min(6, sum(is_holdout))], , drop=FALSE], 3))
}

# ------------------------- Summary table --------------------------------------
cat("\n================ SUMMARY (holdout species) ================\n")
cat(sprintf("Engine A (Rphylopars BM) | RMSE %.4f | Acc %.4f\n", rmse_A, acc_A))
if (use_rMIDAS) cat(sprintf("Engine B (rMIDAS)        | RMSE %.4f | Acc %.4f\n", rmse_B, acc_B))
cat("===========================================================\n")