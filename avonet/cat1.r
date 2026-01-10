# ==============================================================================
# AVONET: Mixed-type imputation engines (species-wise holdout)
# - Targets: 4 continuous traits + 1 categorical variable (cat_var)
# - Engine A: Rphylopars(BM) on continuous + dummy-coded categoricals -> softmax probs
# - Engine B: Phylo-DAE (BM diffusion) mixed head -> continuous + categorical probs
# ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(Rphylopars)
  library(dplyr)
  library(Matrix)
  library(here)
  library(torch)
})

# ------------------------- USER SETTINGS --------------------------------------
trait_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
cat_var <- "Diet"      # if not present, we auto-pick a sensible one

holdout_frac <- 0.20
seed_split   <- 555

max_num_pred <- 20
max_cat_pred <- 6
max_levels_cat_pred <- 10

# Phylo-DAE settings
dae_hidden     <- 192
dae_Kdiff      <- 50
dae_alpha      <- 0.08
dae_dropout    <- 0.05
dae_lr         <- 1e-3
dae_wd         <- 1e-6
dae_max_epochs <- 6000
dae_eval_every <- 200
dae_patience   <- 15

# Denoising corruption
corrupt_x_p <- 0.40
corrupt_y_p <- 0.30
gamma_obs   <- 0.05

min_count_level <- 20

# ------------------------- Helpers -------------------------------------------
softmax_rows <- function(mat) {
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
to_num_mat <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) x <- as.matrix(x)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

# ---- Key: fit phylopars safely by relabeling tree to numeric IDs -------------
fit_phylopars_mixed_safe <- function(T_in, species_names, tree_pruned, model = "BM") {
  stopifnot(length(species_names) == nrow(T_in))
  stopifnot(all(tree_pruned$tip.label == species_names))  # required for clean mapping
  
  n <- length(species_names)
  ids <- seq_len(n)                           # 1..n
  tip_ids_chr <- as.character(ids)
  
  tree_id <- tree_pruned
  tree_id$tip.label <- tip_ids_chr
  
  df_ph <- as.data.frame(T_in, check.names = FALSE)
  df_ph[] <- lapply(df_ph, function(v) { v <- as.numeric(v); storage.mode(v) <- "double"; v })
  
  # phylopars requires first column named 'species'
  df_ph <- cbind(species = ids, df_ph)
  df_ph <- as.data.frame(df_ph, check.names = FALSE)
  
  # Fit
  fit <- phylopars(trait_data = df_ph, tree = tree_id, model = model, pheno_error = FALSE)
  
  # Extract and map back
  T_hat <- fit$anc_recon
  
  # Ensure rows are in 1..n order
  if (!is.null(rownames(T_hat))) {
    T_hat <- T_hat[tip_ids_chr, , drop = FALSE]
  } else {
    # fallback: assume already in input order
    T_hat <- T_hat
  }
  
  rownames(T_hat) <- species_names
  T_hat
}

# ------------------------- Load AVONET + Tree --------------------------------
tree <- read.tree(here("avonet","Stage2_Hackett_MCC_no_neg.tre"))
avonet <- read.csv(here("avonet","AVONET3_BirdTree.csv"))
avonet$Species_Key <- gsub(" ", "_", avonet$Species3)

stopifnot(all(trait_cols %in% names(avonet)))

if (!(cat_var %in% names(avonet))) {
  cand <- c("Trophic.Niche", "Trophic.Level", "Habitat", "Primary.Lifestyle", "Species.Status")
  pick <- cand[cand %in% names(avonet)][1]
  if (is.na(pick)) {
    cand2 <- names(avonet)[!sapply(avonet, is.numeric)]
    cand2 <- setdiff(cand2, c("Species3","Species_Key"))
    pick <- cand2[1]
  }
  message("cat_var '", cat_var, "' not found. Using cat_var = '", pick, "' instead.")
  cat_var <- pick
}

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
stopifnot(all(df_truth$Species_Key == tree_pruned$tip.label))

n <- nrow(df_truth)
cat(sprintf("Aligned dataset: %d species | targets: %d continuous + 1 categorical (%s)\n",
            n, length(trait_cols), cat_var))

# ------------------------- Build predictor pool -------------------------------
df_all <- avonet[match(df_truth$Species_Key, avonet$Species_Key), , drop=FALSE]

num_cols <- names(df_all)[sapply(df_all, is.numeric)]
num_cols <- setdiff(num_cols, trait_cols)
if (length(num_cols) > 0) {
  num_miss <- sapply(df_all[, num_cols, drop=FALSE], function(x) mean(is.na(x)))
  num_cols <- num_cols[order(num_miss)]
  num_cols <- head(num_cols, max_num_pred)
}

cat_cols <- names(df_all)[!sapply(df_all, is.numeric)]
cat_cols <- setdiff(cat_cols, c("Species3","Species_Key", cat_var))

cat_ok <- character(0)
for (cc in cat_cols) {
  x <- as.character(df_all[[cc]])
  x[x == ""] <- NA
  miss <- mean(is.na(x))
  levs <- length(unique(x[!is.na(x)]))
  if (miss < 0.5 && levs >= 2 && levs <= max_levels_cat_pred) cat_ok <- c(cat_ok, cc)
}
if (length(cat_ok) > 0) {
  cat_miss <- sapply(df_all[, cat_ok, drop=FALSE], function(x) mean(is.na(x) | x==""))
  cat_ok <- cat_ok[order(cat_miss)]
  cat_ok <- head(cat_ok, max_cat_pred)
}

cat(sprintf("Using predictors: %d numeric + %d categorical (capped)\n",
            length(num_cols), length(cat_ok)))
if (length(cat_ok) > 0) print(cat_ok)

# ------------------------- Prepare targets ------------------------------------
X_raw  <- as.matrix(df_truth[, trait_cols])
X_log  <- log(X_raw)
X_truth <- to_num_mat(scale(X_log))
p <- ncol(X_truth)

y_truth_raw <- as_unknown_factor(df_truth$.cat)
y_chr <- as.character(y_truth_raw)
tab <- table(y_chr)
rare <- names(tab)[tab < min_count_level & names(tab) != "Unknown"]
if (length(rare) > 0) y_chr[y_chr %in% rare] <- "Other"
y_truth <- factor(y_chr)
y_levels <- levels(y_truth)
K <- length(y_levels)

cat(sprintf("Target categorical '%s': %d classes (after rare->Other, min_count=%d)\n",
            cat_var, K, min_count_level))
print(sort(table(y_truth), decreasing = TRUE))

# ------------------------- Holdout split --------------------------------------
set.seed(seed_split)
holdout_species <- sample.int(n, size = floor(holdout_frac * n), replace = FALSE)
is_holdout <- rep(FALSE, n); is_holdout[holdout_species] <- TRUE
cat(sprintf("Holdout species: %d (%.1f%%)\n", sum(is_holdout), 100*mean(is_holdout)))

# ------------------------- Predictors -----------------------------------------
Z_num <- NULL
if (length(num_cols) > 0) {
  Z_num <- sapply(num_cols, function(cc) zscore(df_all[[cc]]))
  Z_num <- to_num_mat(Z_num)
  colnames(Z_num) <- paste0("NUM__", num_cols)
}

Z_cat <- NULL
if (length(cat_ok) > 0) {
  oh_list <- lapply(cat_ok, function(cc) onehot_from_factor(df_all[[cc]], prefix=paste0("CAT__", cc)))
  Z_cat <- to_num_mat(do.call(cbind, oh_list))
}

Ycat_oh_truth <- to_num_mat(onehot_from_factor(y_truth, prefix=paste0("TGT__", cat_var)))
tgt_oh_cols <- colnames(Ycat_oh_truth)

# ==============================================================================
# ENGINE A: Rphylopars (BM) mixed (SAFE RELABEL)
# ==============================================================================
cat("\n--- Engine A: Rphylopars(BM) mixed (continuous + dummy categorical) ---\n")

T_all <- cbind(
  X_truth,
  if (!is.null(Z_num)) Z_num else NULL,
  if (!is.null(Z_cat)) Z_cat else NULL,
  Ycat_oh_truth
)
T_all <- to_num_mat(T_all)

T_in <- T_all
T_in[is_holdout, 1:p] <- NA
T_in[is_holdout, tgt_oh_cols] <- NA

T_hat <- fit_phylopars_mixed_safe(
  T_in = T_in,
  species_names = df_truth$Species_Key,
  tree_pruned = tree_pruned,
  model = "BM"
)

# Continuous
X_hat_A <- as.matrix(T_hat[, colnames(X_truth), drop=FALSE])
rmse_A <- rmse(X_truth[is_holdout, ], X_hat_A[is_holdout, ])

# Categorical probs via softmax over dummy scores
S_hat <- as.matrix(T_hat[, tgt_oh_cols, drop=FALSE])
P_hat_A <- softmax_rows(S_hat)
colnames(P_hat_A) <- sub(paste0("^TGT__", cat_var, "="), "", tgt_oh_cols)

for (lv in setdiff(y_levels, colnames(P_hat_A))) P_hat_A <- cbind(P_hat_A, setNames(rep(0, n), lv))
P_hat_A <- P_hat_A[, y_levels, drop=FALSE]
P_hat_A <- P_hat_A / rowSums(P_hat_A)

y_hat_A <- factor(y_levels[max.col(P_hat_A, ties.method="first")], levels=y_levels)
acc_A <- acc(y_truth[is_holdout], y_hat_A[is_holdout])

cat(sprintf("Engine A | continuous RMSE (holdout species): %.4f\n", rmse_A))
cat(sprintf("Engine A | categorical accuracy (holdout species): %.4f\n", acc_A))
cat("\nEngine A: example categorical probabilities (first 6 holdout species):\n")
print(round(P_hat_A[which(is_holdout)[1:min(6, sum(is_holdout))], , drop=FALSE], 3))

# ==============================================================================
# ENGINE B: Phylo-DAE (BM diffusion) mixed head
# (Same as your block; left intact.)
# ==============================================================================

cat("\n--- Engine B: Phylo-DAE (BM diffusion) mixed-type ---\n")

device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

C_bm <- vcv.phylo(tree_pruned, corr = TRUE)
C_bm <- C_bm + diag(1e-6, nrow(C_bm))
rs <- rowSums(C_bm) + 1e-8
Dinv <- diag(1 / sqrt(rs))
A_op <- Dinv %*% C_bm %*% Dinv
t_A <- torch_tensor(A_op, dtype = torch_float(), device = device)

X_in2 <- X_truth
X_in2[is_holdout, ] <- NA
Mx <- !is.na(X_in2)
X_fill <- X_in2; X_fill[is.na(X_fill)] <- 0

Z_num_in <- to_num_mat(Z_num); if (!is.null(Z_num_in)) Z_num_in[is.na(Z_num_in)] <- 0
Z_cat_in <- to_num_mat(Z_cat); if (!is.null(Z_cat_in)) Z_cat_in[is.na(Z_cat_in)] <- 0

y_id <- as.integer(y_truth) - 1L
My <- rep(1, n); My[is_holdout] <- 0
unk_id <- match("Unknown", y_levels) - 1L
if (is.na(unk_id)) unk_id <- 0L

y_fill <- y_id
y_fill[is_holdout] <- unk_id

Yoh <- model.matrix(~ factor(y_fill, levels = 0:(K-1)) - 1)
storage.mode(Yoh) <- "double"

Mx_mat <- to_num_mat(Mx * 1)
My_vec <- matrix(as.double(My), ncol = 1)

Z_list <- list(
  to_num_mat(X_fill),
  Mx_mat,
  Z_num_in,
  Z_cat_in,
  to_num_mat(Yoh),
  to_num_mat(My_vec)
)
Z_list <- Z_list[!sapply(Z_list, is.null)]
Z_mat <- to_num_mat(do.call(cbind, Z_list))

t_Z  <- torch_tensor(Z_mat, dtype=torch_float(), device=device)
t_Xt <- torch_tensor(X_truth, dtype=torch_float(), device=device)
t_y  <- torch_tensor(y_id, dtype=torch_long(), device=device)
t_My <- torch_tensor(matrix(as.double(!is_holdout), ncol=1), dtype=torch_float(), device=device)

d_in <- ncol(Z_mat)
x_cols  <- 1:p
mx_cols <- (p+1):(2*p)
d_num <- if (is.null(Z_num_in)) 0 else ncol(Z_num_in)
d_cat <- if (is.null(Z_cat_in)) 0 else ncol(Z_cat_in)
y_start <- 2*p + d_num + d_cat + 1
y_end   <- y_start + K - 1
my_col  <- y_end + 1

PhyloDAE_Mixed <- nn_module(
  initialize = function(d_in, p_cont, K_cat, hidden, Kdiff, alpha, dropout) {
    self$Kdiff <- Kdiff
    self$alpha <- alpha
    
    self$enc <- nn_linear(d_in, hidden)
    self$ln1 <- nn_layer_norm(hidden)
    self$act <- nn_gelu()
    self$drop <- nn_dropout(dropout)
    
    self$hx1 <- nn_linear(hidden, hidden)
    self$lnx <- nn_layer_norm(hidden)
    self$hx2 <- nn_linear(hidden, p_cont)
    
    self$hy1 <- nn_linear(hidden, hidden)
    self$lny <- nn_layer_norm(hidden)
    self$hy2 <- nn_linear(hidden, K_cat)
    
    self$mask_token_x <- nn_parameter(torch_zeros(p_cont))
  },
  forward = function(Z, A) {
    h0 <- self$drop(self$act(self$ln1(self$enc(Z))))
    h <- h0
    for (k in 1:self$Kdiff) {
      h <- self$alpha * h0 + (1 - self$alpha) * torch_matmul(A, h)
    }
    out_x <- self$hx2(self$drop(self$act(self$lnx(self$hx1(h)))))
    logits_y <- self$hy2(self$drop(self$act(self$lny(self$hy1(h)))))
    list(out_x = out_x, logits_y = logits_y)
  }
)

model <- PhyloDAE_Mixed(d_in=d_in, p_cont=p, K_cat=K,
                        hidden=dae_hidden, Kdiff=dae_Kdiff,
                        alpha=dae_alpha, dropout=dae_dropout)$to(device=device)
opt <- optim_adamw(model$parameters, lr=dae_lr, weight_decay=dae_wd)

best_state <- NULL
best_score <- -Inf
pat_left <- dae_patience

for (epoch in 1:dae_max_epochs) {
  model$train()
  opt$zero_grad()
  Zc <- t_Z$clone()
  
  Mx_now <- Zc[, mx_cols]
  u <- torch_rand(c(n, p), device=device)
  corrupt_x <- (u < corrupt_x_p) * (Mx_now > 0.5)
  
  if (as.numeric(corrupt_x$sum()$cpu()) > 0) {
    token_mat <- model$mask_token_x$unsqueeze(1)$transpose(1,2)$expand(c(n, p))
    Zc[, x_cols]  <- torch_where(corrupt_x > 0, token_mat, Zc[, x_cols])
    Zc[, mx_cols] <- torch_where(corrupt_x > 0, torch_zeros_like(Zc[, mx_cols]), Zc[, mx_cols])
  }
  
  My_obs <- Zc[, my_col]
  u2 <- torch_rand(c(n), device=device)
  corrupt_y <- (u2 < corrupt_y_p) * (My_obs > 0.5)
  
  if (as.numeric(corrupt_y$sum()$cpu()) > 0) {
    unk_oh <- torch_zeros(c(n, K), device=device)
    unk_oh[, (unk_id + 1)] <- 1
    cy2 <- corrupt_y$unsqueeze(2)$expand(c(n, K))
    Zc[, y_start:y_end] <- torch_where(cy2 > 0, unk_oh, Zc[, y_start:y_end])
    Zc[, my_col] <- torch_where(corrupt_y > 0, torch_zeros_like(Zc[, my_col]), Zc[, my_col])
  }
  
  out <- model(Zc, t_A)
  pred_x <- out$out_x
  logits <- out$logits_y
  
  loss <- torch_tensor(0, device=device)
  
  if (as.numeric(corrupt_x$sum()$cpu()) > 0) {
    diff <- pred_x - t_Xt
    mse_corrupt <- (diff^2) * corrupt_x
    loss <- loss + mse_corrupt$sum() / (corrupt_x$sum() + 1e-8)
  }
  
  obs_mask <- (Mx_now > 0.5)
  diff2 <- pred_x - t_Xt
  mse_obs <- (diff2^2) * obs_mask
  loss <- loss + gamma_obs * (mse_obs$sum() / (obs_mask$sum() + 1e-8))
  
  if (as.numeric(corrupt_y$sum()$cpu()) > 0) {
    idx <- which(as.array(corrupt_y$cpu()) == 1)
    loss <- loss + nnf_cross_entropy(logits[idx, ], t_y[idx])
  }
  
  loss$backward()
  nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
  opt$step()
  
  if (epoch %% dae_eval_every == 0) {
    model$eval()
    with_no_grad({
      out0 <- model(t_Z, t_A)
      X_hat_B <- as.matrix(out0$out_x$cpu())
      L_hat   <- as.matrix(out0$logits_y$cpu())
      P_hat_B <- softmax_rows(L_hat)
      colnames(P_hat_B) <- y_levels
      y_hat_B <- factor(y_levels[max.col(P_hat_B, ties.method="first")], levels=y_levels)
    })
    
    rmse_B <- rmse(X_truth[is_holdout, ], X_hat_B[is_holdout, ])
    acc_B  <- acc(y_truth[is_holdout], y_hat_B[is_holdout])
    score <- acc_B - 0.02 * rmse_B
    
    cat(sprintf("epoch %4d | loss %.4f | holdout RMSE %.4f | holdout Acc %.4f\n",
                epoch, loss$item(), rmse_B, acc_B))
    
    if (score > best_score + 1e-6) {
      best_score <- score
      best_state <- model$state_dict()
      pat_left <- dae_patience
    } else {
      pat_left <- pat_left - 1
    }
    if (pat_left <= 0) break
  }
}

if (!is.null(best_state)) model$load_state_dict(best_state)

model$eval()
with_no_grad({
  out0 <- model(t_Z, t_A)
  X_hat_B <- as.matrix(out0$out_x$cpu())
  L_hat   <- as.matrix(out0$logits_y$cpu())
  P_hat_B <- softmax_rows(L_hat)
  colnames(P_hat_B) <- y_levels
  y_hat_B <- factor(y_levels[max.col(P_hat_B, ties.method="first")], levels=y_levels)
})

rmse_B <- rmse(X_truth[is_holdout, ], X_hat_B[is_holdout, ])
acc_B  <- acc(y_truth[is_holdout], y_hat_B[is_holdout])

cat(sprintf("\nEngine B (Phylo-DAE) | continuous RMSE (holdout species): %.4f\n", rmse_B))
cat(sprintf("Engine B (Phylo-DAE) | categorical accuracy (holdout species): %.4f\n", acc_B))
cat("\nEngine B: example categorical probabilities (first 6 holdout species):\n")
print(round(P_hat_B[which(is_holdout)[1:min(6, sum(is_holdout))], , drop=FALSE], 3))

cat("\n================ SUMMARY (holdout species) ================\n")
cat(sprintf("Engine A (Rphylopars BM) | RMSE %.4f | Acc %.4f\n", rmse_A, acc_A))
cat(sprintf("Engine B (Phylo-DAE)     | RMSE %.4f | Acc %.4f\n", rmse_B, acc_B))
cat("===========================================================\n")