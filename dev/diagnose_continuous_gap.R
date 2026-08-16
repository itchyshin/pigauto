# Diagnose: why does raw Rphylopars beat pigauto's continuous output on AVONET300?
# (bench_external_comparators.R, 5 seeds: rphylopars wins all 4 continuous traits)
#
# Decomposition on ONE fixed mask (seed 2027 = that bench's rep 1), all arms
# scored with the bench's exact metric (z-RMSE, train-visible mean/sd):
#
#   A  pigauto full 7-trait fit      -> blend z-RMSE  AND  baseline z-RMSE
#   B  phylopars standalone, 4 cont traits, all observed cells      (bench arm)
#   C  phylopars standalone, 4 cont traits, pigauto's val/test ALSO masked
#   D  pigauto continuous-only 4-trait fit -> baseline + blend z-RMSE
#
# Reads:
#   gate effect        = A.blend    - A.baseline   (does the GNN/gate hurt?)
#   mixed-type path    = A.baseline - D.baseline   (threshold-joint cost on cont cols)
#   calibration tax    = C          - B            (fewer visible cells, same solver)
#   solver wrapper     = D.baseline - C            (pigauto joint machinery vs raw)
#   total              = A.blend    - B
#
# Single mask -> directions and rough magnitudes, not intervals.

suppressPackageStartupMessages({ library(pigauto); library(ape) })
try(torch::torch_set_num_threads(1L), silent = TRUE)
t0 <- proc.time()[["elapsed"]]
say <- function(...) { cat(sprintf("[%7.1fs] ", proc.time()[["elapsed"]] - t0),
                           ..., "\n", sep = ""); flush.console() }

SEED <- 2027L; MISS_FRAC <- 0.30
CONT <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

e <- new.env(parent = emptyenv())
utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300",   package = "pigauto", envir = e)
df <- e$avonet300; tree <- e$tree300
rownames(df) <- df$Species_Key; df$Species_Key <- NULL
stopifnot(all(rownames(df) == tree$tip.label))

# identical mask construction to bench_external_comparators.R
set.seed(SEED)
mask <- matrix(FALSE, nrow(df), ncol(df), dimnames = list(rownames(df), names(df)))
for (v in names(df)) {
  obs_idx <- which(!is.na(df[[v]]))
  mask[sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC)), v] <- TRUE
}
df_miss <- df
for (v in names(df)) df_miss[[v]][mask[, v]] <- NA
df_miss4 <- df_miss[, CONT]

train_stats <- lapply(CONT, function(v) {
  obs <- df[[v]][!mask[, v]]
  c(mean = mean(obs, na.rm = TRUE), sd = stats::sd(obs, na.rm = TRUE))
})
names(train_stats) <- CONT

zrmse <- function(pred_df) vapply(CONT, function(v) {
  idx <- which(mask[, v]); s <- train_stats[[v]]
  sqrt(mean(((as.numeric(pred_df[[v]][idx]) - df[[v]][idx]) / s["sd"])^2,
            na.rm = TRUE))
}, numeric(1))

decode_baseline <- function(res, traits) {
  # decode res$baseline$mu for the named continuous traits to raw scale
  tm_all <- res$fit$trait_map
  out <- as.data.frame(matrix(NA_real_, nrow(df), length(traits),
                              dimnames = list(rownames(df), traits)))
  sp <- rownames(res$data$X_scaled)
  for (v in traits) {
    tm <- tm_all[[v]]
    m <- res$baseline$mu[, tm$latent_cols[1]] * tm$sd + tm$mean
    if (isTRUE(tm$log_transform)) m <- expm1(pmax(m, -20))
    out[sp, v] <- m
  }
  out
}

# ---- A: pigauto full 7-trait ------------------------------------------------
say("A: pigauto full 7-trait fit ...")
resA <- pigauto::impute(df_miss, tree, seed = SEED, verbose = FALSE,
                        n_imputations = 1L)
A_blend <- zrmse(resA$completed)
A_base  <- zrmse(decode_baseline(resA, CONT))
gates <- resA$fit$r_cal_gnn
say(sprintf("A done. r_cal_gnn range: %.3f..%.3f",
            min(gates, na.rm = TRUE), max(gates, na.rm = TRUE)))

# which dispatch path fired? Re-derive the conditions.
tm_all <- resA$fit$trait_map
types <- vapply(tm_all, `[[`, character(1), "type")
n_binord <- sum(types %in% c("binary", "ordinal"))
say(sprintf("trait types: %s", paste(types, collapse = ", ")))
say(sprintf("dispatch conditions: bin+ord cols = %d, Rphylopars = %s => %s",
            n_binord, requireNamespace("Rphylopars", quietly = TRUE),
            if (n_binord >= 1) "THRESHOLD-JOINT path" else "continuous joint"))
# log-transform flags matter for decode comparability
say(sprintf("log_transform flags (cont): %s", paste(vapply(CONT, function(v)
  paste0(v, "=", isTRUE(tm_all[[v]]$log_transform)), character(1)), collapse = ", ")))

# pigauto's additional val/test masking, mapped back to (species, trait)
spl <- resA$fit$splits
X <- resA$data$X_scaled; nr <- nrow(X)
sp_names <- rownames(X)
vt_idx <- c(spl$val_idx, spl$test_idx)
vt_col <- ((vt_idx - 1L) %/% nr) + 1L
vt_row <- ((vt_idx - 1L) %% nr) + 1L
df_missC <- df_miss4
n_extra <- 0L
for (v in CONT) {
  lc <- tm_all[[v]]$latent_cols[1]
  hit <- vt_row[vt_col == lc]
  if (length(hit)) {
    df_missC[sp_names[hit], v] <- NA
    n_extra <- n_extra + length(hit)
  }
}
say(sprintf("C mask: %d extra cells hidden (pigauto val+test) on top of %d genuine",
            n_extra, sum(mask[, CONT])))

# ---- B: phylopars standalone, all observed ----------------------------------
say("B: phylopars standalone (bench arm) ...")
fitB <- Rphylopars::phylopars(
  data.frame(species = rownames(df_miss4), df_miss4, stringsAsFactors = FALSE),
  tree = tree, model = "BM",
  phylo_correlated = TRUE, pheno_correlated = TRUE, REML = TRUE)
B_pred <- as.data.frame(fitB$anc_recon[rownames(df), CONT, drop = FALSE])
B_z <- zrmse(B_pred)
say("B done.")

# ---- C: phylopars standalone, pigauto-visible cells only --------------------
say("C: phylopars standalone with pigauto's val/test also hidden ...")
fitC <- Rphylopars::phylopars(
  data.frame(species = rownames(df_missC), df_missC, stringsAsFactors = FALSE),
  tree = tree, model = "BM",
  phylo_correlated = TRUE, pheno_correlated = TRUE, REML = TRUE)
C_pred <- as.data.frame(fitC$anc_recon[rownames(df), CONT, drop = FALSE])
C_z <- zrmse(C_pred)
say("C done.")

# ---- D: pigauto continuous-only 4-trait -------------------------------------
say("D: pigauto continuous-only 4-trait fit ...")
resD <- pigauto::impute(df_miss4, tree, seed = SEED, verbose = FALSE,
                        n_imputations = 1L)
D_blend <- zrmse(resD$completed)
D_base  <- zrmse(decode_baseline(resD, CONT))
say("D done.")

# ---- report -----------------------------------------------------------------
tab <- rbind(A_blend = A_blend, A_base = A_base, D_blend = D_blend,
             D_base = D_base, C_taxed_phylopars = C_z, B_phylopars = B_z)
say("z-RMSE table:")
print(round(tab, 3))
say("Decomposition (positive = that layer HURTS):")
dec <- rbind(gate_effect    = A_blend - A_base,
             mixedtype_path = A_base - D_base,
             calib_tax      = C_z - B_z,
             solver_wrapper = D_base - C_z,
             TOTAL          = A_blend - B_z)
print(round(dec, 3))
saveRDS(list(tab = tab, dec = dec, gates = gates, seed = SEED,
             n_extra_masked = n_extra),
        "dev/diagnose_continuous_gap.rds")
say("wrote dev/diagnose_continuous_gap.rds")
