## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4
)

## ----section1-diag, eval = FALSE----------------------------------------------
# library(pigauto)
# data(avonet300, tree300)
# df <- avonet300
# rownames(df) <- df$Species_Key
# df$Species_Key <- NULL
# sum(is.na(df))                  # if 0, there's nothing for impute() to do

## ----section1-fix, eval = FALSE-----------------------------------------------
# set.seed(1L)
# hide   <- sample(which(!is.na(df$Mass)), 30L)
# df_obs <- df
# df_obs$Mass[hide] <- NA          # hide 30 mass values
# 
# result <- impute(df_obs, tree300)
# 
# result$completed$Mass[hide]      # pigauto's imputations
# df$Mass[hide]                    # held-out truth, for comparison
# sum(result$imputed_mask[, "Mass"])  # 30

## ----section2-diag, eval = FALSE----------------------------------------------
# library(pigauto)
# data(avonet300, tree300)
# df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
# 
# table(df$Migration)              # check the marginal distribution
# hide <- sample(which(!is.na(df$Migration)), 30L)
# df$Migration[hide] <- NA
# result <- impute(df, tree300, verbose = FALSE)
# table(result$completed$Migration[hide])

## ----section2-fix, eval = FALSE-----------------------------------------------
# set.seed(1L)
# hide   <- sample(which(!is.na(df$Migration)), 30L)
# df_obs <- df
# df_obs$Migration[hide] <- NA
# 
# # Default settings: prone to majority-class collapse on imbalanced K = 3
# result <- impute(df_obs, tree300, verbose = FALSE)
# table(result$completed$Migration[hide], df$Migration[hide])
# 
# # Recommended for K = 3 ordinal: more draws + mode pooling
# result_mode <- impute(df_obs, tree300, n_imputations = 20L,
#                       pool_method = "mode", verbose = FALSE)
# table(result_mode$completed$Migration[hide], df$Migration[hide])

## ----section3-diag, eval = FALSE----------------------------------------------
# library(pigauto)
# data(avonet300, tree300)
# df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
# fit <- impute(df, tree300, verbose = FALSE)$fit
# 
# # Per-latent-column calibrated gates (since v0.9.1.9002):
# fit$r_cal_bm        # r assigned to the BM baseline
# fit$r_cal_gnn       # r assigned to the GNN delta
# fit$r_cal_mean      # r assigned to the grand mean

## ----section4-diag, eval = FALSE----------------------------------------------
# library(pigauto)
# data(avonet300, tree300)
# df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
# 
# fit <- impute(df, tree300, phylo_signal_gate = TRUE)$fit
# fit$phylo_signal_per_trait

## ----section5-diag, eval = FALSE----------------------------------------------
# # After running impute(), check whether any imputation exceeds the
# # observed maximum by an unrealistic factor:
# predicted_mass <- result$completed$Mass[result$imputed_mask[, "Mass"]]
# obs_max <- max(df$Mass, na.rm = TRUE)
# sum(predicted_mass > 5 * obs_max)

## ----section5-fix, eval = FALSE-----------------------------------------------
# result <- impute(df_obs, tree300,
#                  clamp_outliers = TRUE,
#                  clamp_factor   = 5,    # Tukey-style outlier cap
#                  verbose = FALSE)

