## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  eval     = FALSE
)

## ----impute-------------------------------------------------------------------
# library(pigauto)
# data(avonet300, tree300)
# 
# traits <- avonet300[, c("Mass", "Wing.Length",
#                         "Beak.Length_Culmen", "Tarsus.Length")]
# rownames(traits) <- avonet300$Species_Key
# 
# set.seed(1)
# traits$Mass[sample(nrow(traits), 60)] <- NA
# traits$Wing.Length[sample(nrow(traits), 45)] <- NA
# 
# mi <- multi_impute(traits, tree300, m = 20L,
#                    draws_method = "posterior", seed = 1L)
# mi

## ----diagnostics--------------------------------------------------------------
# mi$posterior$diagnostics
# attr(mi$posterior$diagnostics, "converged")

## ----gls----------------------------------------------------------------------
# fits <- with_imputations(mi, function(d) {
#   d$species <- rownames(d)
#   nlme::gls(log(Wing.Length) ~ log(Mass), data = d,
#             correlation = ape::corBrownian(phy = tree300, form = ~species),
#             method = "ML")
# })
# pool_mi(fits)

## ----mixed--------------------------------------------------------------------
# add_diet <- function(d) {
#   d$Trophic.Level <- avonet300$Trophic.Level[
#     match(rownames(d), avonet300$Species_Key)]
#   d
# }
# 
# fits_lmer <- with_imputations(mi, function(d) {
#   lme4::lmer(log(Wing.Length) ~ log(Mass) + (1 | Trophic.Level),
#              data = add_diet(d))
# })
# pool_mi(fits_lmer)
# 
# fits_tmb <- with_imputations(mi, function(d) {
#   glmmTMB::glmmTMB(log(Wing.Length) ~ log(Mass) + (1 | Trophic.Level),
#                    data = add_diet(d))
# })
# pool_mi(fits_tmb)

## ----drmtmb-------------------------------------------------------------------
# fits_drm <- with_imputations(mi, function(d) {
#   drmTMB::drmTMB(drmTMB::bf(log(Wing.Length) ~ log(Mass),
#                             sigma ~ log(Mass)),
#                  data = d)
# })
# pool_mi(fits_drm)

## ----gllvmtmb-----------------------------------------------------------------
# to_long <- function(d) {
#   tr <- names(d)
#   data.frame(species = factor(rep(rownames(d), times = length(tr))),
#              trait   = factor(rep(tr, each = nrow(d)), levels = tr),
#              value   = log(unlist(d, use.names = FALSE)))
# }
# 
# fits_gllvm <- with_imputations(mi, function(d) {
#   gllvmTMB::gllvmTMB(value ~ 0 + trait +
#                        latent(0 + trait | species, d = 1, unique = FALSE),
#                      data = to_long(d), trait = "trait", unit = "species",
#                      family = gaussian(), REML = FALSE,
#                      engine = "tmb", silent = TRUE)
# })
# pool_mi(fits_gllvm)

## ----brms---------------------------------------------------------------------
# fit_b <- brms::brm_multiple(log(Wing.Length) ~ log(Mass),
#                             data = mi$datasets, chains = 2, refresh = 0)
# summary(fit_b)

## ----mcmcglmm-----------------------------------------------------------------
# fits_mc <- lapply(mi$datasets, function(d) {
#   MCMCglmm::MCMCglmm(log(Wing.Length) ~ log(Mass), data = d, verbose = FALSE)
# })
# sol <- coda::as.mcmc(do.call(rbind, lapply(fits_mc, function(f) f$Sol)))
# summary(sol)

