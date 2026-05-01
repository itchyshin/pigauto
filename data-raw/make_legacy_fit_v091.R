# data-raw/make_legacy_fit_v091.R
# One-shot builder for the safety-floor backward-compat test fixture.
# Produces a minimal v0.9.1-format pigauto_fit (no safety-floor slots)
# in inst/extdata/legacy_fit_v091.rds.  Re-run this script only if the
# pigauto_fit class layout changes in a way that requires a new fixture.
#
# The model_state is serialised via save_pigauto() / torch::torch_serialize()
# so the fixture is portable across R sessions.

library(devtools)
devtools::load_all(".")

set.seed(2026L)
data("avonet300", package = "pigauto")
data("tree300",   package = "pigauto")

df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL
df$Mass[sample(300, 30)] <- NA_real_

# Use safety_floor = FALSE to match v0.9.1 behaviour exactly, then strip
# the new safety-floor slots to simulate a pre-v0.9.1.9002 fit.
res <- pigauto::impute(df, tree300,
                       safety_floor  = FALSE,
                       epochs        = 30L,
                       n_imputations = 1L,
                       verbose       = FALSE,
                       seed          = 2026L)

fit_v091 <- res$fit
fit_v091$r_cal_bm              <- NULL
fit_v091$r_cal_gnn             <- NULL
fit_v091$r_cal_mean            <- NULL
fit_v091$mean_baseline_per_col <- NULL
fit_v091$safety_floor          <- NULL

# Serialise torch model_state so the fixture is portable across sessions.
# save_pigauto() converts state_dict tensors to raw via torch_serialize().
out_path <- file.path("inst", "extdata", "legacy_fit_v091.rds")
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
save_pigauto(fit_v091, out_path, compress = TRUE)

cat("Wrote fixture:", out_path, "\n")
cat("fit class:", class(fit_v091), "\n")
cat("fit slots:", paste(names(fit_v091), collapse = ", "), "\n")
cat("fixture size:", format(file.size(out_path), big.mark = ","), "bytes\n")
