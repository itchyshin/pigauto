# PanTHERIA: log_transform TRUE (pigauto default, imputes log(log)) vs FALSE (imputes on the analysis' log scale).
# Same masks, same pigauto default residual_prior "sep"; only log_transform differs.
a <- commandArgs(TRUE); dT <- a[1]; dF <- a[2]
k <- c("dataset","arm","seed","trait")
cT <- read.csv(file.path(dT,"coverage_table.csv")); cF <- read.csv(file.path(dF,"coverage_table.csv"))
cT <- cT[cT$dataset == "pantheria", ]
m <- merge(cT[,c(k,"model_coverage","model_mean_width","model_median_width","split_coverage")], cF[,c(k,"model_coverage","model_mean_width","model_median_width")], by=k, suffixes=c(".T",".F"))
cat("trait-cells:", nrow(m), "\n")
print(aggregate(cbind(model_coverage.T, model_coverage.F, split_coverage) ~ trait, m, function(x) round(mean(x),3)))
cat("mean coverage T", round(mean(m$model_coverage.T),3), "F", round(mean(m$model_coverage.F),3), "split", round(mean(m$split_coverage),3), "\n")
cat("range T", round(range(m$model_coverage.T),3), "F", round(range(m$model_coverage.F),3), "\n")
cat("model below split: T", sum(m$model_coverage.T < m$split_coverage), " F", sum(m$model_coverage.F < m$split_coverage), "of", nrow(m), "\n")
print(aggregate(cbind(mean_w = model_mean_width.F/model_mean_width.T, med_w = model_median_width.F/model_median_width.T) ~ trait, m, function(x) round(mean(x),3)))
k2 <- c("dataset","arm","seed","response","predictor")
sT <- read.csv(file.path(dT,"slope_table.csv")); sF <- read.csv(file.path(dF,"slope_table.csv")); sT <- sT[sT$dataset == "pantheria", ]
s <- merge(sT[,c(k2,"ref_slope","ref_se","mi_slope","mi_se","diff_ref_se","within_5pct")], sF[,c(k2,"mi_slope","mi_se","diff_ref_se","within_5pct")], by=k2, suffixes=c(".T",".F"))
s$mask <- paste(s$arm, s$seed %% 100)
print(s[order(s$response, s$arm, s$seed), c("response","predictor","mask","ref_slope","mi_slope.T","mi_slope.F","diff_ref_se.T","diff_ref_se.F","mi_se.T","mi_se.F")], digits=3, row.names=FALSE)
cat("within 5%: T", sum(s$within_5pct.T), " F", sum(s$within_5pct.F), "of", nrow(s), "\n")
cat("mean |diff| in ref SEs: T", round(mean(abs(s$diff_ref_se.T)),2), " F", round(mean(abs(s$diff_ref_se.F)),2), "; max T", round(max(abs(s$diff_ref_se.T)),2), " F", round(max(abs(s$diff_ref_se.F)),2), "\n")
print(aggregate(cbind(T = diff_ref_se.T, F = diff_ref_se.F) ~ response + arm, s, function(x) round(range(x),2)))
