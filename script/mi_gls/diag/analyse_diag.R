#!/usr/bin/env Rscript
# script/mi_gls/diag/analyse_diag.R
#
# Turns the collected diagnosis CSVs into the tables of
# docs/dev-log/mi-posterior/diagnosis.md. Every number in that document is
# printed by this script into evidence/diagnosis/diag_tables.md.
#
# Usage (cwd = worktree root):
#   Rscript script/mi_gls/diag/analyse_diag.R docs/dev-log/mi-posterior/evidence/diagnosis

args <- commandArgs(trailingOnly = TRUE)
ev <- args[[1L]]
camp <- utils::read.csv(file.path(ev, "campaign_cells.csv"))
orc  <- if (file.exists(file.path(ev, "oracle_long.csv"))) utils::read.csv(file.path(ev, "oracle_long.csv")) else NULL
rr   <- if (file.exists(file.path(ev, "rerun_summary.csv"))) utils::read.csv(file.path(ev, "rerun_summary.csv"), check.names = FALSE) else NULL
rdl  <- if (file.exists(file.path(ev, "rerun_diag_long.csv"))) utils::read.csv(file.path(ev, "rerun_diag_long.csv")) else NULL

out <- character(0)
say <- function(...) out <<- c(out, paste0(...))
f3 <- function(x) formatC(x, format = "f", digits = 3)
f4 <- function(x) formatC(x, format = "f", digits = 4)
mb <- function(d) c(mean = mean(d, na.rm = TRUE),
                    mcse = stats::sd(d, na.rm = TRUE) / sqrt(sum(is.finite(d))),
                    n = sum(is.finite(d)))
cell <- function(d) { z <- mb(d); sprintf("%s (%s)", f4(z[["mean"]]), f4(z[["mcse"]])) }

reg_lab <- c("1" = "1: lambda 1, n 300, MCAR, x", "3" = "3: lambda 1, n 1000, MCAR, x",
             "5" = "5: lambda 1, n 300, MAR, x", "7" = "7: lambda 1, n 1000, MAR, x",
             "9" = "9: lambda 1, n 300, MCAR, both", "11" = "11: lambda 1, n 1000, MCAR, both",
             "13" = "13: lambda 1, n 300, MAR, both", "15" = "15: lambda 1, n 1000, MAR, both",
             "2" = "2: lambda 0.5, n 300, MCAR, x", "4" = "4: lambda 0.5, n 1000, MCAR, x")

# ---- Table 1: oracle arms vs campaign, paired bias --------------------------------
if (!is.null(orc)) {
  say("## Table 1. Paired slope bias (imputed minus complete data, same rep); mean (MCSE)")
  say("")
  say("Campaign = posterior_full from the campaign output, restricted to the SAME reps as the oracle arms.")
  say("")
  say("| regime | analysis | reps | campaign posterior | oracle_true (DGP covariance) | oracle_kl (sampler model, KL pseudo-true) | oracle_nominal (sampler model, nominal) | campaign minus oracle_kl | oracle_kl minus oracle_true |")
  say("|---|---|---|---|---|---|---|---|---|")
  for (r in c(1, 3, 5, 7, 9, 11, 13, 15, 2, 4)) {
    for (ds in c("gls", "phylolm")) {
      s <- orc[orc$regime_id == r & orc$downstream == ds, ]
      if (!nrow(s)) next
      comp <- s[s$arm == "complete", c("rep", "estimate")]
      pd <- function(a) {
        z <- merge(comp, s[s$arm == a, c("rep", "estimate")], by = "rep", suffixes = c(".c", ".a"))
        z <- z[order(z$rep), ]
        z$estimate.a - z$estimate.c
      }
      cc <- camp[camp$regime_id == r & camp$rep %in% comp$rep, ]
      # complete-data slopes recomputed here must equal the campaign's
      chk <- merge(comp, cc[, c("rep", paste0("complete_", ds))], by = "rep")
      stopifnot(max(abs(chk$estimate - chk[[paste0("complete_", ds)]])) < 1e-8)
      cc <- cc[order(cc$rep), ]
      stopifnot(identical(cc$rep, sort(comp$rep)))
      campd <- cc[[paste0("full_", ds)]] - cc[[paste0("complete_", ds)]]
      say(sprintf("| %s | %s | %d | %s | %s | %s | %s | %s | %s |", reg_lab[[as.character(r)]], ds,
                  nrow(comp), cell(campd), cell(pd("oracle_true")), cell(pd("oracle_kl")),
                  cell(pd("oracle_nominal")), cell(campd - pd("oracle_kl")),
                  cell(pd("oracle_kl") - pd("oracle_true"))))
    }
  }
  say("")
  say("Check: the complete-data slopes recomputed by oracle_cell.R equal the campaign's to < 1e-8 in every rep used (stopifnot above).")
  say("")
  kl <- unique(orc[, c("regime_id", "rep", "a_star", "b_star", "lambda_star", "mean_diag_V", "min_diag_V", "kl_conv")])
  say("KL pseudo-true parameters of the sampler's model (raw scale): mean over reps")
  say("")
  say("| regime | reps | a* (Sigma_P scale) | b* (Sigma_E scale) | lambda* = a*/(a*+b*) | mean diag(V_sim) | min diag(V_sim) | optim failures |")
  say("|---|---|---|---|---|---|---|---|")
  for (r in sort(unique(kl$regime_id))) {
    k <- kl[kl$regime_id == r, ]
    say(sprintf("| %d | %d | %s | %s | %s | %s | %s | %d |", r, nrow(k), f3(mean(k$a_star)), f4(mean(k$b_star)),
                f3(mean(k$lambda_star)), f3(mean(k$mean_diag_V)), f3(mean(k$min_diag_V)), sum(k$kl_conv != 0)))
  }
  say("")
}

# ---- Table 2: H1, posterior summaries vs truth ------------------------------------
if (!is.null(rr)) {
  h1 <- rr[rr$tag %in% c("h1", "prior", "twin"), ]
  if (nrow(h1)) {
    say("## Table 2. Posterior means per fit (latent z-scale), averaged over reps (SD over reps)")
    say("")
    say("| arm | regime | reps | lambda_x | lambda_y | corr_P | corr_E | beta_P = SP12/SP22 | beta_E = SE12/SE22 | converged |")
    say("|---|---|---|---|---|---|---|---|---|---|")
    msd <- function(v) sprintf("%s (%s)", f3(mean(v)), f3(stats::sd(v)))
    for (tg in c("h1", "prior", "twin")) for (r in c(1, 3)) {
      s <- h1[h1$tag == tg & h1$regime_id == r, ]
      if (!nrow(s)) next
      say(sprintf("| %s | %d | %d | %s | %s | %s | %s | %s | %s | %d/%d |", tg, r, nrow(s),
                  msd(s$pmean_lambda_x), msd(s$pmean_lambda_y), msd(s$pmean_corr_P),
                  msd(s$pmean_corr_E), msd(s$`pmean_beta_P..SP12.SP22.`),
                  msd(s$`pmean_beta_E..SE12.SE22.`), sum(s$converged), nrow(s)))
    }
    say("")
    say("Truth: rho = 0.7 in both components; lambda = 1 (Sigma_E = 0) nominally; the KL pseudo-true lambda* is in Table 1's second table.")
    say("")
    # reproduction check for the h1 re-runs (frozen code, campaign settings)
    s <- rr[rr$tag == "h1", ]
    if (nrow(s)) {
      z <- merge(s, camp, by = c("regime_id", "rep"), suffixes = c("", ".camp"))
      say(sprintf("Reproduction: h1 re-runs (%d fits) vs campaign: max |slope difference| gls %.2e, phylolm %.2e; max |min ESS difference| %.2e.",
                  nrow(z), max(abs(z$full_gls - z$full_gls.camp)), max(abs(z$full_phylolm - z$full_phylolm.camp)),
                  max(abs(z$min_ess - z$min_ess.camp))))
      say("")
    }
    # H3 checks
    say(sprintf("H3 checks over all %d sampler re-runs: any log-transform = %s; max |observed cell change| in completed data = %.2e; any NA in completed data = %s.",
                nrow(rr), any(rr$h3_log_transform_any), max(rr$h3_max_obs_change), any(rr$h3_any_na)))
    say("")
    s1 <- rr[rr$tag == "h1", ]; s2 <- rr[rr$tag == "prior", ]
    if (nrow(s1) && nrow(s2)) {
      z <- merge(s1, s2, by = c("regime_id", "rep"), suffixes = c(".h1", ".prior"))
      say("Prior arm (S_E = 1e-4 diag(obs var) instead of 0.01) minus campaign prior, same reps: MI slope difference, mean (MCSE)")
      say("")
      say("| regime | reps | gls | phylolm | corr_E default | corr_E prior arm | Sigma_E[1,1] default | Sigma_E[1,1] prior arm |")
      say("|---|---|---|---|---|---|---|---|")
      for (r in sort(unique(z$regime_id))) {
        zz <- z[z$regime_id == r, ]
        say(sprintf("| %d | %d | %s | %s | %s | %s | %s | %s |", r, nrow(zz),
                    cell(zz$full_gls.prior - zz$full_gls.h1), cell(zz$full_phylolm.prior - zz$full_phylolm.h1),
                    f3(mean(zz$pmean_corr_E.h1)), f3(mean(zz$pmean_corr_E.prior)),
                    formatC(mean(zz$`pmean_Sigma_E.1.1..h1`), format = "e", digits = 2),
                    formatC(mean(zz$`pmean_Sigma_E.1.1..prior`), format = "e", digits = 2)))
      }
      say("")
    }
    say("## Table 3. Paired slope bias of the sampler: original DGP vs in-model twin (same seeds, trees, noise, masks)")
    say("")
    say("| arm | regime | analysis | reps | paired bias (MCSE) | campaign bias, same reps |")
    say("|---|---|---|---|---|---|")
    for (tg in c("h1", "twin", "prior")) for (r in c(1, 3)) for (ds in c("gls", "phylolm")) {
      s <- rr[rr$tag == tg & rr$regime_id == r, ]
      if (!nrow(s)) next
      cc <- camp[camp$regime_id == r & camp$rep %in% s$rep, ]
      say(sprintf("| %s | %d | %s | %d | %s | %s |", tg, r, ds, nrow(s),
                  cell(s[[paste0("full_", ds)]] - s[[paste0("complete_", ds)]]),
                  cell(cc[[paste0("full_", ds)]] - cc[[paste0("complete_", ds)]])))
    }
    say("")
  }

  # ---- Table 4: F2 -------------------------------------------------------------
  f2 <- rr[grepl("^f2_", rr$tag), ]
  if (nrow(f2)) {
    say("## Table 4. Non-converged campaign fits re-run at 1x (subset), 2x and 4x chain length")
    say("")
    say("| regime | rep | campaign max R-hat | campaign min ESS | 1x min-ESS parameter | 2x max R-hat | 2x min ESS (param) | 4x max R-hat | 4x min ESS (param) | phylolm slope 1x / 2x / 4x | gls slope 1x / 2x / 4x | pooled phylolm SE 1x |")
    say("|---|---|---|---|---|---|---|---|---|---|---|---|")
    keys <- unique(f2[, c("regime_id", "rep")]); keys <- keys[order(keys$regime_id, keys$rep), ]
    for (i in seq_len(nrow(keys))) {
      r <- keys$regime_id[i]; k <- keys$rep[i]
      c0 <- camp[camp$regime_id == r & camp$rep == k, ]
      g <- function(tg) f2[f2$tag == tg & f2$regime_id == r & f2$rep == k, ]
      x1 <- g("f2_x1"); x2 <- g("f2_x2"); x4 <- g("f2_x4")
      fmt <- function(x, col, fn = f3) if (nrow(x)) fn(x[[col]]) else "-"
      essp <- function(x) if (nrow(x)) sprintf("%.0f (%s)", x$min_ess, x$ess_param) else "-"
      say(sprintf("| %d | %d | %s | %.0f | %s | %s | %s | %s | %s | %s / %s / %s | %s / %s / %s | %s |",
                  r, k, f3(c0$max_rhat), c0$min_ess,
                  if (nrow(x1)) sprintf("%s (%.0f)", x1$ess_param, x1$min_ess) else "-",
                  fmt(x2, "max_rhat"), essp(x2), fmt(x4, "max_rhat"), essp(x4),
                  f3(c0$full_phylolm), fmt(x2, "full_phylolm"), fmt(x4, "full_phylolm"),
                  f3(c0$full_gls), fmt(x2, "full_gls"), fmt(x4, "full_gls"),
                  if (nrow(x1)) f3(x1$full_phylolm_se) else "-"))
    }
    say("")
    for (tg in c("f2_x1", "f2_x2", "f2_x4")) {
      s <- f2[f2$tag == tg, ]
      if (!nrow(s)) next
      z <- merge(s, camp, by = c("regime_id", "rep"), suffixes = c("", ".camp"))
      say(sprintf("%s: %d fits; converged %d; max R-hat range %.3f to %.3f; min ESS range %.0f to %.0f; min-ESS parameter counts: %s; |phylolm slope - campaign| mean %.4f max %.4f; |gls slope - campaign| mean %.4f max %.4f; median wall %.0f s.",
                  tg, nrow(s), sum(s$converged), min(s$max_rhat), max(s$max_rhat), min(s$min_ess), max(s$min_ess),
                  paste(names(table(s$ess_param)), table(s$ess_param), sep = "=", collapse = ", "),
                  mean(abs(z$full_phylolm - z$full_phylolm.camp)), max(abs(z$full_phylolm - z$full_phylolm.camp)),
                  mean(abs(z$full_gls - z$full_gls.camp)), max(abs(z$full_gls - z$full_gls.camp)),
                  stats::median(s$wall_s)))
      say("")
    }
    w <- reshape(f2[, c("tag", "regime_id", "rep", "full_gls", "full_phylolm")], idvar = c("regime_id", "rep"),
                 timevar = "tag", direction = "wide")
    w <- merge(w, camp[, c("regime_id", "rep", "full_gls", "full_phylolm")], by = c("regime_id", "rep"))
    for (ds in c("phylolm", "gls")) {
      a1 <- w[[paste0("full_", ds)]]; a2 <- w[[paste0("full_", ds, ".f2_x2")]]; a4 <- w[[paste0("full_", ds, ".f2_x4")]]
      say(sprintf("%s pooled slope, %d reps: signed change 2x - 1x %s, 4x - 1x %s (mean (MCSE)); mean |2x - 1x| %s, mean |4x - 1x| %s, mean |4x - 2x| %s.",
                  ds, nrow(w), cell(a2 - a1), cell(a4 - a1), f4(mean(abs(a2 - a1))), f4(mean(abs(a4 - a1))), f4(mean(abs(a4 - a2)))))
      say("")
    }
    x1se <- f2[f2$tag == "f2_x4", ]
    say(sprintf("Pooled SE of the slope in these reps (4x runs): phylolm median %s (range %s to %s); gls median %s.",
                f3(stats::median(x1se$full_phylolm_se)), f3(min(x1se$full_phylolm_se)), f3(max(x1se$full_phylolm_se)),
                f3(stats::median(x1se$full_gls_se))))
    say("")
    if (!is.null(rdl)) {
      d1 <- rdl[rdl$tag == "f2_x1", ]
      if (nrow(d1)) {
        say("Per-parameter diagnostics of the 1x re-runs (parameters with ESS <= 400 or R-hat >= 1.05):")
        say("")
        bad <- d1[d1$ess_bulk <= 400 | d1$rhat >= 1.05, ]
        say("| regime | rep | parameter | R-hat | bulk ESS |")
        say("|---|---|---|---|---|")
        for (j in seq_len(nrow(bad))) say(sprintf("| %d | %d | %s | %s | %.0f |", bad$regime_id[j], bad$rep[j],
                                                    bad$parameter[j], f3(bad$rhat[j]), bad$ess_bulk[j]))
        say("")
      }
    }
  }
}
writeLines(out, file.path(ev, "diag_tables.md"))
cat(out, sep = "\n")
