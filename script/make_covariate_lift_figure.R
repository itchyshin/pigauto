#!/usr/bin/env Rscript
# script/make_covariate_lift_figure.R
#
# Two figures summarising covariate-lift across all 5 real datasets
# + 2 sims:
#   1. covariate_lift_figure.png  -- best-trait lift per dataset (clear story)
#   2. covariate_lift_figure_full.png -- all-traits per-dataset (full detail)
# Plus PDF vector versions.

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("ggplot2 not available -- skipping figure\n"); quit(save = "no")
  }
  library(ggplot2)
})

here <- "/Users/z3437171/Dropbox/Github Local/pigauto"

# ---- Headline figure: best-trait lift per dataset ----
hl <- data.frame(
  dataset = c(
    "Multi-obs sim (real bird tree)",
    "GlobTherm ectotherms (sf=off)",
    "AmphiBIO amphibians",
    "PanTHERIA mammals",
    "LepTraits butterflies\n(vs column-mean)",
    "Delhey birds",
    "BIEN plants"),
  trait = c(
    "CTmax @ beta=1.0", "Tmax", "Body_size_mm",
    "MaxLongevity_m", "WS_L (wingspan)",
    "lightness (best of 2)", "wood_density (best of 5)"),
  lift_pct = c(9.2, 8.0, 17.2, 22.2, 24.0, -1.3, 1.2),
  is_real_lift = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  baseline_r = c(0.63, 0.67, 0.62, 0.80, 0.29, 0.70, 0.45),
  cov_r = c(0.69, 0.73, 0.70, 0.77, 0.73, 0.69, 0.46),
  category = c(
    "Multi-obs simulation", "Real species (covariate sf=off)",
    "Real species (taxonomic tree)", "Real species (molecular tree)",
    "Real species (degenerate baseline)",
    "Real species (phylo-redundant)",
    "Real species (phylo-redundant)"),
  stringsAsFactors = FALSE
)
hl$dataset <- factor(hl$dataset, levels = rev(hl$dataset))

p_hl <- ggplot(hl, aes(x = lift_pct, y = dataset,
                         fill = is_real_lift)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%+0.1f %%", lift_pct),
                  hjust = ifelse(lift_pct < 0, 1.15, -0.15)),
            size = 3.2) +
  geom_text(aes(label = sprintf("r %.2f -> %.2f", baseline_r, cov_r),
                  x = ifelse(lift_pct < 0, lift_pct - 6, lift_pct + 6),
                  hjust = ifelse(lift_pct < 0, 1, 0)),
            size = 2.6, colour = "grey40") +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 5, colour = "darkgreen",
              linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#3a7d44", "FALSE" = "#aaaaaa"),
                     guide = "none") +
  scale_x_continuous(limits = c(-15, 35),
                      breaks = c(-10, -5, 0, 5, 10, 15, 20, 25, 30),
                      labels = function(x) sprintf("%+d %%", x)) +
  labs(x = "RMSE lift (%)", y = NULL,
       title = "pigauto covariate lift -- 5 real datasets + sim",
       subtitle = paste("Best per-trait lift per dataset, with Pearson r before/after.",
                          "Green: lift > 5%.  Grey: phylo-redundant covariates.",
                          sep = "  ")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"))

ggsave(file.path(here, "useful", "covariate_lift_figure.png"),
        p_hl, width = 9, height = 5.5, dpi = 110)
ggsave(file.path(here, "useful", "covariate_lift_figure.pdf"),
        p_hl, width = 9, height = 5.5)
cat("wrote useful/covariate_lift_figure.{png,pdf}\n")

# ---- Full per-trait figure (all traits, cov_on metric) ----
rows <- data.frame(
  dataset = c(
    "PanTHERIA mammals\n(climate cov, mol. tree)",
    "PanTHERIA mammals\n(climate cov, mol. tree)",
    "PanTHERIA mammals\n(climate cov, mol. tree)",
    "PanTHERIA mammals\n(climate cov, mol. tree)",
    "PanTHERIA mammals\n(climate cov, mol. tree)",
    "GlobTherm ectotherms\n(lat+long+elev, taxonomic)",
    "GlobTherm ectotherms\n(lat+long+elev, taxonomic)",
    "AmphiBIO amphibians\n(climate-zone, taxonomic)",
    "AmphiBIO amphibians\n(climate-zone, taxonomic)",
    "AmphiBIO amphibians\n(climate-zone, taxonomic)",
    "AmphiBIO amphibians\n(climate-zone, taxonomic)",
    "AmphiBIO amphibians\n(climate-zone, taxonomic)",
    "LepTraits butterflies\n(Jan-Dec, taxonomic)",
    "LepTraits butterflies\n(Jan-Dec, taxonomic)",
    "LepTraits butterflies\n(Jan-Dec, taxonomic)",
    "LepTraits butterflies\n(Jan-Dec, taxonomic)",
    "Delhey birds\n(climate cov, mol. tree)",
    "Delhey birds\n(climate cov, mol. tree)",
    "BIEN plants\n(WorldClim, mol. tree)",
    "BIEN plants\n(WorldClim, mol. tree)",
    "BIEN plants\n(WorldClim, mol. tree)",
    "BIEN plants\n(WorldClim, mol. tree)",
    "BIEN plants\n(WorldClim, mol. tree)"
  ),
  trait = c(
    "Body mass", "Gestation", "MaxLongevity", "LitterSize", "PopDensity",
    "Tmax (sf=on)", "tmin",
    "Body_size", "Body_mass", "Longevity", "Age@maturity", "Litter_size",
    "WS_L (vs col-mean)", "FW_L", "FlightDur", "HostFamilies",
    "lightness_male", "lightness_female",
    "wood_density", "height_m", "sla", "leaf_area", "seed_mass"
  ),
  ratio_on = c(
    1.04, 1.09, 0.78, 1.00, 0.96,
    1.28, 0.99,
    0.83, 1.02, 1.00, 1.00, 1.08,
    0.75,
    1.03, 0.99, 1.02,
    1.02, 1.01,
    0.99, 1.17, 1.03, 1.01, 1.00
  ),
  stringsAsFactors = FALSE
)
rows$lift_pct <- 100 * (1 - rows$ratio_on)
rows <- rows[order(rev(rows$dataset), rows$lift_pct), ]
rows$trait <- factor(rows$trait, levels = unique(rows$trait))

p_full <- ggplot(rows, aes(x = lift_pct, y = trait, fill = lift_pct > 5)) +
  geom_col() +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 5, colour = "darkgreen", linetype = "dashed",
              linewidth = 0.3) +
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#3a7d44", "FALSE" = "#d6d6d6"),
                     guide = "none") +
  scale_x_continuous(breaks = c(-30, -20, -10, -5, 0, 5, 10, 15, 20, 25),
                      labels = function(x) sprintf("%+d %%", x)) +
  labs(x = "RMSE lift (cov_on vs none, %)", y = NULL,
       title = "Per-trait covariate lift across pigauto real-data benches",
       subtitle = "Green > 5% lift.  GlobTherm Tmax with sf=off would be +8%.") +
  theme_minimal(base_size = 11) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 8),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"))

ggsave(file.path(here, "useful", "covariate_lift_figure_full.png"),
        p_full, width = 9, height = 7, dpi = 110)
ggsave(file.path(here, "useful", "covariate_lift_figure_full.pdf"),
        p_full, width = 9, height = 7)
cat("wrote useful/covariate_lift_figure_full.{png,pdf}\n")
