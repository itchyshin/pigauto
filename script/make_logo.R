#!/usr/bin/env Rscript
# script/make_logo.R ── pigauto hex sticker
#
# "pigauto" = pig + auto (automobile).
# A pink pig in driving goggles pilots an open-top cartoon car.
# The pig is large and prominent (fr = 0.28), with ears peeking above the
# car body.  A small 4-tip gold cladogram floats above the pig's head.
# Gold colour unifies: wheels (NN nodes) + cladogram tips (phylo leaves).
#
# Run from the package root:
#   Rscript script/make_logo.R
#
# Writes: man/figures/logo_v{1..4}.png
#         man/figures/logo.png   (copy of v1 if not already present)
#         pkgdown/assets/logo_v{1..4}.png   (same images for gallery page)
#
# Required:
#   install.packages(c("hexSticker", "ggplot2", "ggforce"))
# Optional (nicer font):
#   install.packages(c("showtext", "sysfonts"))

for (pkg in c("hexSticker", "ggplot2", "ggforce")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Please install '", pkg, "': install.packages('", pkg, "')")
}
suppressPackageStartupMessages({
  library(hexSticker)
  library(ggplot2)
  library(ggforce)
})

# ── optional font ─────────────────────────────────────────────────────────────
ttf <- tryCatch({
  if (requireNamespace("showtext", quietly = TRUE) &&
      requireNamespace("sysfonts",  quietly = TRUE)) {
    sysfonts::font_add_google("Fredoka One", "fredoka")
    showtext::showtext_auto()
    "fredoka"
  } else "sans"
}, error = function(e) "sans")

# ── 4 colour palettes ─────────────────────────────────────────────────────────
# Each palette: bg (hex fill), border (hex edge + text), car (body),
# car_d (car outline/dark), glass (window tint), node (wheels + NN gold),
# wh_d (wheel hub dark), pig (face pink), pig_d (face outline),
# snout (snout pink), goggle (lens amber)
PALETTES <- list(
  v1 = list(        # sky navy + warm orange  — the "classic"
    name   = "Classic (navy / orange)",
    bg     = "#1a3a5c", border = "#ff6f3c",
    car    = "#2e6fa3", car_d  = "#0d3352",
    glass  = "#a8d8ea", node   = "#f7c948", wh_d = "#0d3352",
    pig    = "#f4a0b0", pig_d  = "#c45c72", snout = "#e8788e",
    goggle = "#ffeaa7"
  ),
  v2 = list(        # mid purple + gold  — "regal"
    name   = "Regal (purple / gold)",
    bg     = "#3d2d6e", border = "#f7c948",
    car    = "#5a3ea0", car_d  = "#2a1a50",
    glass  = "#8ecae6", node   = "#f7c948", wh_d = "#2a1a50",
    pig    = "#f4a0b0", pig_d  = "#c45c72", snout = "#e8788e",
    goggle = "#ffeaa7"
  ),
  v3 = list(        # mid forest green + mint  — "ecological"
    name   = "Ecological (forest / mint)",
    bg     = "#1e5c40", border = "#95d5b2",
    car    = "#3a8a5f", car_d  = "#1b4332",
    glass  = "#b7e4c7", node   = "#f7c948", wh_d = "#0d2b1c",
    pig    = "#f4a0b0", pig_d  = "#c45c72", snout = "#e8788e",
    goggle = "#ffeaa7"
  ),
  v4 = list(        # slate + hot pink  — "punchy"
    name   = "Punchy (slate / hot-pink)",
    bg     = "#454560", border = "#ff4d6d",
    car    = "#6264a0", car_d  = "#2e2e5a",
    glass  = "#adb5bd", node   = "#f7c948", wh_d = "#2e2e5a",
    pig    = "#f4a0b0", pig_d  = "#c45c72", snout = "#e8788e",
    goggle = "#ffeaa7"
  )
)

# ── subplot: cartoon convertible car with large pig driver ────────────────────
# Coordinate space: [-1, 1] × [-1, 1].
# Car faces RIGHT.  Open-top convertible so pig's head pokes out prominently.
# No roof polygon — pig face is big (fr = 0.28) and fully visible above the
# car body.  No cladogram anywhere.
make_subplot <- function(pal) {

  fc <- c(-0.06, 0.30)   # pig face centre (raised a bit for open-top design)
  fr <- 0.28             # pig face radius — big and prominent

  # Pig ears (two small circles above/behind the face)
  ear_df <- data.frame(
    x = c(fc[1] - 0.195, fc[1] + 0.195),
    y = c(fc[2] + 0.215,  fc[2] + 0.215),
    r = c(0.085, 0.085)
  )
  ear_inner_df <- data.frame(
    x = c(fc[1] - 0.195, fc[1] + 0.195),
    y = c(fc[2] + 0.215,  fc[2] + 0.215),
    r = c(0.050, 0.050)
  )

  # Smile arc
  theta_sm <- seq(238, 302, length.out = 22) * pi / 180
  smile_df <- data.frame(
    x = fc[1] + 0.082 * cos(theta_sm),
    y = (fc[2] - 0.155) + 0.082 * sin(theta_sm)
  )

  # Speed lines (motion blur, left of car)
  speed_df <- data.frame(
    x    = -0.76,
    xend = c(-0.93, -0.89, -0.91, -0.87),
    y    = c(-0.08, -0.02,  0.04,  0.08),
    yend = c(-0.08, -0.02,  0.04,  0.08)
  )

  # Exhaust NN nodes (gold puffs)
  exh_df <- data.frame(
    x = c(-0.80, -0.84, -0.82),
    y = c(-0.086, -0.070, -0.055),
    r = c( 0.022,  0.018,  0.014)
  )

  # Wheels
  wh_df  <- data.frame(x = c(-0.45,  0.43), y = c(-0.28, -0.28), r = 0.180)
  hub_df <- data.frame(x = c(-0.45,  0.43), y = c(-0.28, -0.28), r = 0.070)
  dot_df <- data.frame(x = c(-0.45,  0.43), y = c(-0.28, -0.28), r = 0.025)

  ggplot() +

    # ── speed lines (behind everything) ──────────────────────────────────────
    geom_segment(data = speed_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 colour = pal$node, linewidth = 0.5, alpha = 0.55) +

    # ── wheels (behind car body) ──────────────────────────────────────────────
    geom_circle(data = wh_df,  aes(x0 = x, y0 = y, r = r),
                fill = pal$node, colour = pal$wh_d, linewidth = 0.7) +
    geom_circle(data = hub_df, aes(x0 = x, y0 = y, r = r),
                fill = pal$wh_d, colour = NA) +
    geom_circle(data = dot_df, aes(x0 = x, y0 = y, r = r),
                fill = pal$node, colour = NA) +

    # ── pig ears (behind face) ────────────────────────────────────────────────
    geom_circle(data = ear_df,
                aes(x0 = x, y0 = y, r = r),
                fill = pal$pig, colour = pal$pig_d, linewidth = 0.55) +
    geom_circle(data = ear_inner_df,
                aes(x0 = x, y0 = y, r = r),
                fill = pal$snout, colour = NA) +

    # ── pig face ──────────────────────────────────────────────────────────────
    geom_circle(data = data.frame(x = fc[1], y = fc[2], r = fr),
                aes(x0 = x, y0 = y, r = r),
                fill = pal$pig, colour = pal$pig_d, linewidth = 0.8) +

    # ── driving goggles (left + right) ───────────────────────────────────────
    geom_circle(data = data.frame(x = fc[1] - 0.127, y = fc[2] + 0.018,
                                  r = 0.098),
                aes(x0 = x, y0 = y, r = r),
                fill = pal$goggle, colour = "#2a1800", linewidth = 1.3) +
    geom_circle(data = data.frame(x = fc[1] + 0.127, y = fc[2] + 0.018,
                                  r = 0.098),
                aes(x0 = x, y0 = y, r = r),
                fill = pal$goggle, colour = "#2a1800", linewidth = 1.3) +
    geom_segment(data = data.frame(
                   x    = fc[1] - 0.029, xend = fc[1] + 0.029,
                   y    = fc[2] + 0.018, yend = fc[2] + 0.018),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "#2a1800", linewidth = 1.3) +
    geom_segment(data = data.frame(
                   x    = fc[1] - 0.225, xend = fc[1] - fr * 0.80,
                   y    = fc[2] + 0.018, yend = fc[2] + 0.018),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "#2a1800", linewidth = 0.85) +
    geom_segment(data = data.frame(
                   x    = fc[1] + 0.225, xend = fc[1] + fr * 0.80,
                   y    = fc[2] + 0.018, yend = fc[2] + 0.018),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "#2a1800", linewidth = 0.85) +

    # ── pupils + highlights ───────────────────────────────────────────────────
    geom_circle(data = data.frame(x = c(fc[1] - 0.127, fc[1] + 0.127),
                                  y = c(fc[2] + 0.018, fc[2] + 0.018),
                                  r = 0.040),
                aes(x0 = x, y0 = y, r = r), fill = "#1a1a1a", colour = NA) +
    geom_circle(data = data.frame(x = c(fc[1] - 0.098, fc[1] + 0.156),
                                  y = c(fc[2] + 0.040, fc[2] + 0.040),
                                  r = 0.018),
                aes(x0 = x, y0 = y, r = r), fill = "white", colour = NA) +

    # ── snout + nostrils ──────────────────────────────────────────────────────
    geom_circle(data = data.frame(x = fc[1], y = fc[2] - 0.140, r = 0.118),
                aes(x0 = x, y0 = y, r = r),
                fill = pal$snout, colour = pal$pig_d, linewidth = 0.55) +
    geom_circle(data = data.frame(x = c(fc[1] - 0.047, fc[1] + 0.047),
                                  y = c(fc[2] - 0.140, fc[2] - 0.140),
                                  r = 0.036),
                aes(x0 = x, y0 = y, r = r),
                fill = "#4a1a1a", colour = NA, alpha = 0.70) +

    # ── smile ─────────────────────────────────────────────────────────────────
    geom_path(data = smile_df, aes(x = x, y = y),
              colour = pal$pig_d, linewidth = 0.75) +

    # ── car body (covers lower part of pig) ───────────────────────────────────
    geom_rect(
      data = data.frame(xmin = -0.72, xmax = 0.72,
                        ymin = -0.22, ymax =  0.12),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = pal$car, colour = pal$car_d, linewidth = 0.6) +

    # ── windshield strip (low glass bar across the car body top) ─────────────
    geom_rect(
      data = data.frame(xmin = -0.38, xmax = 0.35,
                        ymin =  0.07, ymax =  0.13),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = pal$glass, colour = pal$car_d, linewidth = 0.4, alpha = 0.30) +

    # ── front bumper ──────────────────────────────────────────────────────────
    geom_rect(
      data = data.frame(xmin = 0.72, xmax = 0.76,
                        ymin = -0.12, ymax = 0.06),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = pal$car_d, colour = NA) +

    # ── headlight (gold rectangle) ────────────────────────────────────────────
    geom_rect(
      data = data.frame(xmin = 0.67, xmax = 0.75,
                        ymin =  0.01, ymax = 0.08),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = pal$node, colour = NA) +

    # ── rear bumper ───────────────────────────────────────────────────────────
    geom_rect(
      data = data.frame(xmin = -0.76, xmax = -0.72,
                        ymin = -0.12,  ymax =  0.06),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = pal$car_d, colour = NA) +

    # ── exhaust pipe ──────────────────────────────────────────────────────────
    geom_rect(
      data = data.frame(xmin = -0.77, xmax = -0.67,
                        ymin = -0.11,  ymax = -0.07),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = pal$car_d, colour = NA) +

    # ── exhaust NN nodes (gold) ───────────────────────────────────────────────
    geom_circle(data = exh_df,
                aes(x0 = x, y0 = y, r = r),
                fill = pal$node, colour = NA, alpha = 0.75) +

    # ── mini cladogram above pig head (4 tips, symmetric) ────────────────────
    # Root at (cx, 0.635); tips at y = 0.835.  Gold colour = leaf/node motif.
    {
      cx <- fc[1]   # same horizontal centre as pig face = -0.06
      tr_segs <- data.frame(
        #           1-stem       2-root-H       3-Lnode-V    4-Lfork-H
        x    = c(cx,       cx-0.18,  cx-0.18,  cx-0.28,
        #           5-tip1-V     6-tip2-V     7-Rnode-V    8-Rfork-H
                 cx-0.28,  cx-0.08,  cx+0.18,  cx+0.08,
        #           9-tip3-V     10-tip4-V
                 cx+0.08,  cx+0.28),
        xend = c(cx,       cx+0.18,  cx-0.18,  cx-0.08,
                 cx-0.28,  cx-0.08,  cx+0.18,  cx+0.28,
                 cx+0.08,  cx+0.28),
        y    = c(0.635,    0.695,    0.695,    0.755,
                 0.755,    0.755,    0.695,    0.755,
                 0.755,    0.755),
        yend = c(0.695,    0.695,    0.755,    0.755,
                 0.835,    0.835,    0.755,    0.755,
                 0.835,    0.835)
      )
      tr_tips <- data.frame(
        x = c(cx-0.28, cx-0.08, cx+0.08, cx+0.28),
        y = rep(0.835, 4),
        r = 0.026
      )
      list(
        geom_segment(data = tr_segs,
                     aes(x = x, xend = xend, y = y, yend = yend),
                     colour = pal$node, linewidth = 0.80, alpha = 0.88),
        geom_circle(data = tr_tips,
                    aes(x0 = x, y0 = y, r = r),
                    fill = pal$node, colour = NA, alpha = 0.88)
      )
    } +

    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    )
}

# ── sticker wrapper ───────────────────────────────────────────────────────────
make_sticker <- function(pal, filename) {
  p <- make_subplot(pal)
  hexSticker::sticker(
    subplot  = p,
    package  = "pigauto",
    p_size   = 21,
    p_color  = "white",
    p_family = ttf,
    p_fontface = "bold",
    p_x      = 1.00,
    p_y      = 0.38,
    s_x      = 1.00,
    s_y      = 0.85,
    s_width  = 1.55,
    s_height = 1.55,
    h_fill   = pal$bg,
    h_color  = pal$border,
    h_size   = 1.2,
    filename = filename,
    dpi      = 300
  )
  invisible(filename)
}

# ── generate ──────────────────────────────────────────────────────────────────
dir.create("man/figures",      showWarnings = FALSE, recursive = TRUE)
dir.create("pkgdown/assets",   showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(PALETTES)) {
  pal   <- PALETTES[[i]]
  fname <- sprintf("man/figures/logo_v%d.png", i)
  make_sticker(pal, fname)
  # also copy to pkgdown/assets/ for the gallery page
  file.copy(fname, sprintf("pkgdown/assets/logo_v%d.png", i), overwrite = TRUE)
  message(sprintf("[v%d] %-35s -> %s", i, pal$name, fname))
}

# default logo (v1 unless one already exists and is NOT one of the variants)
if (!file.exists("man/figures/logo.png")) {
  file.copy("man/figures/logo_v1.png", "man/figures/logo.png")
  message("Default logo.png set to v1 (copy to change).")
}

message("\nDone. Open pkgdown/assets/logo_gallery.html in a browser to compare.")
