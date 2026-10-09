# Build a release-style local pkgdown site in a network-restricted audit
# environment. The temporary overrides avoid pkgdown's live CRAN metadata
# lookups and Google Fonts; source configuration remains unchanged.
pkgdown::build_site(
  new_process = FALSE,
  install = TRUE,
  override = list(
    template = list(bslib = list(
      base_font = "system-ui",
      heading_font = "system-ui",
      code_font = "monospace"
    )),
    development = list(mode = "release"),
    home = list(sidebar = FALSE),
    news = list(cran_dates = FALSE)
  )
)
