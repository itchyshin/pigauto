pkgdown::build_site(
  pkg = ".",
  new_process = FALSE,
  install = FALSE,
  quiet = FALSE,
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
