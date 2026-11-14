# The Windows failures named these missing Courier metrics. Fail before check
# unless both metrics resolve and a small T1/TS1 Courier document compiles.
tinytex::tlmgr_install("courier")
metrics <- c("pcrr8t.tfm", "pcrr8c.tfm")
fonts <- lapply(metrics, function(metric) {
  path <- system2("kpsewhich", metric, stdout = TRUE)
  stopifnot(length(path) == 1L, nzchar(path))
  path
})
names(fonts) <- metrics
stopifnot(V8::v8()$eval("1 + 1") == "2")
capture.output(list(courier_metrics = fonts, V8 = V8::engine_info()),
               file = "platform-evidence/manual-tools.txt")
# Probe index processing before a full check. Install makeindex only if the
# executable is observably absent, and keep the failed probe as evidence.
capture.output(list(paths = Sys.which(c("pdflatex", "makeindex", "texi2dvi", "texify")),
                    environment = Sys.getenv(c("MAKEINDEX", "R_MAKEINDEXCMD", "R_TEXI2DVICMD"))),
               file = "platform-evidence/index-tools-before.txt")
original <- getwd()
setwd("platform-evidence")
writeLines("\\indexentry{pigauto}{1}", "index-probe.idx")
probe <- function(label) {
  tryCatch(system2("makeindex", c("-o", "index-probe.ind", "index-probe.idx"),
                   stdout = paste0(label, "-stdout.txt"),
                   stderr = paste0(label, "-stderr.txt")),
           error = function(e) {
             writeLines(conditionMessage(e), paste0(label, "-error.txt"))
             -1L
           })
}
before <- probe("index-before")
writeLines(paste0("exit=", before), "index-before-exit.txt")
if (!nzchar(Sys.which("makeindex"))) {
  tinytex::tlmgr_install("makeindex")
}
after <- probe("index-after")
writeLines(c(paste0("exit=", after), paste0("makeindex=", Sys.which("makeindex"))),
           "index-after-exit.txt")
stopifnot(after == 0L, file.exists("index-probe.ind"),
          file.info("index-probe.ind")$size > 0)
writeLines(c("\\documentclass{article}", "\\usepackage[T1]{fontenc}",
             "\\usepackage{textcomp,courier,url,makeidx}", "\\makeindex",
             "\\begin{document}",
             "{\\fontfamily{pcr}\\selectfont Courier \\textasciigrave}",
             "\\url{https://itchyshin.github.io/pigauto/}\\index{pigauto}",
             "\\printindex", "\\end{document}"), "courier-canary.tex")
# Exercise R's own indexed-PDF orchestration, not only pdflatex directly.
tools::texi2pdf("courier-canary.tex", clean = FALSE, quiet = FALSE, index = TRUE)
setwd(original)
stopifnot(file.exists("platform-evidence/courier-canary.pdf"),
          file.exists("platform-evidence/courier-canary.ind"),
          file.info("platform-evidence/courier-canary.ind")$size > 0)
