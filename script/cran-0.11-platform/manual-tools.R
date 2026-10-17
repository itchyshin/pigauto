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
original <- getwd()
setwd("platform-evidence")
writeLines(c("\\documentclass{article}", "\\usepackage[T1]{fontenc}",
             "\\usepackage{textcomp,courier,url}", "\\begin{document}",
             "{\\fontfamily{pcr}\\selectfont Courier \\textasciigrave}",
             "\\url{https://itchyshin.github.io/pigauto/}", "\\end{document}"),
           "courier-canary.tex")
status <- system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error",
                               "courier-canary.tex"),
                  stdout = "courier-canary-console.txt", stderr = "courier-canary-console.txt")
setwd(original)
stopifnot(status == 0L, file.exists("platform-evidence/courier-canary.pdf"))
