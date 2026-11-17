# Torch factories hide real method closures from covr's namespace scan.
# Register only pigauto's four model methods during coverage runs.
if (identical(tolower(Sys.getenv("R_COVR")), "true")) {
  local({
    covr_ns <- asNamespace("covr")
    required <- c("replacement", "replace", "the")
    stopifnot(all(vapply(required, exists, logical(1), envir = covr_ns,
                         inherits = FALSE)))
    ns <- asNamespace("pigauto")
    methods <- list()
    replacements <- list()
    for (factory in c("GraphTransformerBlock", "ResidualPhyloDAE")) {
      generator <- attr(get(paste0(factory, "_module"), envir = ns)(), "module")
      stopifnot(inherits(generator, "R6ClassGenerator"))
      for (method in c("initialize", "forward")) {
        fun <- generator$public_methods[[method]]
        stopifnot(is.function(fun), inherits(attr(fun, "srcref"), "srcref"))
        methods[[length(methods) + 1L]] <- fun
        # Qualified names also keep Cobertura method summaries file-specific.
        replacement <- get("replacement", envir = covr_ns)(
          paste(factory, method, sep = "$"), env = ns, target_value = fun
        )
        if (is.null(replacement)) {
          # covr >= 4 already traces methods reached through an R6 generator.
          # It returns NULL for traced closures to avoid instrumenting them twice.
          is_traced <- get0("is_traced", envir = covr_ns, inherits = FALSE)
          stopifnot(is.function(is_traced), is_traced(fun))
        } else {
          replacements[[length(replacements) + 1L]] <- replacement
        }
      }
    }
    stopifnot(length(methods) == 4L, all(vapply(replacements,
      function(x) is.function(x$target_value) && is.function(x$new_value) &&
        is.function(x$orig_value), logical(1))))
    if (length(replacements)) {
      state <- get("the", envir = covr_ns)
      # Preserve ordinary counters and let covr restore all original closures.
      state$replacements <- c(state$replacements, replacements)
      invisible(lapply(replacements, get("replace", envir = covr_ns)))
    }
  })
}
