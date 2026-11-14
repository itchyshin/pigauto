args <- commandArgs(TRUE)
root <- normalizePath(args[1])
file <- args[2]
meta <- jsonlite::read_json(file.path(root, 'input.json'), simplifyVector = TRUE)
lib <- normalizePath(meta$library)
stopifnot(identical(Sys.getenv('NOT_CRAN', unset = NA_character_), NA_character_))
stopifnot(identical(normalizePath(.libPaths()[1]), lib))
stopifnot(all(normalizePath(.libPaths()) %in% normalizePath(c(lib, .Library))))
library(pigauto)
stopifnot(identical(normalizePath(find.package('pigauto')), normalizePath(file.path(lib, 'pigauto'))))
stopifnot(as.character(packageVersion('pigauto')) == '0.11.0')
stopifnot(as.character(packageVersion('pigauto')) == read.dcf(file.path(lib, 'pigauto', 'DESCRIPTION'))[1, 'Version'])
torch::torch_set_num_threads(1L)
stopifnot(torch::torch_is_installed(), torch::torch_tensor(3)$item() == 3)
print(list(pigauto = find.package('pigauto'), version = packageVersion('pigauto'),
           library_paths = .libPaths(), torch_installed = torch::torch_is_installed(),
           pandoc = rmarkdown::find_pandoc()))
print(sessionInfo())
output <- file.path(root, 'rendered', sub('.Rmd$', '', file))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
# This observer records authored eval settings and returns each option list unchanged.
knitr::opts_hooks$set(prooftrace = function(options) {
  cat(jsonlite::toJSON(list(label = options$label, eval = options$eval), auto_unbox = TRUE), '\n',
      file = file.path(output, 'chunk-options.jsonl'), append = TRUE)
  options
})
knitr::opts_chunk$set(prooftrace = TRUE)
source <- file.path(root, 'source', 'pigauto', 'vignettes', file)
started <- Sys.time()
startcpu <- proc.time()
result <- tryCatch(rmarkdown::render(source, output_dir = output, intermediates_dir = output,
                                     envir = new.env(parent = globalenv()), clean = FALSE, quiet = FALSE),
                   error = identity)
elapsed <- as.numeric(difftime(Sys.time(), started, units = 'secs'))
cpu <- proc.time() - startcpu
rec <- list(file = file, start_utc = format(started, tz = 'UTC', usetz = TRUE),
            elapsed_seconds = elapsed, user_seconds = unname(cpu[['user.self']]),
            system_seconds = unname(cpu[['sys.self']]), output = if (is.character(result)) result else NULL,
            error = if (inherits(result, 'error')) conditionMessage(result) else NULL)
jsonlite::write_json(rec, file.path(output, 'render-result.json'), auto_unbox = TRUE, pretty = TRUE)
if (inherits(result, 'error')) stop(conditionMessage(result))
stopifnot(file.exists(result))
cat('VIGNETTE_RENDER_COMPLETE ', file, ' elapsed=', elapsed, '\n', sep = '')
