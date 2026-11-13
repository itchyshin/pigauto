# Source identity for validation evidence, independent of the working directory.
validation_source_identity <- function(root = ".") {
  root <- normalizePath(root, mustWork = TRUE)
  rel <- c("DESCRIPTION", "NAMESPACE", list.files(file.path(root, "R"), pattern = "\\.R$", full.names = FALSE))
  rel[-c(1, 2)] <- file.path("R", rel[-c(1, 2)])
  rel <- sort(unique(c(rel, list.files(file.path(root,"script"), pattern="\\.(R|sh|sbatch|json|ya?ml)$", recursive=TRUE, full.names=TRUE))))
  rel <- sub(paste0(root, "/"), "", rel, fixed=TRUE)
  rel <- sort(unique(c(rel, sub(paste0(root,"/"),"",list.files(file.path(root,"data"),recursive=TRUE,full.names=TRUE),fixed=TRUE))))
  files <- setNames(as.character(tools::md5sum(file.path(root,rel))), rel)
  input_paths <- Sys.getenv(c("VALIDATION_BIEN_TRAITS","VALIDATION_BIEN_TREE"))
  inputs <- setNames(vapply(input_paths,function(p) { if(!nzchar(p))return("unavailable");if(!file.exists(p))stop("Configured input missing: ",p);unname(tools::md5sum(p)) },character(1)),names(input_paths))
  revision <- Sys.getenv("VALIDATION_SOURCE_REVISION", "")
  if(!nzchar(revision)) revision <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout=TRUE, stderr=FALSE)[1]), error=function(e) NA_character_)
  list(hash=digest::digest(list(files=files,inputs=inputs), algo="sha256"), files=files, inputs=inputs, revision=revision,
       package_version=read.dcf(file.path(root,"DESCRIPTION"))[1,"Version"])
}

validation_runtime_identity <- function() {
 installed <- utils::installed.packages()
 installed <- installed[order(installed[,"Package"],installed[,"LibPath"]),intersect(c("Package","Version","LibPath","Built"),colnames(installed)),drop=FALSE]
 environment <- Sys.getenv(c("TORCH_HOME","R_LIBS_USER","OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","PIG_TORCH_THREADS"))
 runtime <- list(RNGkind=RNGkind(),R=R.version,OS=Sys.info()[c("sysname","release","machine")],libraries=.libPaths(),packages=installed,external=extSoftVersion(),lapack=La_version(),blas=tryCatch(La_library(),error=function(e)"unknown"),environment=environment)
 runtime$hash <- digest::digest(runtime,algo="sha256")
 runtime
}
