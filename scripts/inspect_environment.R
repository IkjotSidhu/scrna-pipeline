# Read-only package audit. Writes reports only to the supplied output directory.
# Usage: Rscript scripts/inspect_environment.R REPORT_DIR [USER_LIBRARY]
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply an audit output directory")
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
default_paths <- .libPaths()
required <- c("Seurat", "SeuratObject", "Matrix", "yaml", "jsonlite", "digest", "ggplot2",
              "hdf5r", "SoupX", "decontX", "scDblFinder", "SingleCellExperiment",
              "BiocParallel", "future", "renv", "testthat", "BiocManager")
resolve <- function(p) {
  path <- find.package(p, quiet = TRUE)
  if (length(path)) path else NA_character_
}
default_resolution <- setNames(vapply(required, resolve, character(1)), required)
if (length(args) > 1) {
  if (!dir.exists(args[2])) stop("User library does not exist: ", args[2])
  .libPaths(c(args[2], .libPaths()))
}
writeLines(c(R.version.string, paste("Platform:", R.version$platform),
  "Default libraries:", default_paths, "Audited library precedence:", .libPaths(),
  paste("R_HOME:", R.home()), paste("BLAS:", extSoftVersion()["BLAS"]),
  capture.output(Sys.info()[c("sysname", "release", "machine")])), file.path(out, "runtime.txt"))
installed <- installed.packages(lib.loc = .libPaths(), fields = c("Repository", "RemoteType", "RemoteRepo", "RemoteUsername", "RemoteRef", "RemoteSha"))
write.csv(installed, file.path(out, "installed_packages.csv"), row.names = FALSE)
records <- lapply(required, function(p) {
  path <- resolve(p)
  version <- built <- repository <- remote_sha <- NA_character_
  if (!is.na(path)) {
    d <- read.dcf(file.path(path, "DESCRIPTION"))
    val <- function(key) if (key %in% colnames(d)) d[1, key] else NA_character_
    version <- val("Version"); built <- val("Built")
    repository <- val("Repository"); remote_sha <- val("RemoteSha")
  }
  warnings <- character()
  error <- ""
  loaded <- tryCatch(withCallingHandlers({loadNamespace(p); TRUE}, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {error <<- conditionMessage(e); FALSE})
  message(p, " ", version, " load=", loaded, if (nzchar(error)) paste0(" ", error) else "")
  data.frame(package = p, version = version, default_path = default_resolution[p],
    audited_path = path, loaded = loaded, built = built, repository = repository,
    remote_sha = remote_sha, error = error, warnings = paste(unique(warnings), collapse = " | "), row.names = NULL)
})
write.csv(do.call(rbind, records), file.path(out, "required_packages.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Audit saved:", normalizePath(out), "\n")
