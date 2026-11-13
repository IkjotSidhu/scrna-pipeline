#!/usr/bin/env Rscript
# Build the project library from the working Minerva libraries, add missing
# packages, lock it, and restore the lock into a second isolated test library.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/create_environment.R PROJECT_DIR")

project <- normalizePath(args[1], mustWork = TRUE)
external_libraries <- .libPaths()
required <- c(
  "BiocParallel", "decontX", "digest", "future", "ggplot2", "harmony", "hdf5r",
  "jsonlite", "Matrix", "renv", "scDblFinder", "Seurat", "SeuratObject",
  "SingleCellExperiment", "SoupX", "testthat", "yaml"
)

renv::consent(provided = TRUE)
renv::init(
  project = project,
  bare = TRUE,
  bioconductor = "3.20",
  restart = FALSE
)

available <- required[vapply(required, function(p) {
  length(find.package(p, lib.loc = external_libraries, quiet = TRUE)) > 0
}, logical(1))]
message("Hydrating installed dependencies: ", paste(available, collapse = ", "))
renv::hydrate(
  packages = available,
  project = project,
  sources = external_libraries,
  prompt = FALSE
)

if (!requireNamespace("decontX", quietly = TRUE)) {
  message("Installing missing Bioconductor package decontX")
  renv::install("bioc::decontX", project = project, prompt = FALSE)
}

renv::snapshot(
  project = project,
  type = "all",
  packages = required,
  force = TRUE,
  prompt = FALSE
)

lockfile <- file.path(project, "renv.lock")
if (!file.exists(lockfile)) stop("renv.lock was not created")

# Populate the library that normal project runs activate. Hydration and a
# successful lockfile snapshot do not themselves prove this runtime is complete.
project_library <- renv::paths$library(project = project)
renv::restore(
  project = project,
  library = project_library,
  lockfile = lockfile,
  clean = TRUE,
  prompt = FALSE
)

restore_project <- file.path(project, "environment", "restore-test")
dir.create(restore_project, recursive = TRUE, showWarnings = FALSE)
restore_library <- file.path(restore_project, "library")
renv::restore(
  project = restore_project,
  library = restore_library,
  lockfile = lockfile,
  clean = TRUE,
  prompt = FALSE
)

.libPaths(c(restore_library, .Library))
records <- lapply(required, function(package) {
  error <- ""
  loaded <- tryCatch({
    loadNamespace(package)
    TRUE
  }, error = function(e) {
    error <<- conditionMessage(e)
    FALSE
  })
  data.frame(
    package = package,
    version = if (loaded) as.character(packageVersion(package)) else NA_character_,
    loaded = loaded,
    path = if (loaded) find.package(package) else NA_character_,
    error = error
  )
})
audit <- do.call(rbind, records)
write.csv(audit, file.path(project, "environment", "restore-test.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(project, "environment", "restore-test-sessionInfo.txt"))
if (!all(audit$loaded)) stop("Clean restore load test failed: ", paste(audit$package[!audit$loaded], collapse = ", "))
cat("Project environment and clean restore validated\n")
