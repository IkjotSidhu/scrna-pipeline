# Capture installed dependencies without modifying their libraries.
# This baseline is incomplete until DecontX and a clean restore are validated.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript snapshot_baseline.R PROJECT_DIR USER_LIBRARY")
project <- normalizePath(args[1], mustWork = TRUE)
.libPaths(c(normalizePath(args[2], mustWork = TRUE), .libPaths()))
packages <- c("Seurat", "SeuratObject", "Matrix", "yaml", "jsonlite", "digest", "ggplot2",
  "hdf5r", "SoupX", "decontX", "scDblFinder", "SingleCellExperiment", "BiocParallel",
  "future", "renv", "testthat", "BiocManager")
available <- vapply(packages, function(p) length(find.package(p, quiet = TRUE)) > 0, logical(1))
if (any(!available)) message("Baseline excludes missing packages: ", paste(packages[!available], collapse = ", "))
renv::snapshot(project = project, library = .libPaths(), packages = packages[available],
  lockfile = file.path(project, "environment", "renv.baseline.lock"), prompt = FALSE)
