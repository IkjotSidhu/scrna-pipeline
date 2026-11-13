# Explicit environment setup; never called by run.R. Uses a private library.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/setup.R /absolute/path/to/project-library")
lib <- path.expand(args[1])
if (!grepl("^/", lib)) stop("Provide an absolute library path")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages(c("BiocManager", "Seurat", "Matrix", "yaml", "jsonlite", "digest", "ggplot2",
                   "hdf5r", "SoupX", "testthat", "renv"), lib = lib)
BiocManager::install(c("decontX", "scDblFinder", "SingleCellExperiment", "BiocParallel"),
                     lib = lib, ask = FALSE, update = FALSE)
cat("Use SCRNA_R_LIB=", lib, " when launching run.R\n", sep = "")
cat("Generate a lockfile only after validating this environment on the intended machine.\n")
