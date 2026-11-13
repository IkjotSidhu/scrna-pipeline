#!/usr/bin/env Rscript
required <- c(
  "BiocParallel", "decontX", "digest", "future", "ggplot2", "harmony", "hdf5r",
  "jsonlite", "Matrix", "renv", "scDblFinder", "Seurat", "SeuratObject",
  "SingleCellExperiment", "SoupX", "testthat", "yaml"
)

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
print(audit[, c("package", "version", "loaded")], row.names = FALSE)
cat("\nLibrary paths:\n", paste(.libPaths(), collapse = "\n"), "\n", sep = "")
cat("Platform prefix: ", Sys.getenv("RENV_PATHS_PREFIX", unset = "<unset>"), "\n", sep = "")
if (!all(audit$loaded)) {
  stop("Environment validation failed: ", paste(audit$package[!audit$loaded], collapse = ", "))
}
