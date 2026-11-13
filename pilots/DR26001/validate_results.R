#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("SCRNA_R_LIB"), .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript validate_results.R RUN_DIRECTORY")
run_dir <- normalizePath(args[1], mustWork = TRUE)

groups <- c(
  "mouse-ileum-pilot" = "Mouseileum",
  "mouse-duodenum-pilot" = "MouseDuo"
)

checks <- list()
for (group in names(groups)) {
  library_id <- groups[[group]]
  object_path <- file.path(run_dir, "groups", group, "seurat.rds")
  ledger_path <- file.path(run_dir, "libraries", library_id, "cell_ledger.csv")
  stopifnot(file.exists(object_path), file.exists(ledger_path))

  object <- readRDS(object_path)
  ledger <- read.csv(ledger_path, stringsAsFactors = FALSE)
  retained <- ledger[ledger$retained, , drop = FALSE]

  stopifnot(
    ncol(object) == nrow(retained),
    setequal(colnames(object), retained$cell_id),
    all(retained$doublet_class == "singlet"),
    all(retained$ambient_completed),
    all(c("pca", "umap.baseline") %in% SeuratObject::Reductions(object)),
    "baseline_clusters" %in% colnames(object[[]]),
    all(c("feature_id", "symbol", "feature_type", "seurat_id") %in% colnames(object[["RNA"]][[]])),
    isTRUE(object@misc$pipeline$ambient_required),
    isTRUE(object@misc$pipeline$doublet_removal_required)
  )

  corrected <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  observed <- SeuratObject::LayerData(object, assay = "RNA_observed", layer = "counts")
  stopifnot(
    identical(dim(corrected), dim(observed)),
    identical(dimnames(corrected), dimnames(observed)),
    all(is.finite(corrected@x)), all(corrected@x >= 0),
    all(is.finite(observed@x)), all(observed@x >= 0),
    all(abs(observed@x - round(observed@x)) < 1e-8),
    all((observed - corrected)@x >= -1e-7)
  )

  checks[[group]] <- data.frame(
    analysis_group = group,
    input_cells = nrow(ledger),
    retained_cells = nrow(retained),
    predicted_doublets = sum(ledger$doublet_class == "doublet"),
    captures = paste(sort(unique(ledger$capture_id)), collapse = ";"),
    contamination_median = median(ledger$contamination, na.rm = TRUE),
    clusters = length(unique(object$baseline_clusters)),
    validation = "passed"
  )
}

summary <- do.call(rbind, checks)
write.csv(summary, file.path(run_dir, "validation_summary.csv"), row.names = FALSE)
cat("Independent saved-object validation passed\n")
print(summary)
