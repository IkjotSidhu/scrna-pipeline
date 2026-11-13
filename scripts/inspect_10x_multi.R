# Inspect metadata and H5 headers without loading count matrices.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript inspect_10x_multi.R CELLRANGER_DIR REPORT_DIR")
.libPaths(c(Sys.getenv("SCRNA_R_LIB"), .libPaths()))
root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
paths <- list.files(root, pattern = "feature_bc_matrix.h5$", full.names = TRUE, recursive = TRUE)
headers <- lapply(paths, function(p) {
  h <- hdf5r::H5File$new(p, "r")
  on.exit(h$close_all())
  shape <- h[["matrix/shape"]][]
  symbols <- h[["matrix/features/name"]][]
  genome <- if (h$exists("matrix/features/genome")) paste(unique(h[["matrix/features/genome"]][]), collapse = ";") else "unknown"
  data.frame(path = p, n_features = shape[1], n_barcodes = shape[2],
    mt_features = sum(grepl("^mt-", symbols)), genome = genome,
    nnz = h[["matrix/data"]]$dims[1])
})
write.csv(do.call(rbind, headers), file.path(out, "matrix_headers.csv"), row.names = FALSE)
tag_path <- file.path(root, "multi/multiplexing_analysis/cells_per_tag.json")
if (file.exists(tag_path)) {
  tags <- jsonlite::read_json(tag_path, simplifyVector = TRUE)
  print(vapply(tags, length, integer(1)))
  map <- do.call(rbind, lapply(names(tags), function(tag) data.frame(ocm_group = tag, barcode = unlist(tags[[tag]]))))
  write.csv(map, file.path(out, "cell_ocm_map.csv"), row.names = FALSE)
}
for (p in c("config.csv", "_versions", "metadata/meta_sample.csv", "metadata/libraries.csv"))
  if (file.exists(file.path(root, p))) file.copy(file.path(root, p), file.path(out, basename(p)), overwrite = FALSE)
print(do.call(rbind, headers))
