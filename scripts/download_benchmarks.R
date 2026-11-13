#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/download_benchmarks.R PROJECT_DIR")
project <- normalizePath(args[1], mustWork = TRUE)
manifest <- read.csv(file.path(project, "benchmarks", "files.csv"), stringsAsFactors = FALSE)
root <- file.path(project, "benchmarks", "data")
dir.create(root, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  directory <- file.path(root, row$dataset)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(directory, row$file)
  if (!file.exists(destination) || unname(tools::md5sum(destination)) != row$md5) {
    message("Downloading ", row$dataset, "/", row$file)
    utils::download.file(row$url, destination, mode = "wb", quiet = FALSE)
  }
  observed <- unname(tools::md5sum(destination))
  if (!identical(observed, row$md5))
    stop("Checksum mismatch for ", destination, ": expected ", row$md5, ", got ", observed)
}

skin_root <- file.path(root, "GSE142471")
archive <- file.path(skin_root, "GSE142471_RAW.tar")
raw_dir <- file.path(skin_root, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
if (!length(list.files(raw_dir))) utils::untar(archive, exdir = raw_dir)

controls <- c(
  "GSM4230076-Un-Wounded-1" = "GSM4230076_Un-Wounded_1_scRNA-Seq",
  "GSM4230077-Un-Wounded-2" = "GSM4230077_Un-Wounded_2_scRNA-Seq"
)
for (sample_name in names(controls)) {
  control <- controls[[sample_name]]
  sample_dir <- file.path(skin_root, sample_name)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  accession <- sub("_.*", "", control)
  links <- c(
    "matrix.mtx.gz" = file.path(raw_dir, paste0(control, ".mtx.gz")),
    "barcodes.tsv.gz" = file.path(raw_dir, paste0(accession, "_barcodes_", sub("^[^_]+_", "", control), ".tsv.gz")),
    "genes.tsv.gz" = file.path(raw_dir, paste0(accession, "_genes_", sub("^[^_]+_", "", control), ".tsv.gz"))
  )
  for (name in names(links)) {
    target <- links[[name]]
    destination <- file.path(sample_dir, name)
    if (!file.exists(target)) stop("Missing extracted control file: ", target)
    if (!file.exists(destination) && !file.symlink(target, destination))
      stop("Could not link ", destination)
  }
}
writeLines(capture.output(sessionInfo()), file.path(root, "download-sessionInfo.txt"))
cat("Benchmark downloads and checksums completed\n")
