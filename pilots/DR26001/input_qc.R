# Read-only pilot input validation. This is not a completed preprocessing run.
.libPaths(c(Sys.getenv("SCRNA_R_LIB"), .libPaths()))
for (f in c("config.R", "import.R", "report.R")) source(file.path("R", f))
need(c("Matrix", "hdf5r", "jsonlite", "ggplot2"))
root <- "/sc/arion/projects/naiklab/ikjot/DR_scRNA_001/rawdata/DR26001/cellranger_output/SHNA04_DR26001_0_v1"
out <- file.path("pilots", "DR26001", paste0("input-qc-", Sys.getenv("LSB_JOBID", "local")))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
status <- function(value, message = "") jsonlite::write_json(list(status = value, mode = "qc", message = message,
  ambient_completed = FALSE, doublet_removal_completed = FALSE), file.path(out, "status.json"), pretty = TRUE, auto_unbox = TRUE)
status("running")
tryCatch({
  tags <- jsonlite::read_json(file.path(root, "multi/multiplexing_analysis/cells_per_tag.json"), simplifyVector = TRUE)
  assignment <- do.call(rbind, lapply(names(tags), function(t) data.frame(ocm_group = t, barcode = unlist(tags[[t]]))))
  if (anyDuplicated(assignment$barcode)) stop("Some called barcodes belong to multiple OCM groups")
  expected <- list(Mouseileum = c("OB1", "OB2"), MouseDuo = c("OB3", "OB4"))
  summaries <- list(); metrics <- list(); provenance <- list(); all_called <- character()
  for (sample in names(expected)) {
    message("Checking ", sample)
    sample_root <- file.path(root, "per_sample_outs", sample, "count")
    filtered_path <- file.path(sample_root, "sample_filtered_feature_bc_matrix.h5")
    raw_path <- file.path(sample_root, "sample_raw_feature_bc_matrix.h5")
    filtered <- read_10x(filtered_path)
    q <- qc_metrics(filtered$counts, filtered$features, "mouse")
    q$ocm_group <- assignment$ocm_group[match(q$barcode, assignment$barcode)]
    if (anyNA(q$ocm_group) || any(!q$ocm_group %in% expected[[sample]])) stop("Invalid OCM mapping for ", sample)
    raw <- read_10x(raw_path)
    raw_counts <- align_raw(raw, filtered)
    message("Raw/filtered values agree: ", sample)
    if (any(q$barcode %in% all_called)) stop("Called barcodes overlap between biological sample outputs")
    all_called <- c(all_called, q$barcode)
    q$sample_id <- sample
    q$donor_id <- NA_character_
    q$condition <- NA_character_
    write.csv(q, file.path(out, paste0(sample, "_cell_qc.csv")), row.names = FALSE)
    write.csv(filtered$features, file.path(out, paste0(sample, "_features.csv")), row.names = FALSE)
    totals <- sort(Matrix::colSums(raw_counts), decreasing = TRUE)
    rank <- seq_along(totals)
    selected <- unique(as.integer(round(exp(seq(0, log(length(rank)), length.out = 10000)))))
    rank_data <- data.frame(rank = rank[selected], umi = totals[selected])
    rank_data <- rank_data[rank_data$umi > 0, ]
    p <- ggplot2::ggplot(rank_data, ggplot2::aes(rank, umi)) + ggplot2::geom_line() +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + ggplot2::theme_bw() +
      ggplot2::labs(title = paste(sample, "raw barcode rank — both technical replicates"), x = "Barcode rank", y = "UMIs")
    ggplot2::ggsave(file.path(out, paste0(sample, "_barcode_rank.png")), p, width = 8, height = 5)
    for (tag in expected[[sample]]) {
      v <- q[q$ocm_group == tag, ]
      summaries[[tag]] <- data.frame(library_id = tag, sample_id = sample, species = "mouse",
        input_cells = nrow(v), median_counts = median(v$n_counts), median_features = median(v$n_features),
        median_mt_percent = median(v$percent_mt, na.rm = TRUE), ambient_method_planned = "soupx")
    }
    z <- q[, c("n_counts", "n_features", "percent_mt")]
    z$library_id <- q$ocm_group; z$stage <- "original"
    metrics[[sample]] <- z
    provenance[[sample]] <- list(files = as.list(tools::md5sum(c(filtered_path, raw_path))),
      raw_barcodes = ncol(raw_counts), filtered_cells = ncol(filtered$counts), raw_filtered_equal = TRUE)
    rm(raw, raw_counts, filtered)
    gc(verbose = FALSE)
  }
  if (!setequal(all_called, assignment$barcode)) stop("OCM map and called-cell matrices differ")
  summary <- do.call(rbind, summaries)
  write.csv(summary, file.path(out, "capture_summary.csv"), row.names = FALSE)
  write.csv(assignment, file.path(out, "cell_ocm_map.csv"), row.names = FALSE)
  write.csv(data.frame(sample_id = names(expected), donor_id = NA_character_, condition = NA_character_,
    replicate_type = "technical", metadata_status = "donor_and_condition_unconfirmed"), file.path(out, "samples.csv"), row.names = FALSE)
  jsonlite::write_json(list(mode = "qc_only", filter_applied = FALSE, reference = "mm10-2020-A",
    cellranger = "9.0.1", note = "OCM groups are technical replicates, not biological replicates"),
    file.path(out, "resolved_config.json"), pretty = TRUE, auto_unbox = TRUE)
  jsonlite::write_json(list(inputs = provenance, code = as.list(tools::md5sum(c("R/import.R", "R/config.R", "R/report.R", "pilots/DR26001/input_qc.R"))),
    environment = as.list(tools::md5sum("environment/renv.baseline.lock"))), file.path(out, "provenance.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  saveRDS(list(project = "DR26001 — input QC by OCM group", mode = "qc", summary = summary,
    qc = do.call(rbind, metrics), embeddings = list()), file.path(out, "report_data.rds"))
  status("qc_only", "No filtering, ambient correction or doublet removal has run")
  write_report(out)
  print(summary)
  message("INPUT_QC_OUTPUT=", normalizePath(out))
}, error = function(e) {
  status("failed", conditionMessage(e))
  stop(e)
})
