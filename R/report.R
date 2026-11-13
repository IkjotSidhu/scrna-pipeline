html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
html_table <- function(x, limit = 100L) {
  if (!nrow(x)) return("<p>No rows.</p>")
  x <- head(x, limit)
  paste0("<table><thead><tr>", paste0("<th>", html_escape(names(x)), "</th>", collapse = ""),
    "</tr></thead><tbody>", paste(apply(x, 1, function(r)
      paste0("<tr>", paste0("<td>", html_escape(r), "</td>", collapse = ""), "</tr>")), collapse = ""), "</tbody></table>")
}
retention_table <- function(summary) {
  required <- c("library_id", "input_cells", "assessed", "predicted_doublets", "retained")
  for (column in setdiff(required, names(summary))) summary[[column]] <- NA_real_
  out <- summary[, required, drop = FALSE]
  out$assessed_percent <- 100 * out$assessed / out$input_cells
  out$doublet_percent <- 100 * out$predicted_doublets / out$assessed
  out$retained_percent <- 100 * out$retained / out$input_cells
  integer_columns <- c("input_cells", "assessed", "predicted_doublets", "retained")
  out[integer_columns] <- lapply(out[integer_columns], function(x)
    ifelse(is.na(x), "Unavailable", format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)))
  percentage_columns <- c("assessed_percent", "doublet_percent", "retained_percent")
  out[percentage_columns] <- lapply(out[percentage_columns], function(x)
    ifelse(is.na(x), "Unavailable", paste0(format(round(x, 1), nsmall = 1, trim = TRUE), "%")))
  names(out) <- c("Dataset", "Called cells", "Assessed after basic QC", "Predicted doublets",
    "Retained singlets", "Assessed (% called)", "Doublets (% assessed)", "Retained (% called)")
  out
}
qc_policy_table <- function(cfg) {
  if (is.null(cfg$qc)) return(data.frame())
  rules <- data.frame(
    Check = c("Minimum UMIs", "Minimum detected genes", "Mitochondrial RNA", "High detected genes", "Doublets"),
    Setting = c(paste0("< ", cfg$qc$min_counts, " UMIs"),
      paste0("< ", cfg$qc$min_features, " genes"),
      paste0("> ", cfg$qc$max_mt_percent, "% mitochondrial RNA"),
      paste0("> ", cfg$qc$high_features_review, " genes"),
      "scDblFinder prediction"),
    Action = c("Remove", "Remove", "Remove", "Review only", "Remove"),
    stringsAsFactors = FALSE
  )
  rules
}
quality_flag_table <- function(summary) {
  required <- c("library_id", "assessed", "high_features_review", "high_mt")
  if (length(setdiff(required, names(summary)))) return(data.frame())
  out <- summary[, required, drop = FALSE]
  for (column in c("high_features_review", "high_mt")) {
    out[[paste0(column, "_percent")]] <- 100 * out[[column]] / out$assessed
  }
  out <- out[, c("library_id", "high_features_review", "high_features_review_percent", "high_mt", "high_mt_percent")]
  out$high_features_review <- format(round(out$high_features_review), big.mark = ",", scientific = FALSE, trim = TRUE)
  out$high_mt <- format(round(out$high_mt), big.mark = ",", scientific = FALSE, trim = TRUE)
  for (column in c("high_features_review_percent", "high_mt_percent")) {
    out[[column]] <- ifelse(is.na(out[[column]]), "Unavailable",
      paste0(format(round(out[[column]], 1), nsmall = 1, trim = TRUE), "%"))
  }
  names(out) <- c("Dataset", "High-gene review flags", "High-gene flags (% assessed)",
    "High-mito failures", "High-mito failures (% assessed)")
  out
}
qc_comparison_table <- function(qc) {
  required <- c("library_id", "stage", "n_counts", "n_features")
  if (length(setdiff(required, names(qc)))) stop("QC report data are missing required columns", call. = FALSE)
  if (!"dataset_role" %in% names(qc)) qc$dataset_role <- "unspecified"
  if (!"tissue" %in% names(qc)) qc$tissue <- "unspecified"
  selected_stage <- if ("retained" %in% qc$stage) "retained" else "original"
  qc <- qc[qc$stage == selected_stage, , drop = FALSE]
  groups <- split(qc, qc$library_id, drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(x) data.frame(
    dataset = x$library_id[1], role = x$dataset_role[1], tissue = x$tissue[1],
    cells = nrow(x),
    mean_counts = mean(x$n_counts), median_counts = median(x$n_counts),
    mean_features = mean(x$n_features), median_features = median(x$n_features),
    mean_mt_percent = mean(x$percent_mt), median_mt_percent = median(x$percent_mt))))
  out <- out[order(out$dataset), , drop = FALSE]
  numeric <- c("cells", "mean_counts", "median_counts", "mean_features", "median_features")
  out[numeric] <- lapply(out[numeric], function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE))
  out[c("mean_mt_percent", "median_mt_percent")] <- lapply(out[c("mean_mt_percent", "median_mt_percent")],
    function(x) format(round(x, 1), nsmall = 1, trim = TRUE))
  rownames(out) <- NULL
  names(out) <- c("Dataset", "Role", "Tissue", "Final singlets", "Mean counts", "Median counts",
                  "Mean features", "Median features", "Mean mito %", "Median mito %")
  out
}
cellranger_qc_table <- function(metrics) {
  wanted <- c("Cells", "Mean reads per cell", "Median UMI counts per cell",
    "Median genes per cell", "Sequencing saturation", "Confidently mapped reads in cells",
    "Valid barcodes", "Q30 barcodes", "Q30 UMI", "Q30 RNA read")
  if (!nrow(metrics) || !all(c("library_id", "metric_name", "metric_value") %in% names(metrics)))
    return(data.frame())
  x <- metrics[metrics$metric_name %in% wanted, , drop = FALSE]
  if ("library_type" %in% names(x)) x <- x[x$library_type %in% c("", "Gene Expression"), , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  x$metric_rank <- match(x$metric_name, wanted)
  x$category_rank <- ifelse(x$category == "Cells", 0L, 1L)
  x <- x[order(x$library_id, x$metric_rank, x$category_rank), , drop = FALSE]
  x <- x[!duplicated(x[c("library_id", "metric_name")]), c("library_id", "metric_name", "metric_value"), drop = FALSE]
  names(x) <- c("Dataset", "Cell Ranger metric", "Value")
  rownames(x) <- NULL
  x
}
write_reports <- function(run_dir, data_file = "report_data.rds", plot_dir = "plots") {
  need("jsonlite")
  bundle <- readRDS(file.path(run_dir, data_file))
  qc <- bundle$qc
  manifest_path <- file.path(run_dir, "samples.csv")
  legacy_role <- if (grepl("benchmark", bundle$project, ignore.case = TRUE)) "public_benchmark" else "study"
  if ((!"dataset_role" %in% names(qc) || !"tissue" %in% names(qc)) && file.exists(manifest_path)) {
    manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    match_row <- match(qc$library_id, manifest$library_id)
    if (!"dataset_role" %in% names(qc))
      qc$dataset_role <- if ("dataset_role" %in% names(manifest)) manifest$dataset_role[match_row] else legacy_role
    if (!"tissue" %in% names(qc)) qc$tissue <- manifest$tissue[match_row]
  }
  qc$library_id <- factor(qc$library_id)
  comparison <- qc_comparison_table(qc)
  retention <- retention_table(bundle$summary)
  policy <- qc_policy_table(bundle$config)
  flags <- quality_flag_table(bundle$summary)
  utils::write.csv(comparison, file.path(run_dir, "qc_comparison.csv"), row.names = FALSE)
  utils::write.csv(retention, file.path(run_dir, "cell_retention.csv"), row.names = FALSE)
  if (nrow(policy)) utils::write.csv(policy, file.path(run_dir, "qc_policy.csv"), row.names = FALSE)
  if (nrow(flags)) utils::write.csv(flags, file.path(run_dir, "quality_flags.csv"), row.names = FALSE)
  cellranger <- if (!is.null(bundle$cellranger_metrics)) bundle$cellranger_metrics else data.frame()
  cellranger_qc <- cellranger_qc_table(cellranger)
  if (nrow(cellranger_qc)) utils::write.csv(cellranger_qc, file.path(run_dir, "cellranger_qc.csv"), row.names = FALSE)
  status <- jsonlite::read_json(file.path(run_dir, "status.json"), simplifyVector = TRUE)
  header <- function(title) paste0('<!doctype html><html><head><meta charset="utf-8"><title>', title, '</title>',
    '<style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#192b3b}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px;border-bottom:1px solid #ddd;text-align:left}img{max-width:100%}code{background:#eee;padding:2px}h2{margin-top:38px}.note{background:#eff5fa;padding:16px}.warn{background:#fff4dc;padding:16px}</style></head><body>')
  technical <- paste0(header("Single-cell bioinformatics report"),
    '<h1>', html_escape(bundle$project), '</h1><p>Run status: <strong>', html_escape(status$status), '</strong></p>',
    '<div class="note">', if (bundle$mode == "qc") 'Input inspection only. Ambient correction and doublet removal have not run. This is not a final analysis.' else
      'Mandatory ambient correction and doublet removal completed. These are computational estimates, not proof that every contaminant or doublet was removed. Review diagnostics before biological interpretation.', '</div>',
    '<h2>Cell retention through filtering</h2><p>This shows how many Cell Ranger-called cells remain after basic QC, doublet assessment, and final filtering. Download the table as <a href="cell_retention.csv">CSV</a>.</p>',
    html_table(retention),
    '<p>Retained cells passed the configured QC rules and were classified as singlets. Ambient correction and doublet calls are estimates that should be reviewed in the saved cell-level files.</p>',
    if (nrow(policy)) paste0('<h2>QC policy used for this run</h2><p>These are explicit starting rules, saved with the run so they can be reviewed and changed. The high-gene rule is diagnostic only; it does not remove cells. Download the table as <a href="qc_policy.csv">CSV</a>.</p>', html_table(policy)) else '',
    if (nrow(flags)) paste0('<h2>QC flags to review</h2><p>High-gene cells are counted here for review alongside the doublet results. High mitochondrial RNA is an exclusion criterion in this initial policy. Download the table as <a href="quality_flags.csv">CSV</a>.</p>', html_table(flags)) else '',
    '<h2>Final-cell QC comparison</h2><p>All values below are from retained singlets only. Compare each study dataset with a compatible benchmark. Mean values make direct comparisons easy; medians describe the typical cell. Download the same table as <a href="qc_comparison.csv">CSV</a>.</p>',
    html_table(comparison),
    if (nrow(cellranger_qc)) paste0('<h2>10x sequencing and library QC</h2><p>These metrics come directly from the Cell Ranger <code>metrics_summary.csv</code> file. They are available only for captures where that file was supplied. Download the table as <a href="cellranger_qc.csv">CSV</a>.</p>', html_table(cellranger_qc)) else '',
    '<h2>Provenance and review</h2><p>See <a href="resolved_config.json">resolved configuration</a>, ',
    '<a href="samples.csv">sample manifest</a>, <a href="provenance.json">input fingerprints and versions</a>, and <a href="sessionInfo.txt">R session</a>.</p>',
    '<p>Annotation review sheets contain unassigned labels. No Biomni call or automatic cell-type assignment has been made.</p></body></html>')
  biological_summary <- bundle$summary[, intersect(c("library_id", "tissue", "input_cells", "retained"), names(bundle$summary)), drop = FALSE]
  biological <- paste0(header("Single-cell biological review"),
    '<h1>', html_escape(bundle$project), ': biological review</h1>',
    '<div class="warn">No cell-type labels have been reviewed for this run. This report intentionally does not show UMAPs: a cluster map without marker or reference evidence is not an interpretation.</div>',
    '<h2>Samples and retained cells</h2>', html_table(biological_summary),
    '<p>Once marker-based or reference-assisted labels are reviewed, this report will summarize cell types, proportions, and the biological questions supported by those labels. Technical QC is in the <a href="report-bioinformatics.html">bioinformatics report</a>.</p></body></html>')
  index <- paste0(header("Single-cell reports"), '<h1>', html_escape(bundle$project), '</h1><ul>',
    '<li><a href="report-bioinformatics.html">Bioinformatics report</a>: processing summary, QC comparison table, and provenance.</li>',
    '<li><a href="report-biologist.html">Biological review report</a>: retained samples and reviewed biological conclusions when annotations are available.</li>',
    '</ul></body></html>')
  writeLines(technical, file.path(run_dir, "report-bioinformatics.html"))
  writeLines(biological, file.path(run_dir, "report-biologist.html"))
  writeLines(index, file.path(run_dir, "report.html"))
  invisible(file.path(run_dir, c("report.html", "report-bioinformatics.html", "report-biologist.html")))
}
