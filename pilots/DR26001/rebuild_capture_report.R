#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("SCRNA_R_LIB"), .libPaths()))
source(file.path("R", "config.R"))
source(file.path("R", "report.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript rebuild_capture_report.R RUN_DIRECTORY")
run_dir <- normalizePath(args[1], mustWork = TRUE)

ledgers <- lapply(c("Mouseileum", "MouseDuo"), function(id) {
  read.csv(file.path(run_dir, "libraries", id, "cell_ledger.csv"), stringsAsFactors = FALSE)
})
ledger <- do.call(rbind, ledgers)

qc <- list()
for (stage in c("original", "corrected", "retained")) {
  keep <- switch(stage,
    original = rep(TRUE, nrow(ledger)),
    corrected = ledger$preliminary_pass,
    retained = ledger$retained
  )
  prefix <- if (stage == "original") "" else "corrected_"
  block <- data.frame(
    n_counts = ledger[[paste0(prefix, "n_counts")]][keep],
    n_features = ledger[[paste0(prefix, "n_features")]][keep],
    percent_mt = ledger[[paste0(prefix, "percent_mt")]][keep],
    library_id = ledger$capture_id[keep],
    stage = stage
  )
  qc[[stage]] <- block
}

old <- readRDS(file.path(run_dir, "report_data.rds"))
embeddings <- lapply(old$embeddings, function(d) {
  d$capture_id <- ledger$capture_id[match(rownames(d), ledger$cell_id)]
  if (anyNA(d$capture_id)) stop("Could not map all embedding cells to capture IDs")
  d
})

summary <- do.call(rbind, lapply(sort(unique(ledger$capture_id)), function(capture) {
  x <- ledger[ledger$capture_id == capture, , drop = FALSE]
  data.frame(
    capture_id = capture,
    input_cells = nrow(x),
    assessed = sum(x$preliminary_pass),
    predicted_doublets = sum(x$doublet_class == "doublet"),
    retained = sum(x$retained),
    retained_percent = round(mean(x$retained) * 100, 1),
    ambient_method = unique(x$ambient_method),
    contamination_median = median(x$contamination, na.rm = TRUE)
  )
}))

bundle <- list(
  project = paste0(old$project, " — capture-aware report"),
  mode = old$mode,
  summary = summary,
  qc = do.call(rbind, qc),
  embeddings = embeddings
)
saveRDS(bundle, file.path(run_dir, "report_data_capture.rds"))
write_report(run_dir, data_file = "report_data_capture.rds",
             output_file = "report-capture.html", plot_dir = "plots-capture")
cat("Capture-aware report written\n")
