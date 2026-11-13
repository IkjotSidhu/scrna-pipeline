write_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) abort("Could not save ", path)
}
cached <- function(key, cache_dir, compute) {
  path <- file.path(cache_dir, paste0(key, ".rds"))
  if (file.exists(path)) {
    message("Reusing validated cache ", substr(key, 1, 10))
    return(readRDS(path))
  }
  value <- compute()
  atomic_rds(value, path)
  value
}
run_pipeline <- function(config_path, mode, pipeline_root) {
  choice(mode, "mode", c("qc", "run"))
  cfg <- read_config(config_path)
  samples <- read_manifest(cfg)
  # Required dependencies are checked before expensive import/analysis.
  if (mode == "run") {
    need(c("scDblFinder", "SingleCellExperiment", "BiocParallel"))
    if (any(samples$ambient_method == "soupx")) need("SoupX")
    if (any(samples$ambient_method == "decontx")) need("decontX")
    if (cfg$integration$method == "harmony" && cfg$integration$enabled != "never") need("harmony")
  }
  future::plan(future::multisession, workers = cfg$workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  options(future.globals.maxSize = cfg$future_max_gb * 1024^3)
  set.seed(cfg$seed)
  input_hashes <- lapply(seq_len(nrow(samples)), function(i) {
    row <- samples[i, ]
    files <- input_files(row$filtered_matrix)
    if (!is.na(row$raw_matrix)) files <- c(files, input_files(row$raw_matrix))
    if (!is.na(row$cell_metadata)) files <- c(files, cell_metadata = row$cell_metadata)
    if (!is.na(row$cellranger_metrics)) files <- c(files, cellranger_metrics = row$cellranger_metrics)
    as.list(tools::md5sum(files))
  })
  names(input_hashes) <- samples$library_id
  code_files <- file.path(pipeline_root, "R", c("config.R", "import.R", "process.R", "report.R", "engine.R"))
  code_hash <- tools::md5sum(code_files)
  packages <- c("Seurat", "SeuratObject", "Matrix", "SoupX", "decontX", "scDblFinder", "SingleCellExperiment", "BiocParallel", "digest", "yaml", "jsonlite")
  versions <- setNames(lapply(packages, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not installed"), packages)
  environment <- list(R = R.version.string, packages = versions)
  key <- function(...) digest::digest(list(..., code_hash, environment), algo = "sha256")
  run_key <- key(cfg, samples, input_hashes, mode)
  run_dir <- file.path(cfg$output, paste0(mode, "-", substr(run_key, 1, 12)))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(cfg$output, "cache")
  status_path <- file.path(run_dir, "status.json")
  success <- FALSE
  on.exit(if (!success) write_json(list(status = "failed", mode = mode,
    message = "See execution log. Partial files are not final results."), status_path), add = TRUE)
  write_json(list(status = "running", mode = mode), status_path)
  write_json(cfg, file.path(run_dir, "resolved_config.json"))
  utils::write.csv(samples, file.path(run_dir, "samples.csv"), row.names = FALSE)
  write_json(list(run_key = run_key, inputs = input_hashes, code = as.list(code_hash), environment = environment,
    ambient_required = TRUE, doublet_removal_required = TRUE), file.path(run_dir, "provenance.json"))
  writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
  summaries <- list(); report_qc <- list(); report_cellranger <- list(); processed <- list(); process_keys <- list()
  for (i in seq_len(nrow(samples))) {
    row <- samples[i, , drop = FALSE]
    id <- row$library_id
    message("Processing capture ", id, " (", i, "/", nrow(samples), ")")
    library_dir <- file.path(run_dir, "libraries", id)
    dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
    import_key <- key("import", row, input_hashes[[id]])
    dat <- cached(import_key, cache_dir, function() import_library(row))
    utils::write.csv(dat$features, file.path(library_dir, "features.csv"), row.names = FALSE)
    utils::write.csv(dat$qc, file.path(library_dir, "qc_original.csv"), row.names = FALSE)
    write_json(dat$excluded_modalities, file.path(library_dir, "excluded_modalities.json"))
    cellranger <- read_cellranger_metrics(row$cellranger_metrics)
    if (nrow(cellranger)) {
      cellranger$library_id <- id
      utils::write.csv(cellranger, file.path(library_dir, "cellranger_metrics.csv"), row.names = FALSE)
      report_cellranger[[id]] <- cellranger
    }
    q <- dat$qc[, c("n_counts", "n_features", "percent_mt")]
    q$library_id <- id; q$capture_id <- dat$cell_metadata$capture_id; q$stage <- "original"
    q$sample_id <- row$sample_id; q$batch_id <- row$batch_id
    q$dataset_role <- row$dataset_role; q$tissue <- row$tissue
    report_qc[[paste0(id, "_original")]] <- q
    summary <- data.frame(library_id = id, dataset_role = row$dataset_role,
      species = row$species, tissue = row$tissue,
      capture_groups = paste(sort(unique(dat$cell_metadata$capture_id)), collapse = ";"),
      input_cells = ncol(dat$counts), assessed = NA_integer_, predicted_doublets = NA_integer_,
      retained = NA_integer_, high_features_review = NA_integer_, high_mt = NA_integer_,
      ambient_method = row$ambient_method, median_contamination = NA_real_,
      cellranger_metrics_available = !is.na(row$cellranger_metrics))
    if (mode == "run") {
      process_key <- key("process", import_key, cfg[c("seed", "qc", "doublets", "analysis")])
      p <- cached(process_key, cache_dir, function() process_library(dat, cfg))
      processed[[id]] <- p
      process_keys[[id]] <- process_key
      utils::write.csv(p$ledger, file.path(library_dir, "cell_ledger.csv"), row.names = FALSE)
      atomic_rds(p$corrected_counts, file.path(library_dir, "corrected_counts_assessed.rds"))
      atomic_rds(p$object, file.path(library_dir, "singlets.rds"))
      write_json(p$status, file.path(library_dir, "processing_status.json"))
      summary$assessed <- sum(p$ledger$preliminary_pass)
      summary$predicted_doublets <- sum(p$ledger$doublet_class == "doublet")
      summary$retained <- sum(p$ledger$retained)
      summary$high_features_review <- sum(p$ledger$high_features_review, na.rm = TRUE)
      summary$high_mt <- sum(p$ledger$high_mt, na.rm = TRUE)
      summary$median_contamination <- median(p$ledger$contamination, na.rm = TRUE)
      for (stage in c("corrected", "retained")) {
        select <- if (stage == "corrected") p$ledger$preliminary_pass else p$ledger$retained
        q <- p$ledger[select, c("corrected_n_counts", "corrected_n_features", "corrected_percent_mt")]
        names(q) <- c("n_counts", "n_features", "percent_mt")
        q$library_id <- id; q$capture_id <- p$ledger$capture_id[select]; q$stage <- stage
        q$sample_id <- row$sample_id; q$batch_id <- row$batch_id
        q$dataset_role <- row$dataset_role; q$tissue <- row$tissue
        report_qc[[paste0(id, "_", stage)]] <- q
      }
    }
    summaries[[id]] <- summary
    rm(dat)
    gc(verbose = FALSE)
  }
  embeddings <- list()
  if (mode == "run") for (group in unique(samples$analysis_group)) {
    ids <- samples$library_id[samples$analysis_group == group]
    message("Analyzing group ", group)
    group_key <- key("analysis", process_keys[ids], cfg[c("seed", "analysis", "integration", "markers")])
    result <- cached(group_key, cache_dir, function() analyze_group(processed[ids], cfg))
    group_dir <- file.path(run_dir, "groups", group)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    result$object@misc$run_key <- run_key
    atomic_rds(result$object, file.path(group_dir, "seurat.rds"))
    utils::write.csv(result$markers, file.path(group_dir, "markers.csv"), row.names = FALSE)
    utils::write.csv(result$annotation_review, file.path(group_dir, "annotation_review.csv"), row.names = FALSE)
    prefixes <- "baseline"
    if (cfg$integration$enabled) prefixes <- c(prefixes, cfg$integration$method)
    for (prefix in prefixes) {
      coordinates <- SeuratObject::Embeddings(result$object[[paste0("umap.", prefix)]])
      colnames(coordinates) <- c("UMAP_1", "UMAP_2")
      metadata <- result$object[[]][rownames(coordinates), ]
      embeddings[[paste(group, prefix, sep = "_")]] <- data.frame(coordinates,
        library_id = metadata$library_id, capture_id = metadata$capture_id,
        sample_id = metadata$sample_id, tissue = metadata$tissue,
        condition = metadata$condition,
        cluster = metadata[[paste0(prefix, "_clusters")]])
    }
  }
  summary <- do.call(rbind, summaries)
  utils::write.csv(summary, file.path(run_dir, "capture_summary.csv"), row.names = FALSE)
  atomic_rds(list(project = cfg$project, mode = mode, config = cfg, summary = summary,
    qc = do.call(rbind, report_qc),
    cellranger_metrics = if (length(report_cellranger)) do.call(rbind, report_cellranger) else data.frame(),
    embeddings = embeddings), file.path(run_dir, "report_data.rds"))
  write_json(list(status = if (mode == "qc") "qc_only" else "completed", mode = mode,
    ambient_completed = mode == "run", doublet_removal_completed = mode == "run"), status_path)
  write_reports(run_dir)
  success <- TRUE
  message("Results: ", run_dir)
  invisible(run_dir)
}
