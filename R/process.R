make_object <- function(counts, features, row) {
  # Seurat replaces underscores/pipes; precompute a stable reversible mapping.
  features$seurat_id <- make.unique(gsub("[_|]", "-", features$feature_id))
  original <- colnames(counts)
  rownames(counts) <- features$seurat_id
  colnames(counts) <- paste(row$library_id, original, sep = "_")
  obj <- Seurat::CreateSeuratObject(counts, project = row$library_id, min.cells = 0, min.features = 0)
  obj$barcode_original <- original
  for (key in c("library_id", "sample_id", "donor_id", "species", "tissue", "condition", "batch_id", "reference_id", "analysis_group", "dataset_role")) obj[[key]] <- row[[key]]
  rownames(features) <- features$seurat_id
  obj[["RNA"]] <- SeuratObject::AddMetaData(obj[["RNA"]], features)
  obj
}
pca_baseline <- function(obj, settings, seed) {
  SeuratObject::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = min(settings$nfeatures, nrow(obj)), verbose = FALSE)
  variable <- SeuratObject::VariableFeatures(obj)
  dims <- min(settings$npcs, ncol(obj) - 1L, length(variable) - 1L)
  if (dims < 2) abort("Insufficient cells/variable genes for PCA")
  obj <- Seurat::ScaleData(obj, features = variable, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, features = variable, npcs = dims, seed.use = seed, verbose = FALSE)
  obj
}
cluster_reduction <- function(obj, reduction, prefix, settings, seed, umap = TRUE) {
  dims <- seq_len(ncol(SeuratObject::Embeddings(obj[[reduction]])))
  graph <- paste0(prefix, "_snn")
  obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = dims,
    k.param = min(settings$neighbors, ncol(obj) - 1L),
    graph.name = c(paste0(prefix, "_nn"), graph), verbose = FALSE)
  obj <- Seurat::FindClusters(obj, graph.name = graph, resolution = settings$resolution,
    cluster.name = paste0(prefix, "_clusters"), random.seed = seed, verbose = FALSE)
  if (umap) obj <- Seurat::RunUMAP(obj, reduction = reduction, dims = dims,
    reduction.name = paste0("umap.", prefix), n.neighbors = min(settings$umap_neighbors, ncol(obj) - 1L),
    umap.method = "uwot", metric = "cosine", seed.use = seed, verbose = FALSE)
  obj
}
process_library <- function(dat, cfg) {
  set.seed(cfg$seed)
  observed <- dat$counts
  row <- dat$row
  q <- dat$qc
  preliminary <- q$n_features >= cfg$qc$preliminary_min_features & q$n_counts > 0
  if (sum(preliminary) < 50) abort(row$library_id, ": fewer than 50 cells after minimal cleanup; review input/QC before mandatory correction")
  x <- observed[, preliminary, drop = FALSE]
  cell_metadata <- dat$cell_metadata[match(colnames(x), dat$cell_metadata$barcode), , drop = FALSE]
  if (anyNA(cell_metadata$capture_id)) abort("Missing capture_id after barcode matching")
  raw <- if (!is.na(row$raw_matrix)) align_raw(read_10x(row$raw_matrix), dat) else NULL
  if (row$ambient_method == "soupx") {
    need("SoupX")
    provisional <- make_object(x, dat$features, row)
    provisional <- pca_baseline(provisional, cfg$analysis, cfg$seed)
    provisional <- cluster_reduction(provisional, "pca", "provisional", cfg$analysis, cfg$seed, umap = FALSE)
    labels <- setNames(as.character(provisional$provisional_clusters), colnames(x))
    if (length(unique(labels)) < 2) abort("SoupX needs informative clusters; only one provisional cluster in ", row$library_id)
    channel <- SoupX::SoupChannel(tod = raw, toc = x)
    channel <- SoupX::setClusters(channel, labels)
    channel <- SoupX::autoEstCont(channel, doPlot = FALSE)
    corrected <- SoupX::adjustCounts(channel, roundToInt = FALSE)
    contamination <- channel$metaData[colnames(x), "rho"]
  } else if (row$ambient_method == "decontx") {
    need(c("decontX", "SingleCellExperiment"))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = x))
    # Exclude all called cells from the background, including preliminary rejects.
    background <- if (!is.null(raw)) raw[, setdiff(colnames(raw), colnames(observed)), drop = FALSE] else NULL
    if (!is.null(background) && !ncol(background)) abort("No background droplets remain for DecontX")
    sce <- decontX::decontX(sce, background = background, seed = cfg$seed, verbose = FALSE)
    corrected <- decontX::decontXcounts(sce)
    contamination <- sce$decontX_contamination
  } else abort("Mandatory ambient correction has no valid method")
  corrected <- methods::as(corrected, "dgCMatrix")
  if (!identical(dimnames(corrected), dimnames(x)) || any(!is.finite(corrected@x)) || any(corrected@x < 0)) abort("Invalid corrected count matrix")
  if (length(contamination) != ncol(x) || any(!is.finite(contamination)) || any(contamination < 0 | contamination > 1)) abort("Invalid contamination estimates")
  # Use observed integer counts to score doublets; ambient-adjusted values remain fractional.
  need(c("scDblFinder", "SingleCellExperiment", "BiocParallel"))
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = x))
  sce$capture_id <- cell_metadata$capture_id
  rate <- if (cfg$doublets$expected_rate == 0) NULL else cfg$doublets$expected_rate
  sce <- scDblFinder::scDblFinder(sce, samples = "capture_id", dbr = rate, dbr.per1k = 0.008,
    BPPARAM = BiocParallel::SerialParam(RNGseed = cfg$seed), verbose = FALSE)
  score <- sce$scDblFinder.score
  label <- as.character(sce$scDblFinder.class)
  if (length(label) != ncol(x) || anyNA(label) || any(!label %in% c("singlet", "doublet")) || any(!is.finite(score))) abort("Invalid doublet results")
  corrected_qc <- qc_metrics(corrected, dat$features, row$species)
  ledger <- q
  ledger$cell_id <- paste(row$library_id, ledger$barcode, sep = "_")
  ledger$capture_id <- dat$cell_metadata$capture_id[match(ledger$barcode, dat$cell_metadata$barcode)]
  ledger$preliminary_pass <- preliminary
  ledger$ambient_method <- row$ambient_method
  ledger$ambient_completed <- preliminary
  ledger$contamination <- NA_real_
  ledger$doublet_score <- NA_real_
  ledger$doublet_class <- "not_tested"
  ledger$contamination[preliminary] <- contamination
  ledger$doublet_score[preliminary] <- score
  ledger$doublet_class[preliminary] <- label
  for (metric in c("n_counts", "n_features", "percent_mt")) {
    ledger[[paste0("corrected_", metric)]] <- NA_real_
    ledger[[paste0("corrected_", metric)]][preliminary] <- corrected_qc[[metric]]
  }
  ledger$low_counts <- ledger$corrected_n_counts < cfg$qc$min_counts
  ledger$low_features <- ledger$corrected_n_features < cfg$qc$min_features
  # A large number of detected genes is useful doublet/biology context, but
  # scDblFinder is the mandatory multiplet decision. Do not remove cells on
  # this flag alone.
  ledger$high_features_review <- ledger$corrected_n_features > cfg$qc$high_features_review
  ledger$high_mt <- ledger$corrected_percent_mt > cfg$qc$max_mt_percent
  ledger$invalid_qc <- !is.finite(ledger$corrected_percent_mt)
  fail <- ledger$low_counts | ledger$low_features | ledger$high_mt | ledger$invalid_qc
  fail[is.na(fail)] <- TRUE
  ledger$retained <- preliminary & !fail & ledger$doublet_class == "singlet"
  ledger$exclusion_reason <- vapply(seq_len(nrow(ledger)), function(i) {
    if (!preliminary[i]) return("preliminary_qc")
    reasons <- c("low_counts", "low_features", "high_mt", "invalid_qc")
    out <- reasons[vapply(reasons, function(r) isTRUE(ledger[[r]][i]), logical(1))]
    if (ledger$doublet_class[i] == "doublet") out <- c(out, "doublet")
    paste(out, collapse = ";")
  }, character(1))
  if (sum(ledger$retained) < 10) abort(row$library_id, ": fewer than 10 retained singlets; review QC thresholds and correction diagnostics")
  keep <- ledger$barcode[ledger$retained]
  obj <- make_object(corrected[, keep, drop = FALSE], dat$features, row)
  obs <- observed[, keep, drop = FALSE]
  dimnames(obs) <- dimnames(SeuratObject::LayerData(obj, assay = "RNA", layer = "counts"))
  obj[["RNA_observed"]] <- SeuratObject::CreateAssay5Object(counts = obs)
  metadata <- ledger[ledger$retained, , drop = FALSE]
  rownames(metadata) <- metadata$cell_id
  obj <- SeuratObject::AddMetaData(obj, metadata)
  obj@misc$processing <- list(ambient_method = row$ambient_method, ambient_completed = TRUE,
    doublet_method = "scDblFinder", doublets_removed = TRUE, doublet_count_source = "observed",
    RNA_count_source = "ambient_corrected_fractional", seed = cfg$seed)
  # Keep corrected values for all assessed cells for independent correction review.
  list(object = obj, ledger = ledger, corrected_counts = corrected, features = dat$features,
       status = obj@misc$processing)
}
analyze_group <- function(processed, cfg) {
  objects <- lapply(processed, `[[`, "object")
  feature_maps <- lapply(processed, function(p) p$features[order(p$features$feature_id), ])
  if (!all(vapply(feature_maps, identical, logical(1), feature_maps[[1]]))) abort("Feature definitions differ within analysis group; harmonize references explicitly")
  obj <- if (length(objects) == 1) objects[[1]] else merge(objects[[1]], y = objects[-1], merge.data = FALSE)
  obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  obj <- SeuratObject::JoinLayers(obj, assay = "RNA_observed")
  obj <- pca_baseline(obj, cfg$analysis, cfg$seed)
  obj <- cluster_reduction(obj, "pca", "baseline", cfg$analysis, cfg$seed)
  batch_variable <- cfg$integration$batch_variable
  batch <- obj[[batch_variable]][, 1]
  integration_enabled <- switch(cfg$integration$enabled,
    auto = length(unique(batch)) > 1, always = TRUE, never = FALSE)
  integration_status <- if (integration_enabled) "requested" else "single_sequencing_run"
  active <- "baseline_clusters"
  if (integration_enabled) {
    if (length(unique(batch)) < 2)
      abort(cfg$integration$method, " requested but ", batch_variable, " has fewer than two values")
    design <- unique(data.frame(batch = batch, condition = obj$condition))
    if (length(unique(design$condition)) > 1 && all(vapply(split(design$condition, design$batch), function(z) length(unique(z)) == 1, logical(1))))
      warning("Batch and condition are confounded; Harmony is applied by pipeline policy, so between-condition effects require extra care", call. = FALSE)
  }
  if (integration_enabled && cfg$integration$method == "harmony") {
    need("harmony")
    obj <- harmony::RunHarmony(obj, group.by.vars = cfg$integration$batch_variable,
      reduction.use = "pca", reduction.save = "integrated.harmony",
      plot_convergence = FALSE, verbose = FALSE)
    obj <- cluster_reduction(obj, "integrated.harmony", "harmony", cfg$analysis, cfg$seed)
    integration_status <- "harmony_completed_requires_biological_review"
    active <- "harmony_clusters"
  }
  if (integration_enabled && cfg$integration$method == "rpca") {
    # Analyze a copy; retain the unintegrated PCA and embedding in the returned object.
    integrated <- obj
    integrated[["RNA"]] <- split(integrated[["RNA"]], f = batch)
    integrated <- pca_baseline(integrated, cfg$analysis, cfg$seed)
    dims <- seq_len(min(ncol(SeuratObject::Embeddings(integrated[["pca"]])), min(table(batch)) - 2L))
    if (length(dims) < 2) abort("Insufficient cells per batch for RPCA")
    integrated <- Seurat::IntegrateLayers(integrated, method = Seurat::RPCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.rpca", dims = dims,
      k.weight = min(100L, min(table(batch)) - 1L), verbose = FALSE)
    obj[["integrated.rpca"]] <- integrated[["integrated.rpca"]]
    obj <- cluster_reduction(obj, "integrated.rpca", "rpca", cfg$analysis, cfg$seed)
    integration_status <- "rpca_completed_requires_biological_review"
    active <- "rpca_clusters"
  }
  SeuratObject::Idents(obj) <- active
  obj$annotation_broad <- "unassigned"
  obj$annotation_fine <- "unassigned"
  obj@misc$pipeline <- list(integration = integration_status,
    integration_method = if (integration_enabled) cfg$integration$method else "none",
    integration_batch_variable = cfg$integration$batch_variable, active_clusters = active,
    ambient_required = TRUE, doublet_removal_required = TRUE, count_source = "ambient_corrected_fractional")
  markers <- data.frame()
  if (cfg$markers$enabled) {
    markers <- Seurat::FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox",
      min.pct = cfg$markers$min_pct, logfc.threshold = cfg$markers$logfc_threshold, verbose = FALSE)
    if (nrow(markers)) {
      fmap <- obj[["RNA"]][[]]
      markers$feature_id <- fmap[markers$gene, "feature_id"]
      markers$symbol <- fmap[markers$gene, "symbol"]
    }
  }
  review <- data.frame(cluster = sort(unique(as.character(obj[[active]][, 1]))),
    proposed_broad = "", proposed_fine = "", accepted_broad = "", accepted_fine = "",
    evidence = "", reviewer = "", review_date = "")
  list(object = obj, markers = markers, annotation_review = review)
}
