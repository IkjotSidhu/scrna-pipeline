# Configuration and manifests are data, never sourced as executable R.
abort <- function(...) stop(..., call. = FALSE)
need <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) abort("Missing packages: ", paste(missing, collapse = ", "),
                             ". Install separately; analysis never installs packages.")
}
absolute <- function(path, root) {
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(root, path)
  normalizePath(path, mustWork = FALSE)
}
merge_settings <- function(default, supplied, prefix = "config") {
  if (!is.list(supplied) || is.null(names(supplied)) || anyDuplicated(names(supplied)))
    abort(prefix, " must be a mapping with unique keys")
  unknown <- setdiff(names(supplied), names(default))
  if (length(unknown)) abort("Unknown setting: ", prefix, ".", paste(unknown, collapse = ", "))
  for (key in names(supplied)) {
    if (is.list(default[[key]])) default[[key]] <- merge_settings(default[[key]], supplied[[key]], paste0(prefix, ".", key))
    else default[key] <- supplied[key]
  }
  default
}
number <- function(x, name, low = 0, high = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x) || x < low || x > high ||
      (integer && x != floor(x))) abort("Invalid numeric setting: ", name)
}
choice <- function(x, name, choices) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !x %in% choices)
    abort(name, " must be one of: ", paste(choices, collapse = ", "))
}
read_config <- function(path) {
  need(c("yaml", "jsonlite", "digest", "Matrix", "Seurat", "ggplot2"))
  if (utils::packageVersion("Seurat") < "5.1.0") abort("Seurat >= 5.1.0 required")
  path <- normalizePath(path, mustWork = TRUE)
  defaults <- list(project = "scrna", samples = "samples.csv", output = "../results",
    seed = 42L, workers = 1L, future_max_gb = 4,
    qc = list(min_counts = 500, min_features = 200, high_features_review = 10000,
              max_mt_percent = 25, preliminary_min_features = 50),
    ambient = list(method = "auto"),
    doublets = list(method = "scDblFinder", expected_rate = 0.0),
    analysis = list(nfeatures = 2000L, npcs = 30L, resolution = 0.5,
                    neighbors = 20L, umap_neighbors = 30L),
    integration = list(enabled = "auto", method = "harmony", batch_variable = "batch_id"),
    markers = list(enabled = TRUE, min_pct = 0.1, logfc_threshold = 0.25))
  cfg <- merge_settings(defaults, yaml::read_yaml(path))
  for (key in c("project", "samples", "output"))
    if (!is.character(cfg[[key]]) || length(cfg[[key]]) != 1 || !nzchar(cfg[[key]])) abort("Invalid ", key)
  for (key in c("seed", "workers")) number(cfg[[key]], key, if (key == "seed") 0 else 1, integer = TRUE)
  number(cfg$future_max_gb, "future_max_gb", 0.1)
  for (key in names(cfg$qc)) number(cfg$qc[[key]], paste0("qc.", key), 0)
  if (cfg$qc$max_mt_percent > 100 || cfg$qc$high_features_review < cfg$qc$min_features ||
      cfg$qc$preliminary_min_features > cfg$qc$min_features) abort("Inconsistent QC thresholds")
  for (key in c("nfeatures", "npcs", "neighbors", "umap_neighbors"))
    number(cfg$analysis[[key]], paste0("analysis.", key), 2, integer = TRUE)
  number(cfg$analysis$resolution, "analysis.resolution", 0.001)
  choice(cfg$ambient$method, "ambient.method", c("auto", "soupx", "decontx"))
  choice(cfg$doublets$method, "doublets.method", c("scDblFinder"))
  if (is.logical(cfg$integration$enabled) && length(cfg$integration$enabled) == 1 && !is.na(cfg$integration$enabled))
    cfg$integration$enabled <- if (cfg$integration$enabled) "always" else "never"
  choice(cfg$integration$enabled, "integration.enabled", c("auto", "always", "never"))
  choice(cfg$integration$method, "integration.method", c("harmony", "rpca"))
  choice(cfg$integration$batch_variable, "integration.batch_variable",
         c("batch_id", "library_id", "capture_id", "sample_id"))
  number(cfg$doublets$expected_rate, "doublets.expected_rate", 0, 0.5)
  for (v in list(cfg$markers$enabled))
    if (!is.logical(v) || length(v) != 1 || is.na(v)) abort("Expected a boolean setting")
  number(cfg$markers$min_pct, "markers.min_pct", 0, 1)
  number(cfg$markers$logfc_threshold, "markers.logfc_threshold", 0)
  cfg$samples <- absolute(cfg$samples, dirname(path))
  cfg$output <- absolute(cfg$output, dirname(path))
  cfg$config_path <- path
  cfg
}
read_manifest <- function(cfg) {
  tab <- read.csv(cfg$samples, colClasses = "character", check.names = FALSE, na.strings = "",
                  strip.white = TRUE)
  if (!"dataset_role" %in% names(tab)) tab$dataset_role <- "study"
  required <- c("library_id", "sample_id", "donor_id", "species", "tissue", "condition",
                "batch_id", "reference_id", "analysis_group", "filtered_matrix", "raw_matrix",
                "cell_metadata", "ambient_method")
  if (anyDuplicated(names(tab))) abort("Duplicate manifest columns")
  if (length(setdiff(required, names(tab)))) abort("Missing manifest columns: ", paste(setdiff(required, names(tab)), collapse = ", "))
  if (!nrow(tab)) abort("Empty manifest")
  for (col in setdiff(required, c("raw_matrix", "cell_metadata", "ambient_method")))
    if (anyNA(tab[[col]]) || any(!nzchar(trimws(tab[[col]])))) abort("Missing values in ", col)
  if (anyDuplicated(tab$library_id)) abort("library_id must be unique: one row per physical capture")
  for (col in c("library_id", "analysis_group"))
    if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9.-]*$", tab[[col]]))) abort(col, " must use letters, numbers, dots or hyphens")
  if (any(!tab$species %in% c("human", "mouse"))) abort("species must be human or mouse")
  if (anyNA(tab$dataset_role) || any(!tab$dataset_role %in%
      c("study", "public_benchmark", "internal_benchmark")))
    abort("dataset_role must be study, public_benchmark, or internal_benchmark")
  tab$ambient_method[is.na(tab$ambient_method)] <- cfg$ambient$method
  if (any(!tab$ambient_method %in% c("auto", "soupx", "decontx"))) abort("Invalid per-library ambient_method")
  for (col in c("filtered_matrix", "raw_matrix", "cell_metadata", "cellranger_metrics")) {
    if (!col %in% names(tab)) tab[[col]] <- NA_character_
    present <- !is.na(tab[[col]])
    tab[[col]][present] <- vapply(tab[[col]][present], absolute, character(1), root = dirname(cfg$samples))
    if (any(!file.exists(tab[[col]][present]))) abort("Missing input path in ", col)
  }
  if (anyDuplicated(tab$filtered_matrix)) abort("The same filtered input is listed more than once")
  auto <- tab$ambient_method == "auto"
  tab$ambient_method[auto] <- ifelse(is.na(tab$raw_matrix[auto]), "decontx", "soupx")
  if (any(tab$ambient_method == "soupx" & is.na(tab$raw_matrix))) abort("SoupX requested without raw_matrix; supply raw droplets or select decontx")
  for (group in split(tab, tab$analysis_group)) {
    for (col in c("species", "reference_id"))
      if (length(unique(group[[col]])) != 1) abort("Mixed ", col, " in analysis_group ", group$analysis_group[1])
    if (length(unique(group$ambient_method)) > 1)
      warning("Mixed ambient correction within group ", group$analysis_group[1], "; inspect correction effects before comparisons", call. = FALSE)
  }
  for (s in split(tab, tab$sample_id))
    for (col in c("donor_id", "species", "tissue", "condition"))
      if (length(unique(s[[col]])) != 1) abort("Inconsistent ", col, " for sample ", s$sample_id[1])
  tab
}
