# 10x v3 H5 and v2/v3 Matrix Market readers retain IDs and symbols separately.
input_files <- function(path) {
  if (!dir.exists(path)) return(normalizePath(path, mustWork = TRUE))
  choose_file <- function(stems) {
    candidates <- file.path(path, c(stems, paste0(stems, ".gz")))
    found <- candidates[file.exists(candidates)]
    if (length(found) != 1) abort("Expected exactly one of ", paste(stems, collapse = "/"), " (optionally gzipped) in ", path)
    normalizePath(found)
  }
  c(matrix = choose_file("matrix.mtx"), features = choose_file(c("features.tsv", "genes.tsv")),
    barcodes = choose_file("barcodes.tsv"))
}
read_tsv <- function(path) read.delim(path, header = FALSE, quote = "", comment.char = "", colClasses = "character")
read_10x <- function(path) {
  if (dir.exists(path)) {
    files <- input_files(path)
    con <- if (grepl("\\.gz$", files["matrix"])) gzfile(files["matrix"]) else file(files["matrix"])
    counts <- methods::as(Matrix::readMM(con), "CsparseMatrix")
    f <- read_tsv(files["features"])
    if (ncol(f) < 2) abort("Feature table needs ID and symbol columns")
    features <- data.frame(feature_id = f[[1]], symbol = f[[2]],
                           feature_type = if (ncol(f) >= 3) f[[3]] else "Gene Expression")
    barcodes <- read_tsv(files["barcodes"])[[1]]
  } else {
    need("hdf5r")
    h <- hdf5r::H5File$new(path, mode = "r")
    on.exit(h$close_all(), add = TRUE)
    required <- c("matrix/data", "matrix/indices", "matrix/indptr", "matrix/shape",
                  "matrix/barcodes", "matrix/features/id", "matrix/features/name", "matrix/features/feature_type")
    if (!all(vapply(required, h$exists, logical(1)))) abort("Expected 10x v3 H5 schema: ", path)
    shape <- as.integer(h[["matrix/shape"]][])
    counts <- Matrix::sparseMatrix(i = as.integer(h[["matrix/indices"]][]),
      p = as.integer(h[["matrix/indptr"]][]), x = as.numeric(h[["matrix/data"]][]),
      dims = shape, index1 = FALSE)
    features <- data.frame(feature_id = h[["matrix/features/id"]][], symbol = h[["matrix/features/name"]][],
                           feature_type = h[["matrix/features/feature_type"]][])
    barcodes <- h[["matrix/barcodes"]][]
  }
  if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes)) abort("Input dimensions do not match feature/barcode tables")
  if (anyNA(barcodes) || any(!nzchar(barcodes)) || anyDuplicated(barcodes)) abort("Invalid/duplicate barcodes")
  excluded <- table(features$feature_type[features$feature_type != "Gene Expression"])
  keep <- features$feature_type == "Gene Expression"
  features <- features[keep, , drop = FALSE]
  counts <- counts[keep, , drop = FALSE]
  if (!nrow(counts) || !ncol(counts)) abort("Empty Gene Expression matrix")
  if (anyNA(features$feature_id) || any(!nzchar(features$feature_id)) || anyDuplicated(features$feature_id)) abort("Gene Expression feature IDs must be unique and nonempty")
  if (anyNA(features$symbol) || any(!nzchar(features$symbol))) abort("Missing gene symbols; provide annotated 10x features")
  if (any(!is.finite(counts@x)) || any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)) abort("Expected finite, nonnegative, integer observed UMI counts")
  rownames(counts) <- features$feature_id
  colnames(counts) <- barcodes
  list(counts = counts, features = features, excluded_modalities = as.list(excluded))
}
qc_metrics <- function(counts, features, species) {
  mt <- grepl(if (species == "human") "^MT-" else "^mt-", features$symbol)
  if (!any(mt)) abort("No mitochondrial symbols found for ", species, "; check species and feature annotation")
  n <- Matrix::colSums(counts)
  detected <- Matrix::colSums(counts > 0)
  data.frame(barcode = colnames(counts), n_counts = n, n_features = detected,
    percent_mt = ifelse(n > 0, Matrix::colSums(counts[mt, , drop = FALSE]) / n * 100, NA_real_),
    genes_per_umi = ifelse(n > 0, detected / n, NA_real_), row.names = colnames(counts))
}
import_library <- function(row) {
  dat <- read_10x(row$filtered_matrix)
  dat$qc <- qc_metrics(dat$counts, dat$features, row$species)
  if (is.na(row$cell_metadata)) {
    dat$cell_metadata <- data.frame(
      barcode = colnames(dat$counts),
      capture_id = row$library_id,
      stringsAsFactors = FALSE
    )
  } else {
    meta <- read.csv(row$cell_metadata, colClasses = "character", check.names = FALSE,
                     na.strings = "", strip.white = TRUE)
    required <- c("barcode", "capture_id")
    if (length(setdiff(required, names(meta))))
      abort("cell_metadata needs columns: ", paste(required, collapse = ", "))
    if (anyNA(meta$barcode) || anyNA(meta$capture_id) || any(!nzchar(meta$barcode)) ||
        any(!nzchar(meta$capture_id)) || anyDuplicated(meta$barcode))
      abort("Invalid barcode/capture_id values in ", row$cell_metadata)
    if (!setequal(meta$barcode, colnames(dat$counts)))
      abort("cell_metadata must contain every filtered barcode exactly once for ", row$library_id)
    dat$cell_metadata <- meta[match(colnames(dat$counts), meta$barcode), , drop = FALSE]
  }
  dat$row <- row
  dat
}
align_raw <- function(raw, filtered) {
  # Cell Ranger sample-level filtered matrices can contain a documented subset
  # of the parent raw matrix's gene features (for example, 10x Flex data).
  # Align by stable feature ID rather than requiring identical feature sets or
  # trusting row order. Every filtered feature must still be present in raw.
  missing_features <- setdiff(rownames(filtered$counts), rownames(raw$counts))
  if (length(missing_features))
    abort("Raw matrix is missing ", length(missing_features), " filtered gene IDs")
  if (!all(colnames(filtered$counts) %in% colnames(raw$counts))) abort("Filtered barcodes missing from raw matrix")
  raw$counts <- raw$counts[rownames(filtered$counts), , drop = FALSE]
  delta <- raw$counts[, colnames(filtered$counts), drop = FALSE] - filtered$counts
  if (any(delta@x != 0)) abort("Raw and filtered counts disagree; provide matched unmodified Cell Ranger outputs")
  attr(raw$counts, "raw_features_excluded") <- nrow(raw$features) - nrow(filtered$features)
  raw$counts
}
read_cellranger_metrics <- function(path) {
  empty <- data.frame(category = character(), library_type = character(), grouped_by = character(),
    group_name = character(), metric_name = character(), metric_value = character(), stringsAsFactors = FALSE)
  if (is.na(path)) return(empty)
  tab <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")
  tidy_names <- c("Category", "Library Type", "Grouped By", "Group Name", "Metric Name", "Metric Value")
  if (all(tidy_names %in% names(tab))) {
    tab <- tab[, tidy_names, drop = FALSE]
    names(tab) <- names(empty)
    return(tab)
  }
  if (nrow(tab) != 1L) abort("Unsupported Cell Ranger metrics format in ", path)
  data.frame(category = "Library", library_type = "Gene Expression", grouped_by = "",
    group_name = "", metric_name = names(tab), metric_value = as.character(tab[1, ]),
    stringsAsFactors = FALSE)
}
