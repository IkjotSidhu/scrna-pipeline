source(file.path("..", "R", "config.R"))
source(file.path("..", "R", "import.R"))

make_10x_fixture <- function(root) {
  dir.create(root, recursive = TRUE)
  counts <- Matrix::Matrix(matrix(c(1, 0, 2, 3, 1, 0), nrow = 2), sparse = TRUE)
  Matrix::writeMM(counts, file.path(root, "matrix.mtx"))
  write.table(
    data.frame(
      c("ENSMUSG1", "ENSMUSG2"),
      c("mt-Nd1", "Actb"),
      c("Gene Expression", "Gene Expression")
    ),
    file.path(root, "features.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  write.table(
    c("AA-1", "BB-1", "CC-1"), file.path(root, "barcodes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
}

testthat::test_that("barcode metadata is aligned to matrix column order", {
  root <- tempfile("scrna-fixture-")
  make_10x_fixture(root)
  metadata <- file.path(root, "cells.csv")
  write.csv(
    data.frame(
      barcode = c("CC-1", "AA-1", "BB-1"),
      capture_id = c("OB2", "OB1", "OB1")
    ),
    metadata, row.names = FALSE
  )
  row <- data.frame(
    library_id = "lib1", sample_id = "s1", donor_id = "d1", species = "mouse",
    tissue = "gut", condition = "pilot", batch_id = "b1", reference_id = "mm10",
    analysis_group = "g1", filtered_matrix = root, raw_matrix = NA_character_,
    cell_metadata = metadata, ambient_method = "decontx"
  )
  dat <- import_library(row)
  testthat::expect_identical(dat$cell_metadata$barcode, colnames(dat$counts))
  testthat::expect_identical(dat$cell_metadata$capture_id, c("OB1", "OB1", "OB2"))
})

testthat::test_that("incomplete barcode metadata is rejected", {
  root <- tempfile("scrna-fixture-")
  make_10x_fixture(root)
  metadata <- file.path(root, "cells.csv")
  write.csv(data.frame(barcode = c("AA-1", "BB-1"), capture_id = "OB1"),
            metadata, row.names = FALSE)
  row <- data.frame(
    library_id = "lib1", sample_id = "s1", donor_id = "d1", species = "mouse",
    tissue = "gut", condition = "pilot", batch_id = "b1", reference_id = "mm10",
    analysis_group = "g1", filtered_matrix = root, raw_matrix = NA_character_,
    cell_metadata = metadata, ambient_method = "decontx"
  )
  testthat::expect_error(import_library(row), "every filtered barcode exactly once")
})

testthat::test_that("raw feature supersets align by stable gene ID", {
  filtered_counts <- methods::as(Matrix::Matrix(
    matrix(c(2, 0, 1, 3), nrow = 2,
      dimnames = list(c("gene_b", "gene_a"), c("cell_1", "cell_2"))),
    sparse = TRUE), "dgCMatrix")
  raw_counts <- methods::as(Matrix::Matrix(
    matrix(c(3, 1, 9, 0, 2, 8), nrow = 3,
      dimnames = list(c("gene_a", "gene_b", "raw_only"), c("cell_2", "cell_1"))),
    sparse = TRUE), "dgCMatrix")
  filtered <- list(counts = filtered_counts,
    features = data.frame(feature_id = rownames(filtered_counts)))
  raw <- list(counts = raw_counts,
    features = data.frame(feature_id = rownames(raw_counts)))
  aligned <- align_raw(raw, filtered)
  testthat::expect_identical(rownames(aligned), rownames(filtered_counts))
  testthat::expect_identical(colnames(aligned), colnames(raw_counts))
  testthat::expect_equal(attr(aligned, "raw_features_excluded"), 1L)
})

testthat::test_that("raw matrices missing filtered features are rejected", {
  filtered_counts <- methods::as(Matrix::Matrix(matrix(1, nrow = 2,
    dimnames = list(c("gene_a", "gene_b"), "cell_1")), sparse = TRUE), "dgCMatrix")
  raw_counts <- methods::as(Matrix::Matrix(matrix(1, nrow = 1,
    dimnames = list("gene_a", "cell_1")), sparse = TRUE), "dgCMatrix")
  filtered <- list(counts = filtered_counts,
    features = data.frame(feature_id = rownames(filtered_counts)))
  raw <- list(counts = raw_counts,
    features = data.frame(feature_id = rownames(raw_counts)))
  testthat::expect_error(align_raw(raw, filtered), "missing 1 filtered gene IDs")
})

testthat::test_that("tidy Cell Ranger metrics are preserved", {
  path <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    Category = c("Cells", "Library"), `Library Type` = "Gene Expression",
    `Grouped By` = "", `Group Name` = "",
    `Metric Name` = c("Median genes per cell", "Sequencing saturation"),
    `Metric Value` = c("1,234", "65.2%"), check.names = FALSE), path, row.names = FALSE)
  observed <- read_cellranger_metrics(path)
  testthat::expect_identical(observed$metric_name, c("Median genes per cell", "Sequencing saturation"))
  testthat::expect_identical(observed$metric_value, c("1,234", "65.2%"))
})

testthat::test_that("legacy one-row Cell Ranger metrics are converted to tidy form", {
  path <- tempfile(fileext = ".csv")
  write.csv(data.frame(`Estimated Number of Cells` = "4000", `Reads per Cell` = "20000",
    check.names = FALSE), path, row.names = FALSE)
  observed <- read_cellranger_metrics(path)
  testthat::expect_identical(observed$metric_name, c("Estimated Number of Cells", "Reads per Cell"))
  testthat::expect_identical(observed$metric_value, c("4000", "20000"))
})
