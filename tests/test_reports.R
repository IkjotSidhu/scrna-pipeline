source(file.path("..", "R", "config.R"))
source(file.path("..", "R", "report.R"))

testthat::test_that("analyst and biologist reports render from saved data", {
  root <- tempfile("scrna-report-")
  dir.create(root)
  summary <- data.frame(library_id = c("ileum", "duodenum"),
    tissue = c("ileum", "duodenum"), input_cells = c(100L, 120L), retained = c(80L, 90L))
  qc <- expand.grid(library_id = c("OB1", "OB2"), stage = c("original", "retained"), cell = 1:4)
  qc$n_counts <- seq_len(nrow(qc)) * 100
  qc$n_features <- seq_len(nrow(qc)) * 20
  qc$percent_mt <- seq_len(nrow(qc))
  qc$dataset_role <- "study"
  qc$tissue <- "small_intestine"
  embedding <- data.frame(UMAP_1 = c(-1, 0, 1, 2), UMAP_2 = c(0, 1, 0, -1),
    library_id = c("ileum", "ileum", "duodenum", "duodenum"),
    capture_id = c("OB1", "OB2", "OB3", "OB4"), sample_id = c("s1", "s1", "s2", "s2"),
    tissue = c("ileum", "ileum", "duodenum", "duodenum"), condition = "control",
    cluster = c("0", "0", "1", "1"))
  config <- list(qc = list(min_counts = 500, min_features = 200, high_features_review = 10000, max_mt_percent = 25))
  summary$assessed <- c(100L, 120L)
  summary$predicted_doublets <- c(5L, 6L)
  summary$high_features_review <- c(2L, 3L)
  summary$high_mt <- c(4L, 5L)
  saveRDS(list(project = "fixture", mode = "run", config = config, summary = summary, qc = qc,
    cellranger_metrics = data.frame(library_id = "ileum", category = "Library",
      library_type = "Gene Expression", grouped_by = "", group_name = "",
      metric_name = "Median genes per cell", metric_value = "1,000"),
    embeddings = list(small_intestine_baseline = embedding)), file.path(root, "report_data.rds"))
  jsonlite::write_json(list(status = "completed"), file.path(root, "status.json"), auto_unbox = TRUE)
  write_reports(root)
  testthat::expect_true(all(file.exists(file.path(root,
    c("report.html", "report-bioinformatics.html", "report-biologist.html")))))
  testthat::expect_match(paste(readLines(file.path(root, "report.html")), collapse = " "),
                         "Bioinformatics report")
  comparison <- read.csv(file.path(root, "qc_comparison.csv"), check.names = FALSE)
  testthat::expect_equal(nrow(comparison), 2L)
  testthat::expect_true(all(c("Mean counts", "Median features", "Role") %in% names(comparison)))
  testthat::expect_true(file.exists(file.path(root, "cellranger_qc.csv")))
  testthat::expect_true(file.exists(file.path(root, "cell_retention.csv")))
  testthat::expect_true(all(file.exists(file.path(root, c("qc_policy.csv", "quality_flags.csv")))))
  analyst <- paste(readLines(file.path(root, "report-bioinformatics.html")), collapse = " ")
  testthat::expect_match(analyst, "Final-cell QC comparison")
  testthat::expect_match(analyst, "Cell retention through filtering")
  testthat::expect_match(analyst, "QC policy used for this run")
  testthat::expect_match(analyst, "high-gene rule is diagnostic only", ignore.case = TRUE)
  testthat::expect_match(analyst, "10x sequencing and library QC")
  testthat::expect_false(grepl("qc_n_counts.png|qc_n_features.png", analyst))
  testthat::expect_false(grepl("<img|UMAP", analyst))
  biological <- paste(readLines(file.path(root, "report-biologist.html")), collapse = " ")
  testthat::expect_false(grepl("<img", biological))
})
