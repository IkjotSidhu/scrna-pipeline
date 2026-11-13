source(file.path("..", "R", "config.R"))

make_manifest_fixture <- function(root) {
  dir.create(file.path(root, "ileum"), recursive = TRUE)
  dir.create(file.path(root, "duodenum"), recursive = TRUE)
  data.frame(
    library_id = c("ileum", "duodenum"), sample_id = c("s1", "s2"),
    donor_id = "d1", species = "mouse", tissue = c("ileum", "duodenum"),
    condition = "control", batch_id = "b1", reference_id = "mm10",
    analysis_group = "small-intestine", filtered_matrix = c("ileum", "duodenum"),
    raw_matrix = NA_character_, cell_metadata = NA_character_, ambient_method = "decontx"
  )
}

testthat::test_that("compatible tissues can share one analysis group", {
  root <- tempfile("scrna-config-")
  dir.create(root)
  manifest <- make_manifest_fixture(root)
  write.csv(manifest, file.path(root, "samples.csv"), row.names = FALSE, na = "")
  yaml::write_yaml(list(samples = "samples.csv", output = "results"), file.path(root, "project.yml"))
  cfg <- read_config(file.path(root, "project.yml"))
  observed <- read_manifest(cfg)
  testthat::expect_identical(unique(observed$analysis_group), "small-intestine")
  testthat::expect_setequal(observed$tissue, c("ileum", "duodenum"))
  testthat::expect_identical(cfg$integration$enabled, "auto")
  testthat::expect_identical(cfg$integration$method, "harmony")
})

testthat::test_that("species cannot be mixed within an analysis group", {
  root <- tempfile("scrna-config-")
  dir.create(root)
  manifest <- make_manifest_fixture(root)
  manifest$species[2] <- "human"
  write.csv(manifest, file.path(root, "samples.csv"), row.names = FALSE, na = "")
  yaml::write_yaml(list(samples = "samples.csv", output = "results"), file.path(root, "project.yml"))
  cfg <- read_config(file.path(root, "project.yml"))
  testthat::expect_error(read_manifest(cfg), "Mixed species")
})
