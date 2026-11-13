#!/usr/bin/env Rscript
if (nzchar(Sys.getenv("SCRNA_R_LIB"))) .libPaths(c(Sys.getenv("SCRNA_R_LIB"), .libPaths()))
testthat::test_dir("tests", reporter = "summary", stop_on_failure = TRUE)
