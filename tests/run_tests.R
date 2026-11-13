#!/usr/bin/env Rscript
testthat::test_dir("tests", reporter = "summary", stop_on_failure = TRUE)
