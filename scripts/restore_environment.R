#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/restore_environment.R PROJECT_DIR")
project <- normalizePath(args[1], mustWork = TRUE)
lockfile <- file.path(project, "renv.lock")
if (!file.exists(lockfile)) stop("Missing lockfile: ", lockfile)
renv::consent(provided = TRUE)
library <- renv::paths$library(project = project)
dir.create(library, recursive = TRUE, showWarnings = FALSE)
message("Restoring explicit project library: ", library)
renv::restore(project = project, library = library, lockfile = lockfile,
              clean = TRUE, prompt = FALSE)
cat("Active project library restored at ", library, "\n", sep = "")
