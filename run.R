#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
root <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
if (nzchar(Sys.getenv("SCRNA_R_LIB"))) .libPaths(c(Sys.getenv("SCRNA_R_LIB"), .libPaths()))
for (file in c("config.R", "import.R", "process.R", "report.R", "qc_visual_report.R", "downstream_diagnostics.R", "engine.R")) source(file.path(root, "R", file))
if (length(args) != 2 || !args[1] %in% c("qc", "run", "report", "qc-plots", "diagnostics")) {
  cat("Usage: Rscript run.R qc|run config/project.yml\n       Rscript run.R report|qc-plots|diagnostics /path/to/run-directory\n")
  quit(status = 2)
}
if (args[1] == "report") {
  reports <- write_reports(normalizePath(args[2], mustWork = TRUE))
  cat("Reports:", paste(reports, collapse = ", "), "\n")
} else if (args[1] == "qc-plots") {
  report <- write_qc_visual_report(normalizePath(args[2], mustWork = TRUE))
  cat("QC visual report:", report, "\n")
} else if (args[1] == "diagnostics") {
  report <- write_downstream_diagnostics(normalizePath(args[2], mustWork = TRUE))
  cat("Downstream diagnostics:", report, "\n")
} else {
  captured_warnings <- character()
  run_dir <- withCallingHandlers(
    run_pipeline(args[2], args[1], root),
    warning = function(w) captured_warnings <<- c(captured_warnings, conditionMessage(w))
  )
  writeLines(unique(captured_warnings), file.path(run_dir, "warnings.log"))
}
