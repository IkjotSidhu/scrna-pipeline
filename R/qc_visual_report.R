save_qc_plot <- function(plot, path, width = 10, height = 6) {
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = 150, bg = "white", limitsize = FALSE)
}

add_qc_groups <- function(qc, run_dir) {
  manifest <- utils::read.csv(file.path(run_dir, "samples.csv"), stringsAsFactors = FALSE)
  index <- match(qc$library_id, manifest$library_id)
  for (column in c("sample_id", "batch_id")) {
    if (!column %in% names(qc) || anyNA(qc[[column]])) qc[[column]] <- manifest[[column]][index]
  }
  if (anyNA(qc$sample_id) || anyNA(qc$batch_id)) abort("Could not map all cells to sample_id and batch_id from samples.csv")
  qc
}

render_qc_visuals <- function(qc, cfg, group, group_title, output_dir, report_path) {
  qc$plot_group <- factor(qc[[group]], levels = sort(unique(qc[[group]])))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  theme <- ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"), strip.text = ggplot2::element_text(face = "bold"), legend.position = "none")
  r2 <- vapply(split(qc, qc$plot_group), function(x) if (nrow(x) < 3 || length(unique(x$n_counts)) < 2) NA_real_ else summary(stats::lm(n_features ~ log10(n_counts), data = x))$r.squared, numeric(1))
  labels <- data.frame(plot_group = factor(names(r2), levels = levels(qc$plot_group)), x = Inf, y = -Inf, label = paste0("R² = ", format(round(r2, 2), nsmall = 2)))
  scatter <- ggplot2::ggplot(qc, ggplot2::aes(x = n_counts, y = n_features)) + ggplot2::geom_point(alpha = 0.22, size = 0.45, colour = "#286090") + ggplot2::geom_vline(xintercept = cfg$min_counts, linetype = "dashed", colour = "#b23a48") + ggplot2::geom_hline(yintercept = cfg$min_features, linetype = "dashed", colour = "#b23a48") + ggplot2::geom_hline(yintercept = cfg$high_features_review, linetype = "dotted", colour = "#8a5a00") + ggplot2::geom_label(data = labels, ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = -0.5, size = 3, linewidth = 0) + ggplot2::scale_x_log10(labels = scales::label_comma()) + ggplot2::scale_y_continuous(labels = scales::label_comma()) + ggplot2::facet_wrap(~plot_group, scales = "free") + theme + ggplot2::labs(title = paste("UMIs versus detected genes by", group_title), subtitle = "Dashed red lines: removal thresholds. Dotted amber line: high-gene review flag only.", x = "UMIs per cell (log scale)", y = "Detected genes per cell")
  save_qc_plot(scatter, file.path(output_dir, "umi-vs-genes.png"), 11, 5)
  histogram <- function(metric, title, x_label, threshold, filename, log_x = FALSE) {
    p <- ggplot2::ggplot(qc, ggplot2::aes(x = .data[[metric]])) + ggplot2::geom_histogram(bins = 60, fill = "#5b8db8", colour = "white", linewidth = 0.2) + ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", colour = "#b23a48", linewidth = 0.7) + ggplot2::facet_wrap(~plot_group, scales = "free_y") + theme + ggplot2::labs(title = paste(title, "by", group_title), subtitle = "Dashed red line: configured removal threshold.", x = x_label, y = "Cells")
    if (log_x) p <- p + ggplot2::scale_x_log10(labels = scales::label_comma()) else p <- p + ggplot2::scale_x_continuous(labels = scales::label_comma())
    save_qc_plot(p, file.path(output_dir, filename), 11, 5)
  }
  histogram("n_counts", "UMI distribution", "UMIs per cell (log scale)", cfg$min_counts, "umi-distribution.png", TRUE)
  histogram("n_features", "Detected-gene distribution", "Detected genes per cell", cfg$min_features, "gene-distribution.png")
  histogram("percent_mt", "Mitochondrial RNA distribution", "Mitochondrial RNA (%)", cfg$max_mt_percent, "mito-distribution.png")
  violin <- function(metric, title, y_label, threshold, filename) {
    p <- ggplot2::ggplot(qc, ggplot2::aes(x = plot_group, y = .data[[metric]], fill = plot_group)) + ggplot2::geom_violin(trim = FALSE, colour = "#34536f", alpha = 0.8) + ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", colour = "#b23a48", linewidth = 0.7) + theme + ggplot2::labs(title = paste(title, "by", group_title), subtitle = "Dashed red line: configured removal threshold. No boxplot is overlaid.", x = group_title, y = y_label)
    save_qc_plot(p, file.path(output_dir, filename), 10, 5)
  }
  violin("n_counts", "UMIs per cell", "UMIs per cell", cfg$min_counts, "umi-violin.png")
  violin("n_features", "Detected genes per cell", "Detected genes per cell", cfg$min_features, "gene-violin.png")
  violin("percent_mt", "Mitochondrial RNA per cell", "Mitochondrial RNA (%)", cfg$max_mt_percent, "mito-violin.png")
  image <- function(filename, alt) paste0('<figure><img src="', basename(output_dir), '/', filename, '" alt="', html_escape(alt), '"><figcaption>', html_escape(alt), '</figcaption></figure>')
  header <- paste0('<!doctype html><html><head><meta charset="utf-8"><title>Exploratory QC plots</title><style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#192b3b}img{max-width:100%;border:1px solid #d7e0e8}figure{margin:28px 0}figcaption{color:#52616d;margin-top:7px}.note{background:#eff5fa;padding:16px}h2{margin-top:42px}</style></head><body>')
  html <- paste0(header, '<h1>Exploratory QC plots by ', html_escape(group_title), '</h1><div class="note">These plots use ambient-corrected, assessed cells before final filtering so that low-quality cells and potential outliers remain visible. They are a diagnostic companion to the table-first <a href="report-bioinformatics.html">bioinformatics report</a>.</div><h2>Relationship between sequencing depth and complexity</h2><p>R² comes from a linear model of detected genes against log10 UMIs within each ', html_escape(group_title), '. It describes the trend only; it is not a QC pass/fail score.</p>', image("umi-vs-genes.png", "UMIs versus detected genes per cell"), '<h2>Distributions</h2>', image("umi-distribution.png", "UMI distribution"), image("gene-distribution.png", "Detected-gene distribution"), image("mito-distribution.png", "Mitochondrial RNA distribution"), '<h2>Simple violin plots</h2><p>Each violin shows one metric only. There are deliberately no overlaid boxplots or competing axes.</p>', image("umi-violin.png", "UMIs per cell"), image("gene-violin.png", "Detected genes per cell"), image("mito-violin.png", "Mitochondrial RNA per cell"), '</body></html>')
  writeLines(html, report_path); report_path
}

write_qc_visual_report <- function(run_dir, data_file = "report_data.rds") {
  need(c("ggplot2", "scales")); bundle <- readRDS(file.path(run_dir, data_file))
  if (is.null(bundle$config$qc)) abort("This run has no saved QC configuration")
  qc <- bundle$qc; stage <- if ("corrected" %in% qc$stage) "corrected" else "original"
  qc <- qc[qc$stage == stage & is.finite(qc$n_counts) & is.finite(qc$n_features) & is.finite(qc$percent_mt) & qc$n_counts > 0, , drop = FALSE]
  if (!nrow(qc)) abort("No QC cells are available for visual review")
  qc <- add_qc_groups(qc, run_dir)
  sample_report <- render_qc_visuals(qc, bundle$config$qc, "sample_id", "biological sample", file.path(run_dir, "qc-visuals-by-sample"), file.path(run_dir, "report-qc-by-sample.html"))
  run_report <- render_qc_visuals(qc, bundle$config$qc, "batch_id", "sequencing run", file.path(run_dir, "qc-visuals-by-sequencing-run"), file.path(run_dir, "report-qc-by-sequencing-run.html"))
  path <- file.path(run_dir, "report-qc-visuals.html")
  writeLines(paste0('<!doctype html><html><head><meta charset="utf-8"><title>Exploratory QC plots</title></head><body><h1>Exploratory QC plots</h1><ul><li><a href="', basename(sample_report), '">QC plots by biological sample</a></li><li><a href="', basename(run_report), '">QC plots by sequencing run</a></li></ul></body></html>'), path)
  invisible(c(sample_report, run_report, path))
}
