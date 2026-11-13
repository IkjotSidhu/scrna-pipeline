adjusted_rand_index <- function(a, b) {
  tab <- table(a, b); choose2 <- function(x) x * (x - 1) / 2
  total <- sum(choose2(tab)); rows <- sum(choose2(rowSums(tab))); cols <- sum(choose2(colSums(tab)))
  n <- sum(tab); expected <- rows * cols / choose2(n); maximum <- (rows + cols) / 2
  if (maximum == expected) return(NA_real_)
  (total - expected) / (maximum - expected)
}

cluster_resolution_summary <- function(obj, graph_name, resolutions, seed) {
  previous <- NULL
  out <- vector("list", length(resolutions))
  for (i in seq_along(resolutions)) {
    resolution <- resolutions[[i]]
    set.seed(seed)
    clustered <- Seurat::FindClusters(obj, graph.name = graph_name, resolution = resolution,
      random.seed = seed, verbose = FALSE)
    labels <- clustered$seurat_clusters; sizes <- table(labels)
    column <- paste0("resolution_", format(resolution, nsmall = 1, trim = TRUE))
    obj[[column]] <- labels
    row <- data.frame(resolution = resolution, clusters = length(sizes), smallest_cluster = min(sizes),
      median_cluster_size = stats::median(as.integer(sizes)), clusters_under_20 = sum(sizes < 20),
      adjacent_ari = if (is.null(previous)) NA_real_ else adjusted_rand_index(previous, labels))
    previous <- labels; out[[i]] <- row
  }
  list(summary = do.call(rbind, out), object = obj)
}

write_downstream_diagnostics <- function(run_dir) {
  need(c("ggplot2", "scales", "clustree"))
  bundle <- readRDS(file.path(run_dir, "report_data.rds"))
  settings <- bundle$config$analysis
  group_root <- file.path(run_dir, "groups")
  groups <- list.dirs(group_root, recursive = FALSE, full.names = FALSE)
  if (!length(groups)) abort("No analyzed groups found for downstream diagnostics")
  plot_dir <- file.path(run_dir, "downstream-diagnostics")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  theme <- ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"))
  resolutions <- seq(0.1, 2.0, by = 0.1)
  all_resolution <- list(); figures <- character()
  for (group in groups) {
    obj <- readRDS(file.path(group_root, group, "seurat.rds"))
    pcs <- SeuratObject::Embeddings(obj[["pca"]])
    variance <- apply(pcs, 2, stats::var); variance <- variance / sum(variance) * 100
    pc_data <- data.frame(pc = seq_along(variance), variance = variance, cumulative = cumsum(variance), group = group)
    variance_plot <- ggplot2::ggplot(pc_data, ggplot2::aes(pc, variance)) + ggplot2::geom_line(colour = "#286090") + ggplot2::geom_point(colour = "#286090") + ggplot2::geom_vline(xintercept = settings$npcs, linetype = "dashed", colour = "#b23a48") + ggplot2::scale_x_continuous(breaks = pc_data$pc) + theme + ggplot2::labs(title = paste(group, "PC variance"), subtitle = "Dashed red line: PCs used in the current analysis.", x = "Principal component", y = "Variance explained (%)")
    cumulative_plot <- ggplot2::ggplot(pc_data, ggplot2::aes(pc, cumulative)) + ggplot2::geom_line(colour = "#286090") + ggplot2::geom_point(colour = "#286090") + ggplot2::geom_vline(xintercept = settings$npcs, linetype = "dashed", colour = "#b23a48") + ggplot2::scale_x_continuous(breaks = pc_data$pc) + ggplot2::scale_y_continuous(limits = c(0, 100)) + theme + ggplot2::labs(title = paste(group, "cumulative PC variance"), subtitle = "This supports review of the elbow; it is not a standalone rule for choosing PCs.", x = "Principal component", y = "Cumulative variance explained (%)")
    prefix <- gsub("[^A-Za-z0-9.-]", "-", group)
    ggplot2::ggsave(file.path(plot_dir, paste0(prefix, "-pc-variance.png")), variance_plot, width = 9, height = 4.5, dpi = 150, bg = "white")
    ggplot2::ggsave(file.path(plot_dir, paste0(prefix, "-pc-cumulative.png")), cumulative_plot, width = 9, height = 4.5, dpi = 150, bg = "white")
    active_clusters <- obj@misc$pipeline$active_clusters
    graph_name <- paste0(sub("_clusters$", "", active_clusters), "_snn")
    if (!graph_name %in% names(obj@graphs)) abort("Active clustering graph is missing: ", graph_name)
    sweep <- cluster_resolution_summary(obj, graph_name, resolutions, bundle$config$seed)
    resolution <- sweep$summary
    resolution$group <- group; all_resolution[[group]] <- resolution
    tree <- clustree::clustree(sweep$object[[]], prefix = "resolution_") +
      ggraph::scale_edge_colour_gradientn(colours = viridis::viridis(256), guide = "none") +
      ggraph::scale_edge_alpha(guide = "none") +
      ggplot2::labs(title = paste(group, "cluster tree"),
        subtitle = "Each node is a cluster; edges show how cells split as resolution increases from 0.1 to 2.0.") +
      ggplot2::theme(legend.position = "bottom")
    ggplot2::ggsave(file.path(plot_dir, paste0(prefix, "-clustree.png")), tree,
      width = 16, height = 10, dpi = 150, bg = "white", limitsize = FALSE)
    figures <- c(figures, paste0('<h2>', html_escape(group), '</h2><img src="downstream-diagnostics/', prefix, '-pc-variance.png" alt="PC variance"><img src="downstream-diagnostics/', prefix, '-pc-cumulative.png" alt="Cumulative PC variance"><img src="downstream-diagnostics/', prefix, '-clustree.png" alt="Clustering tree across resolutions">'))
  }
  resolution <- do.call(rbind, all_resolution)
  resolution$adjacent_ari <- ifelse(is.na(resolution$adjacent_ari), "—", format(round(resolution$adjacent_ari, 3), nsmall = 3))
  names(resolution) <- c("Resolution", "Clusters", "Smallest cluster", "Median cluster size", "Clusters under 20 cells", "Agreement with prior resolution (ARI)", "Analysis group")
  utils::write.csv(resolution, file.path(run_dir, "resolution_sweep.csv"), row.names = FALSE)
  header <- '<!doctype html><html><head><meta charset="utf-8"><title>Downstream analysis diagnostics</title><style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#192b3b}img{max-width:100%;border:1px solid #d7e0e8;margin:16px 0}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px;border-bottom:1px solid #ddd;text-align:left}.note{background:#eff5fa;padding:16px}.warn{background:#fff4dc;padding:16px}h2{margin-top:42px}</style></head><body>'
  html <- paste0(header, '<h1>PC and clustering-resolution review</h1><div class="note">Current settings remain <strong>', settings$npcs, ' PCs</strong> and <strong>resolution ', settings$resolution, '</strong>. This report does not change them automatically.</div><h2>Workflow order</h2><ol><li>Build PCA from normalized expression and review its variance and loadings.</li><li>When multiple sequencing runs are present, run Harmony on that PCA representation.</li><li>Build the neighbor graph from the active representation: Harmony when correction ran, otherwise PCA.</li><li>Sweep clustering resolution on that active graph and choose the lowest setting that preserves marker-supported populations.</li></ol><div class="warn">ARI measures how similar neighboring resolution clusterings are; it does not prove biological correctness. A small cluster is not automatically bad, especially for a reproducible rare population.</div><h2>Resolution sweep</h2><p>This compares every resolution from 0.1 through 2.0 on the active clustering graph for each analysis group. Download the table as <a href="resolution_sweep.csv">CSV</a>.</p>', html_table(resolution), paste(figures, collapse = ""), '<h2>Evidence</h2><p>Seurat recommends an elbow plot as a fast dimensionality heuristic and encourages checking more than one PC setting. Its clustering tutorial describes resolution as the granularity control: larger values produce more clusters. [<a href="https://satijalab.org/seurat/articles/pbmc3k_tutorial">Seurat guided clustering tutorial</a>; <a href="https://satijalab.org/seurat/reference/elbowplot">ElbowPlot reference</a>]</p><p>The broader workflow should remain data- and biology-aware, as outlined by <a href="https://doi.org/10.15252/msb.20188746">Luecken &amp; Theis (2019)</a>.</p></body></html>')
  path <- file.path(run_dir, "report-downstream-diagnostics.html"); writeLines(html, path); invisible(path)
}
