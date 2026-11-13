# Design decisions and ideas under review

## Analysis groups and batch correction

`analysis_group` defines which compatible inputs are merged for one final biological analysis. An analysis group may contain multiple tissues, samples, libraries, and physical captures, provided species and reference annotation agree. Tissue remains cell-level metadata, so ileum and duodenum can be analyzed jointly and compared explicitly.

Harmony batch correction runs automatically when an analysis group contains more than one sequencing run, represented by `batch_id`. The uncorrected PCA, graph, clusters, and UMAP always remain available for audit. Harmony is the default method and RPCA is an alternative. `batch_variable` selects `batch_id`, `library_id`, `capture_id`, or `sample_id`. DR26001 joins ileum and duodenum in one analysis group and has one sequencing run, so it remains unintegrated.

## Mitochondrial expression should be interpreted in context

A universal mitochondrial cutoff can remove valid metabolically active cells, including mature intestinal enterocytes. The current cutoff remains a configurable development filter, but the intended design is a two-pass review:

1. Use permissive initial flags and retain a cell-level reason ledger.
2. Assign broad cell classes using marker and reference evidence.
3. Compare mitochondrial percentage within cell type, tissue, and capture against healthy reference distributions.
4. Review mitochondrial percentage together with UMI/gene complexity, stress and apoptosis programs, ambient contamination, doublet score, and the shape of each capture's distribution.
5. Distinguish a biological high-mito state from probable damage; record the final decision and sensitivity analysis rather than silently changing a threshold.

This feature is a design requirement and is not implemented yet.

## Two audiences, two reports

Every completed run produces:

- `report-bioinformatics.html`: capture-level QC, ambient correction, doublets, exclusion reasons, technical embeddings, configuration, versions, and provenance.
- `report-biologist.html`: samples, retained cells, tissue/condition/cluster structure, annotation status, and focused biological review questions.

`report.html` is a small index linking to both. The reports read saved results and can be rebuilt without rerunning preprocessing.
