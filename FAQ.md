# FAQ

## Why are there two input modes?

The pipeline begins with Cell Ranger filtered matrices in every case. When matching raw droplets are also available, it uses SoupX; otherwise it uses DecontX. Ambient-RNA correction is mandatory in a full run.

## Are doublets removed?

Yes. scDblFinder is run per physical capture. Predicted doublets are excluded from the final Seurat object, while their score and call remain in the cell ledger.

## How are QC thresholds chosen?

The example settings are starting values, not universal biological rules. The pipeline shows retention percentages, distributions, and cell-level flags so each project can document a justified choice. High detected-gene count is a review flag; it is not a standalone doublet filter.

## When does Harmony run?

Harmony runs automatically when an `analysis_group` includes more than one sequencing run (`batch_id`). The unintegrated PCA and clustering remain in the Seurat object for audit.

## How are PCs and clustering resolution chosen?

Use the PCA variance plots to review the dimensions entering the active representation. Then inspect the clustree sweep from 0.1 through 2.0 and select the lowest resolution that preserves marker-supported populations. Do not choose a resolution from UMAP appearance alone.

## Which report should I read first?

Start with `report-bioinformatics.html` for cell retention and sequencing QC. Use the sample- and sequencing-run visual reports to inspect distributions and outliers. Review the PCA/clustree report before finalizing clustering settings. The biological report becomes useful after annotation review.

## What should be version controlled?

Commit code, configuration templates, manifests without sensitive paths, documentation, and `renv.lock`. Do not commit Cell Ranger matrices, raw droplets, result folders, RDS objects, or sample-identifying metadata.
