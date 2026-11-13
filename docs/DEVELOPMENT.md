# Development and learning notes

## Current status

The modular draft exists, import fixtures pass, and DR26001 completed the SoupX/scDblFinder preprocessing path on Minerva. Saved objects passed independent consistency checks. Optional integration and biological annotation have not been validated.

Local development uses R 4.5.2/Seurat 5.4.0. The supported Minerva runtime is R 4.4.1. Its project-local `renv` library was built on Linux and does not reuse compiled macOS packages or modify the existing lab R library.

The Minerva `renv.lock` restored into a second isolated library and loaded all required packages, including Seurat 5.5.0, SoupX 1.6.2, scDblFinder 1.20.2, and DecontX 1.4.1. `.Rprofile` pins a common prefix only on Rocky Linux so login and compute nodes activate the same project library despite different OS minor releases. The macOS development library is not part of the supported runtime.

## What the architecture teaches

1. **Separate biology from technology.** `donor_id` and `sample_id` identify biological units; `library_id` identifies an input matrix, and barcode-level `capture_id` identifies the droplet loading partition. Ambient contamination and droplet doublets are capture-level problems. Biological replication is donor/sample-level.
2. **Configuration is data.** YAML is read once and validated. Unknown settings are errors. This avoids notebooks overwriting values and sourced R configuration executing unexpected code.
3. **Preserve original measurements.** Corrected counts are model estimates. Keeping observed counts makes corrections auditable and allows later sensitivity analyses.
4. **Required stages are contracts.** A failed correction or doublet step stops a full run. Successful software execution does not guarantee perfect contaminant/doublet detection.
5. **Reports consume results.** Rebuilding a plot should not rerun correction. The report command reads saved summaries.
6. **Cache dependencies, not file existence.** Draft cache keys include inputs, settings, code and package versions. The implementation is a lightweight RDS cache, not `targets`; fine-grained invalidation remains to be tested and refined.
7. **Test in layers.** Parsing catches syntax errors; synthetic fixtures test barcode/feature alignment and exclusion behavior; a real pilot tests practical and biological behavior.

## Next validation steps

- Confirm donor and condition metadata for DR26001.
- Build fixtures covering MTX/H5 parity, repeated barcodes across captures, duplicate gene symbols, missing raw droplets, species mismatch and forbidden stage disabling.
- Verify feature metadata survive Seurat merge/layer operations.
- Verify cache invalidation, resume behavior and failed-stage reporting.
- Run the filtered-only DecontX path on a suitable pilot dataset; DR26001 validated SoupX plus scDblFinder.
- Review contamination estimates, raw/corrected marker profiles, doublet scores and cell-retention counts before treating results as usable.
- Test optional RPCA independently, including small batches and confounded designs.
- Add automated end-to-end fixtures and document supported resource profiles across study sizes.

## Current scope limits

10x v3 H5 and legacy/current MTX layouts are drafted. Only Gene Expression is analyzed; excluded modality types are recorded, not analyzed. Multiplexed captures need an explicit barcode-to-sample mapping extension and are not supported yet. Large-scale disk-backed execution, pseudobulk/composition inference, reference annotation and Biomni are later milestones. The current draft holds multiple processed libraries in memory and has not been benchmarked for large studies.

## References

- Seurat: https://satijalab.org/seurat/
- SoupX: https://github.com/constantAmateur/SoupX
- DecontX: https://bioconductor.org/packages/decontX/
- scDblFinder: https://bioconductor.org/packages/scDblFinder/
