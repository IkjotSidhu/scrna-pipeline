# Evidence basis for QC thresholds — draft for review

## What the literature supports

Review and best-practice articles consistently recommend inspecting library size / UMI counts, detected genes, mitochondrial proportion, and doublet evidence. They do **not** establish one set of numerical cutoffs that applies across tissues, chemistries, sample handling, sequencing depths, or species.

The recurring recommendation is to assess these metrics within a library or batch, use robust distribution-based outlier detection, and document the final decision. A threshold that is suitable for PBMCs can remove valid epithelial or tissue-resident populations in skin or gut.

## Sources to use for pipeline decisions

| Source | Type | Relevant principle for this pipeline |
|---|---|---|
| Luecken MD, Theis FJ. *Current best practices in single-cell RNA-seq analysis: a tutorial.* Molecular Systems Biology (2019). doi:10.15252/msb.20188746 | Review/tutorial | QC is based on counts, detected genes, mitochondrial proportion, and context-specific inspection rather than a universal cutoff. |
| Heumos L, et al. *Best practices for single-cell analysis across modalities.* Nature Reviews Genetics (2023). doi:10.1038/s41576-023-00586-w | Review | QC should account for platform, modality, and experimental design; decisions should be recorded and reproducible. |
| Amezquita RA, et al. *Orchestrating single-cell analysis with Bioconductor.* Nature Methods (2020). doi:10.1038/s41592-019-0654-x | Best-practice workflow | Library-level QC and robust outlier handling should be performed before downstream analysis. |
| Lun ATL, McCarthy DJ, Marioni JC. *A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.* F1000Research (2016). doi:10.12688/f1000research.9501.2 | Methods workflow | Use library-specific distributions and robust statistics such as median absolute deviation (MAD) to flag low library size, low feature count, and high mitochondrial proportion. |
| Ilicic T, et al. *Classification of low quality cells from single-cell RNA-seq data.* Genome Biology (2016). doi:10.1186/s13059-016-0888-1 | Methods study | Low quality is multivariate; no single metric should be interpreted alone. |

## Consequence for our pipeline

The current values (`500` counts, `200` genes, `25%` mitochondrial reads) are **provisional guardrails**, not literature-established universal standards. The report must distinguish:

1. Configured hard rules used for a run.
2. Library-specific distributional flags.
3. Final removals and their reason.
4. Cell types or tissues where high mitochondrial proportion may be biologically plausible.

## Decision still needed

Before replacing the provisional defaults, decide whether the first implementation should:

- retain conservative hard guardrails and add MAD-based review flags; or
- use MAD-based thresholds as the primary filter, with only a minimal absolute floor for near-empty droplets.

Either implementation should preserve every cell's QC values and exclusion reason in the cell ledger.
