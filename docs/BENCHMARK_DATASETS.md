# Benchmark dataset shortlist

The benchmark set should test both software behavior and biological plausibility. Processed reference atlases are useful for annotation and cell-type QC, but they cannot validate SoupX unless raw droplets are also available.

## Mouse small intestine

**Primary pipeline benchmark candidate:** 10x Genomics, “40k Mixture of Cells Dissociated from 4 Fixed Mouse Tissues using gentleMACS,” Cell Ranger 7.1.0. It contains 8,916 healthy C57BL/6 small-intestine cells in a four-sample Flex experiment. This is an official, well-documented Cell Ranger dataset, but its fixed-RNA chemistry differs from DR26001's Universal 3' chemistry. Download only the Cell Ranger output and sample-level matrix needed for development before considering FASTQs.

Source: https://www.10xgenomics.com/datasets/40k-mixture-of-cells-dissociated-from-4-fixed-mouse-tissues-using-gentlemacs-dissociator-multiplexed-samples-4-probe-barcodes-1-standard

**Biological reference candidate:** MIST ileum CITE-seq, a large healthy C57BL/6 distal-small-intestine atlas with epithelial and lamina-propria cells. Its size makes it better as an annotation/reference benchmark than a routine end-to-end test.

Source: https://gutpath.org/ileum_cite-seq

## Mouse skin

**Selected public candidate:** GSE142471, “Single cell RNA-Seq adult mouse skin during homeostasis and wound repair.” It contains two independent unwounded adult mouse back-skin controls generated with 10x Chromium 3' v2 and Cell Ranger 2.1.0. The controls provide biological replication and standard MTX/barcode/gene files. Raw droplets are not deposited, so this benchmark exercises import, DecontX fallback, QC, clustering, and annotation rather than SoupX.

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142471

The Tabula Muris droplet release was downloaded and inspected during candidate validation. Contrary to secondary catalog descriptions, Figshare release 5968960 contains no skin matrices and its annotations contain no skin cells. It is recorded as rejected in `benchmarks/datasets.csv` and is not part of the reproducible download manifest.

## Internal Naik lab controls

Past lab runs should be registered in `benchmarks/internal_datasets.csv`, using `internal_datasets.template.csv` as the schema. Comparisons should use the closest stratum by species, tissue, assay chemistry, cell versus nucleus, fresh/frozen/fixed preparation, reference build, and major protocol era. The registry should contain study-safe identifiers and paths, without patient identifiers.

Internal historical distributions are the primary QC comparator because they reflect the lab's handling and sequencing practices. Public datasets provide external biological context. Current-study captures provide the most direct technical peer comparison.

The initial read-only inventory on Minerva found two useful candidates that still need experimental metadata review before they can enter the benchmark registry:

- DR26001 contains current mouse ileum and duodenum sample-level raw and filtered Cell Ranger matrices. It is suitable for testing the complete SoupX path and joint-tissue analysis, but it should be treated as a development dataset until its biological condition and preparation metadata are finalized.
- `Daniel-Redone-Analysis-NYU-Transfer/2024-06-05-Daniel-Single-Cell` contains historical ileum and duodenum Seurat objects, including control-versus-SFB analyses. These processed objects may provide cell-type and mitochondrial-QC context, but they cannot independently benchmark ambient correction without the original raw droplet matrices.

Do not infer control status from a directory or object name. Add a historical dataset only after confirming tissue, condition, preparation, chemistry, reference build, and which processing steps were already applied.

## Acceptance criteria

For each benchmark, record source URL, license, assay chemistry, organism/strain, tissue, health state, Cell Ranger/reference version, file checksums, and whether raw droplets exist. A benchmark passes only when expected broad compartments are recovered, mandatory stages complete, cell/feature identities remain aligned, and the two reports render from saved outputs.

## Execution stages

The public benchmark has two deliberately separate stages:

1. `bsub < scripts/submit_benchmark_qc.sh` validates imports, original-count QC metrics, and report rendering without claiming that processing is complete.
2. `bsub < scripts/submit_benchmark_run.sh` performs the mandatory processing path. The 10x small-intestine sample uses SoupX with raw droplets; the two filtered-only skin controls use DecontX. All captures then undergo scDblFinder doublet detection, corrected-count QC, clustering, and dual-report generation.

The second stage is the benchmark used for biological and technical comparisons. The QC-only output is a software smoke test and must not be treated as a processed reference.

The first complete public run (`run-067b5fe04510`) processed 19,805 called cells in about 17 minutes on Minerva and peaked at 7.3 GB RAM. Because the workload used one CPU effectively, the benchmark submission requests one core and 12 GB for subsequent runs. Larger studies should size resources from their own pilot rather than copying this benchmark request.
