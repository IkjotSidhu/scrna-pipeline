# scRNA-seq pipeline

![Illustrated path from droplets to quality-controlled single-cell profiles](assets/scrna-pipeline-banner.png)

A practical Seurat pipeline for 10x Genomics single-cell RNA-seq data. It starts from Cell Ranger count matrices and produces an auditable Seurat object, QC reports, clustering diagnostics, marker tables, and annotation review sheets.

It is designed for mouse or human skin, gut, and PBMC datasets, but uses a general sample manifest so other 10x studies can be added without rewriting the workflow.

## What it does

1. Reads Cell Ranger filtered matrices and raw droplets when available.
2. Removes ambient RNA with SoupX when raw droplets exist, otherwise DecontX.
3. Calls and removes doublets with scDblFinder.
4. Applies transparent cell-level QC rules and saves a ledger for every called cell.
5. Runs Harmony automatically when an analysis group has more than one sequencing run.
6. Generates table-first QC reports, visual QC reports by sample and sequencing run, and PCA/clustree diagnostics.
7. Exports markers and a review sheet for broad and fine cell-type annotation.

## Quick start

```bash
git clone https://github.com/IkjotSidhu/scrna-pipeline.git
cd scrna-pipeline
R -e 'renv::restore()'
```

### Minerva setup

On Minerva, first make the environment from the locked package versions:

```bash
git clone https://github.com/IkjotSidhu/scrna-pipeline.git
cd scrna-pipeline
bash scripts/minerva_restore.sh
```

The setup loads the required `R/4.4.1` module, clears inherited `R_LIBS*`
paths, and uses this clone's `renv` library and cache. This prevents packages
from a non-editable shared library or a user's unrelated R library from being
silently mixed into the run. Do not set `SCRNA_R_LIB` for this pipeline.

1. Copy and edit [config/project.yml](config/project.yml).
2. Copy and edit [config/samples.csv](config/samples.csv), with one row per Cell Ranger matrix.
3. Inspect the input first:

```bash
Rscript run.R qc config/project.yml
```

4. After reviewing the input report, run the full pipeline:

```bash
Rscript run.R run config/project.yml
```

## Everyday commands

```bash
Rscript run.R report /path/to/results/run-<id>
Rscript run.R qc-plots /path/to/results/run-<id>
Rscript run.R diagnostics /path/to/results/run-<id>
```

For Minerva, run the setup once, then copy [scripts/submit_lsf.sh](scripts/submit_lsf.sh), fill in its placeholders, and submit it with `bsub`. The job template applies the same clean environment automatically.

## Input manifest

`config/samples.csv` records sample metadata and paths. Required fields include `library_id`, `sample_id`, `species`, `tissue`, `condition`, `batch_id`, `analysis_group`, and `filtered_matrix`.

- Use `raw_matrix` whenever available. It enables SoupX.
- Use `cell_metadata` with `barcode,capture_id` for multiplexed libraries.
- Keep biological samples in `sample_id`; put sequencing runs in `batch_id`.
- Libraries in the same `analysis_group` are analyzed together only when their species and reference are compatible.

## Results

```text
results/run-<id>/
├── report.html                         # report index
├── report-bioinformatics.html          # QC, retention, Cell Ranger metrics
├── report-qc-by-sample.html            # visual QC by biological sample
├── report-qc-by-sequencing-run.html    # visual QC by sequencing run
├── report-downstream-diagnostics.html  # PCA and clustree review
├── groups/<analysis-group>/seurat.rds  # final Seurat object
├── groups/<analysis-group>/markers.csv
├── groups/<analysis-group>/annotation_review.csv
└── libraries/<library>/cell_ledger.csv # every called cell and its QC decision
```

The final object retains ambient-corrected values in `RNA` and observed retained-cell counts in `RNA_observed`.

## Repository map

| Folder | Contents |
|---|---|
| `config/` | Project and manifest templates |
| `R/` | Pipeline implementation |
| `scripts/` | Generic LSF template and benchmark downloader |
| `benchmarks/` | Public benchmark registry and example configuration |
| `docs/` | Benchmark notes and QC evidence |
| `tests/` | Configuration, import, and report tests |

Read [FAQ.md](FAQ.md) for ambient RNA, doublets, mitochondrial RNA, batch correction, and clustering. Benchmark sources are listed in [docs/BENCHMARK_DATASETS.md](docs/BENCHMARK_DATASETS.md).

## Branches

`dev` is the active development branch. `main` is the reviewed release branch.
