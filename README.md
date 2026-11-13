# Single-cell RNA-seq pipeline

Development draft with a validated preprocessing pilot and clean environment restore. Cell-type annotation and downstream biological inference remain under development.

The pipeline supports local execution and Minerva/LSF. Configure input and output paths in your project YAML and manifest; no site-specific paths are required by the core workflow.

## Agreed requirements

- Seurat-based R analysis of human or mouse 10x count matrices.
- Skin, gut and PBMC studies; separate analysis groups for incompatible tissues/species/references.
- **Ambient RNA correction and doublet removal are mandatory for every full run.**
- SoupX with matched raw droplets; DecontX when only filtered counts are available. Methods and their diagnostics must be recorded per capture.
- scDblFinder predicts doublets per physical capture; predicted doublets are excluded from final analysis objects.
- Input-only QC inspection is allowed but is explicitly not a completed preprocessing run.
- Compatible tissues and samples can share an `analysis_group` for one joint final object.
- Harmony batch correction runs automatically when an analysis group includes multiple sequencing runs (`batch_id`); the unintegrated baseline is retained for audit. RPCA remains an alternative method.
- Cell-type annotation will combine marker/reference evidence and optional Biomni assistance. The current draft only exports a review sheet; it does not call Biomni.
- Minerva/LSF and local execution. Do not run heavy analysis on login nodes.

## Files and responsibilities

| Path | Purpose |
|---|---|
| `run.R` | Command-line entry point |
| `config/project.yml` | One validated configuration; illustrative QC defaults |
| `config/samples.csv` | One row per input matrix; optional barcode metadata identifies multiple captures within it |
| `R/config.R` | Validate configuration, metadata and mandatory-stage policy |
| `R/import.R` | Read 10x matrices and preserve gene IDs/symbols |
| `R/process.R` | Ambient correction, doublets, filtering and Seurat analysis |
| `R/engine.R` | Execute stages, fingerprint inputs and record provenance |
| `R/report.R` | Render reports from saved summaries |
| `R/qc_visual_report.R` | Render sample- and sequencing-run QC plot reports |
| `R/downstream_diagnostics.R` | Render PCA and clustree diagnostics across resolutions 0.1–2.0 |
| `scripts/submit_lsf.sh` | LSF template; account/resources require confirmation |
| `scripts/setup.R` | Explicit dependency installation into a private library |
| `docs/DEVELOPMENT.md` | Validation status, limitations and next steps |
| `docs/DESIGN_DECISIONS.md` | Joint analysis, mitochondrial-QC, and report design |
| `docs/BENCHMARK_DATASETS.md` | Curated normal-tissue benchmark candidates and provenance |

## Intended commands (after validation and environment setup)

```bash
Rscript run.R qc config/project.yml
Rscript run.R run config/project.yml
Rscript run.R report /absolute/path/to/completed-run
Rscript run.R qc-plots /absolute/path/to/completed-run
Rscript run.R diagnostics /absolute/path/to/completed-run
```

Configuration paths are resolved relative to the YAML file; matrix paths are resolved relative to the sample manifest. Analysis never installs dependencies automatically.

`ambient.method: auto` selects a method according to raw-matrix availability; it never skips correction. `none` is not valid for ambient correction or doublet detection. Disabling doublet removal is not supported.

For multiplexed data, `cell_metadata` points to a CSV containing at least `barcode,capture_id`. It must cover every filtered barcode exactly once. scDblFinder splits and thresholds by `capture_id`; technical capture groups remain distinct from biological `sample_id` and `donor_id`.

Thresholds in the example configuration are starting values to review, not universal recommendations for skin, gut or PBMC. Counts in the `RNA` assay are ambient-adjusted estimates; original retained-cell UMI counts are kept in `RNA_observed`. Doublet scoring currently uses observed integer counts. Final QC uses corrected counts. These choices must be checked on the pilot.

DR26001 has been processed as a development pilot; see `docs/DR26001_LESSONS.md`. This validates the implemented preprocessing path on that dataset, not every supported input type or biological use case.

## Git workflow

`main` contains reviewed, working releases. Development continues on `dev`; changes should be tested there before being merged into `main`.
