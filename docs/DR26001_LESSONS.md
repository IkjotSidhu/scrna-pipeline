# DR26001 pilot: results and lessons

Run date: 2026-10-07  
LSF job: `274070519`  
Run ID: `run-9fa8c41ccb14`  
Status: mandatory ambient correction and doublet removal completed; saved objects independently validated.

## Dataset structure

DR26001 is mouse gut 10x Universal 3' on-chip multiplexing processed with Cell Ranger 9.0.1 and reference `mm10-2020-A`.

| Biological sample label | OCM capture groups | Called cells |
|---|---|---:|
| Mouseileum | OB1, OB2 | 6,434 |
| MouseDuo | OB3, OB4 | 9,541 |

The OCM groups are technical replicates. Donor identity and condition remain unconfirmed, so the pipeline records placeholders and does not treat the four OCM groups as four biological replicates. Ileum and duodenum are separate analysis groups.

## Mandatory preprocessing results

| Capture | Input | SoupX rho | scDblFinder doublets | Retained |
|---|---:|---:|---:|---:|
| OB1 | 3,023 | 1.4% | 132 | 2,558 |
| OB2 | 3,411 | 1.4% | 168 | 2,767 |
| OB3 | 4,294 | 1.0% | 252 | 3,108 |
| OB4 | 5,247 | 1.0% | 315 | 4,112 |

Total retained: 12,545 of 15,975 called cells. SoupX ran on each Cell Ranger biological-sample output because its raw droplet matrix contains both technical OCM groups. scDblFinder used barcode-level `capture_id` and thresholded OB1–OB4 separately with `dbr.per1k = 0.008`, avoiding the incorrect default shipped in early scDblFinder 1.20 releases.

Pilot QC settings were `min_counts=500`, `min_features=200`, `high_features_review=10000`, and `max_mt_percent=25` on corrected counts. They are pilot choices, not validated universal gut thresholds. The high-gene threshold is now a review flag, not an exclusion rule; scDblFinder is used for mandatory doublet removal. The input distributions contain a distinct high-mitochondrial tail, especially in duodenum; this accounts for many exclusions and requires biological review.

## Independent validation

The validation script loaded the saved Seurat objects and confirmed:

- retained cells exactly match the cell ledgers;
- every retained cell is an scDblFinder singlet and completed ambient correction;
- original retained-cell counts remain finite, nonnegative integers;
- corrected counts are finite, nonnegative, aligned with observed counts, and never exceed them;
- feature IDs and symbols survive import and Seurat processing;
- PCA, UMAP, clusters, and required-stage provenance are present.

Final exploratory objects contain 5,325 ileum cells in 16 clusters and 7,220 duodenum cells in 20 clusters. Technical replicates are broadly intermixed in the unintegrated UMAPs. This is a diagnostic observation, not a formal batch-effect metric.

## Computational profile

The full run used R 4.4.1, Seurat 5.5.0, SoupX 1.6.2, and scDblFinder 1.20.2. It completed in 9 minutes 54 seconds with 7.4 GB peak memory. The pilot requested 2 cores and 32 GB total; future jobs of this size can request less memory after a safety margin is selected.

## Development lessons

1. A Cell Ranger sample directory is not necessarily a single capture. OCM barcode groups must remain available as barcode-level metadata.
2. Ambient correction and doublet detection can require different grouping units. Here, SoupX used the per-sample raw matrix, while scDblFinder split by OCM capture group.
3. Input QC should precede filtering. The pilot showed low median mitochondrial fractions alongside a large high-mitochondrial tail that a median-only summary would conceal.
4. Reports must display capture-level distributions. The first report grouped by input matrix and hid OB1–OB4 differences; the capture-aware report fixes this.
5. Successful execution is insufficient validation. A separate script verified the contents and relationships of the saved objects.
6. Small fixtures are valuable. A two-gene Matrix Market fixture found a connection-cleanup bug before it reached an LSF job.

## Remaining limits

- Cell-type labels are unassigned; marker analysis was intentionally disabled for this preprocessing pilot.
- Donor and condition metadata are unresolved.
- The filtered-only DecontX path is not yet runnable in the Minerva project environment because DecontX is missing.
- The baseline `renv` lockfile has not passed a clean restore and is not the final reproducible environment.
- SoupX's global contamination estimate is constant within each biological-sample output. Residual tissue-specific ambient markers still need review during annotation.
