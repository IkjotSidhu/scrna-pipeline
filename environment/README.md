# Minerva environment baseline

Audited on 2026-10-06 under module `R/4.4.1`, platform `x86_64-pc-linux-gnu`.

The existing user library `/sc/arion/projects/naiklab/ikjot/R` was inspected, not modified. The audit placed it first in `.libPaths()`. Minerva's default order places two system libraries before it, so explicit library precedence is necessary during development.

The project pins `RENV_PATHS_PREFIX=linux-rocky-9` on Rocky Linux in `.Rprofile` and the LSF environment scripts. Minerva login and compute nodes currently use different Rocky Linux 9 minor releases; without this pin, `renv` 1.3 creates separate libraries for each minor release. The library still remains separated by R version and CPU platform, and non-Rocky development machines retain their native prefix.

## Findings

| Package | Version | Loads | Source library |
|---|---|---|---|
| Seurat | 5.5.0 | Yes | User |
| SeuratObject | 5.4.0 | Yes | User |
| Matrix | 1.7-5 | Yes | User |
| SoupX | 1.6.2 | Yes | System R 4.4.1 site library |
| scDblFinder | 1.20.2 | Yes | User |
| SingleCellExperiment | 1.28.1 | Yes | User |
| BiocParallel | 1.40.2 | Yes | User |
| renv | 1.3.0 | Yes | User |
| decontX | Missing | No | — |

Other checked dependencies loaded without captured warnings. Detailed versions, paths, build information and the package inventory are in `initial-audit/`.

## Reproducibility status

`renv.baseline.lock` records the starting user environment and intentionally excludes the then-missing DecontX. The root `renv.lock` is the supported pipeline environment. It restored into a second empty library on Minerva and loaded every direct dependency, including DecontX 1.4.1.

Completed environment checks:

1. DecontX was added to a private project environment without modifying the user's library.
2. The dependency graph and final `renv.lock` were created for R 4.4.1/Bioconductor 3.20.
3. The lock restored into `environment/restore-test/library` and all direct dependencies loaded.
4. Synthetic import fixtures and the DR26001 SoupX/scDblFinder pilot passed on compute nodes.
5. Each run records package versions, code and input hashes, seeds, effective settings, and thread limits.

Still required: test restoration on a second Minerva date/node after cache changes, validate the DecontX algorithm path on real filtered-only data, and add system-library inventory checks for packages linked to HDF5/BLAS.

A lockfile records R package versions and sources. It does not by itself freeze the OS, external libraries, datasets or code, and it does not guarantee bit-identical floating-point results across platforms.
