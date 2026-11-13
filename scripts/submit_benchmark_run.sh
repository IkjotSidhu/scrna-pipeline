#!/usr/bin/env bash
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 1
#BSUB -W 4:00
#BSUB -R "rusage[mem=12000]"
#BSUB -R "span[hosts=1]"
#BSUB -J scrna-benchmark-run
#BSUB -oo logs/benchmark-run-%J.out
#BSUB -eo logs/benchmark-run-%J.err
set -euo pipefail

module purge
module load R/4.4.1
export RENV_PATHS_PREFIX=linux-rocky-9
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1

pipeline_dir="/sc/arion/projects/naiklab/ikjot/Development/scrna-pipeline"
cd "$pipeline_dir"

# Unlike the qc smoke test, run enforces ambient-RNA correction and doublet
# removal before clustering and report generation.
Rscript run.R run benchmarks/config/project.yml
