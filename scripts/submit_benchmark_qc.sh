#!/usr/bin/env bash
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 2
#BSUB -W 1:00
#BSUB -R "rusage[mem=16000]"
#BSUB -R "span[hosts=1]"
#BSUB -J scrna-benchmark-qc
#BSUB -oo logs/benchmark-qc-%J.out
#BSUB -eo logs/benchmark-qc-%J.err
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
Rscript run.R qc benchmarks/config/project.yml
