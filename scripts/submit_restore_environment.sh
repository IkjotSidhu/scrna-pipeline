#!/usr/bin/env bash
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 2
#BSUB -W 1:00
#BSUB -R "rusage[mem=8000]"
#BSUB -R "span[hosts=1]"
#BSUB -J scrna-renv-restore
#BSUB -oo logs/renv-restore-%J.out
#BSUB -eo logs/renv-restore-%J.err
set -euo pipefail

module purge
module load R/4.4.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export RENV_PATHS_PREFIX=linux-rocky-9

pipeline_dir="/sc/arion/projects/naiklab/ikjot/Development/scrna-pipeline"
cd "$pipeline_dir"
R_LIBS=/sc/arion/projects/naiklab/ikjot/R \
  Rscript --vanilla scripts/restore_environment.R "$pipeline_dir"
Rscript scripts/validate_environment.R
Rscript tests/run_tests.R
