#!/usr/bin/env bash
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 2
#BSUB -W 2:00
#BSUB -R "rusage[mem=16000]"
#BSUB -R "span[hosts=1]"
#BSUB -J scrna-environment
#BSUB -oo logs/environment-%J.out
#BSUB -eo logs/environment-%J.err
set -euo pipefail

module purge
module load R/4.4.1
export R_LIBS="/sc/arion/projects/naiklab/ikjot/R:/hpc/packages/minerva-rocky9/rpackages/4.4.1/site-library:/hpc/packages/minerva-rocky9/rpackages/bioconductor/3.20"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export RENV_PATHS_PREFIX=linux-rocky-9

pipeline_dir="/sc/arion/projects/naiklab/ikjot/Development/scrna-pipeline"
cd "$pipeline_dir"
Rscript --vanilla scripts/create_environment.R "$pipeline_dir"
