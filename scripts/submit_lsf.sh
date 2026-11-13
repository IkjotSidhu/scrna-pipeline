#!/usr/bin/env bash
# Replace every placeholder before submission. Resource values require confirmation.
#BSUB -P <ACCOUNT>
#BSUB -q <QUEUE>
#BSUB -n <CORES>
#BSUB -W <HH:MM>
#BSUB -R "rusage[mem=<MEM_MB_PER_CORE>]"
#BSUB -R "span[hosts=1]"
#BSUB -J scrna-pilot
#BSUB -oo logs/scrna-%J.out
#BSUB -eo logs/scrna-%J.err
set -euo pipefail
# Run mkdir -p logs before bsub. Confirm the R module with your Minerva environment.
module purge
module load R/4.4.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLAS_NUM_THREADS=1
export SCRNA_R_LIB="/REPLACE/project/R-library"
PIPELINE_DIR="/REPLACE/project/scrna-pipeline"
CONFIG_PATH="/REPLACE/project/scrna-pipeline/config/project.yml"
cd "$PIPELINE_DIR"
# Start with qc. Switch to run only after inspecting the input report.
Rscript run.R qc "$CONFIG_PATH"
