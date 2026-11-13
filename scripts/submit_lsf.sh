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
PIPELINE_DIR="/REPLACE/project/scrna-pipeline"
CONFIG_PATH="/REPLACE/project/scrna-pipeline/config/project.yml"
export SCRNA_PIPELINE_DIR="$PIPELINE_DIR"
# Run mkdir -p logs before bsub. This loads R/4.4.1 and isolates renv from
# inherited user/site R libraries.
source "$PIPELINE_DIR/scripts/minerva_r_env.sh"
cd "$PIPELINE_DIR"
# Start with qc. Switch to run only after inspecting the input report.
Rscript run.R qc "$CONFIG_PATH"
