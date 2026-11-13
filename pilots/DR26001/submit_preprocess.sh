#!/usr/bin/env bash
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 2
#BSUB -W 2:00
#BSUB -R "rusage[mem=16000]"
#BSUB -R "span[hosts=1]"
#BSUB -J scrna-DR26001-process
#BSUB -oo logs/DR26001-process-%J.out
#BSUB -eo logs/DR26001-process-%J.err
set -euo pipefail

cd /sc/arion/projects/naiklab/ikjot/Development/scrna-pipeline
source scripts/minerva_environment.sh
Rscript run.R run pilots/DR26001/config/project.yml
