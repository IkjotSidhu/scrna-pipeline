#!/usr/bin/env bash
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 2
#BSUB -W 2:00
#BSUB -R "rusage[mem=16000]"
#BSUB -R "span[hosts=1]"
#BSUB -J scrna-DR26001-qc
#BSUB -oo logs/DR26001-qc-%J.out
#BSUB -eo logs/DR26001-qc-%J.err
set -euo pipefail
cd /sc/arion/projects/naiklab/ikjot/Development/scrna-pipeline
source scripts/minerva_environment.sh
Rscript pilots/DR26001/input_qc.R
