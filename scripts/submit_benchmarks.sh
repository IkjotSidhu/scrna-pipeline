#!/usr/bin/env bash
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 1
#BSUB -W 1:00
#BSUB -R "rusage[mem=4000]"
#BSUB -R "span[hosts=1]"
#BSUB -J scrna-benchmarks
#BSUB -oo logs/benchmarks-%J.out
#BSUB -eo logs/benchmarks-%J.err
set -euo pipefail

module purge
module load R/4.4.1
pipeline_dir="/sc/arion/projects/naiklab/ikjot/Development/scrna-pipeline"
cd "$pipeline_dir"
Rscript --vanilla scripts/download_benchmarks.R "$pipeline_dir"
