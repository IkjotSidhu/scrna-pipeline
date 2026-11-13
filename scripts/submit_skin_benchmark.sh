#!/usr/bin/env bash
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 1
#BSUB -W 1:00
#BSUB -R "rusage[mem=4000]"
#BSUB -R "span[hosts=1]"
#BSUB -J scrna-skin-benchmark
#BSUB -oo logs/skin-benchmark-%J.out
#BSUB -eo logs/skin-benchmark-%J.err
set -euo pipefail

pipeline_dir="/sc/arion/projects/naiklab/ikjot/Development/scrna-pipeline"
destination="$pipeline_dir/benchmarks/data/GSE142471"
archive="$destination/GSE142471_RAW.tar"
mkdir -p "$destination/raw"
if [[ ! -s "$archive" ]]; then
  curl --fail --location --retry 3 \
    --output "$archive" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142471/suppl/GSE142471_RAW.tar"
fi
md5sum "$archive" > "$archive.md5"
tar -xf "$archive" -C "$destination/raw"
printf 'Downloaded and extracted GSE142471\n'
