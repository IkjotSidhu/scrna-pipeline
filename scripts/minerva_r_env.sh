#!/usr/bin/env bash
# Source this file before restoring or running the pipeline on Minerva.
# It deliberately prevents package discovery from inherited user/site libraries.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Source this file: source scripts/minerva_r_env.sh" >&2
  exit 2
fi

SCRNA_PIPELINE_DIR="${SCRNA_PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

module purge
module load R/4.4.1

# A module-compatible R package may still be binary-incompatible with another
# R build. Start clean; renv supplies only this project's locked packages.
unset R_LIBS R_LIBS_USER R_LIBS_SITE SCRNA_R_LIB
export RENV_CONFIG_USER_PROFILE=FALSE
export RENV_PATHS_PREFIX="${RENV_PATHS_PREFIX:-linux-rocky-9}"
export RENV_PATHS_CACHE="${SCRNA_RENV_CACHE:-${SCRNA_PIPELINE_DIR}/.renv-cache}"

# Prevent threaded BLAS from oversubscribing an LSF allocation.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLAS_NUM_THREADS=1
