#!/usr/bin/env bash
# Build the isolated, lockfile-defined Minerva environment once after cloning
# or whenever renv.lock changes.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SCRNA_PIPELINE_DIR="$PIPELINE_DIR"
source "$PIPELINE_DIR/scripts/minerva_r_env.sh"
cd "$PIPELINE_DIR"
R -e 'renv::restore(prompt = FALSE)'
