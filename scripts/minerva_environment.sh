# Source this in an interactive Minerva shell or LSF job before pipeline runs.
# R loads the project .Rprofile, which activates the locked renv library.
module purge
module load R/4.4.1 || return 1
unset SCRNA_R_LIB
export RENV_PATHS_PREFIX=linux-rocky-9
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1
