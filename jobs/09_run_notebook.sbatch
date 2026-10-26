#!/bin/bash
#SBATCH -J run_notebook
#SBATCH -o logs/%x.%A_%a.o
#SBATCH -e logs/%x.%A_%a.e
#SBATCH -p compute
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 6:0:0
#SBATCH --array=0-2

# Runs a v3 notebook top to bottom, in place, with the "Python (azimuth)" kernel. Rewrites its outputs, its
# data/processed files and figures/<notebook>/*.svg.
# Array index: 0 = pbmc, 1 = brain, 2 = bmmc, 3 = bmmc_all (longer walltime; see jobs/run_notebooks.sh).
set -euo pipefail

PROJECT=$SLURM_SUBMIT_DIR  # submit from the repo root
NOTEBOOKS=(pbmc_v3 brain_v3 bmmc_s1d1_s1d3_v3 bmmc_all_v3)
NOTEBOOK=${NOTEBOOKS[$SLURM_ARRAY_TASK_ID]}

module load miniforge
mamba activate multiome
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=azimuth \
    --ExecutePreprocessor.timeout=-1 "$PROJECT/notebooks/$NOTEBOOK.ipynb"
echo "$(ls "$PROJECT/figures/$NOTEBOOK" | wc -l) SVGs in figures/$NOTEBOOK"
