#!/bin/bash
#SBATCH -J azimuth
#SBATCH -o logs/%x.%A_%a.o
#SBATCH -e logs/%x.%A_%a.e
#SBATCH -p compute
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=160G
#SBATCH -t 12:0:0
#SBATCH --array=3

# Azimuth cell-type annotation of a notebook's final cells, outside the notebook (scripts/run_azimuth.py).
# Reads data/scratch/processed/<dataset>/clustered.h5ad, written by the notebook, and writes azimuth.h5ad next to it.
# Array index: 0 = pbmc, 1 = brain, 2 = bmmc, 3 = bmmc_all.
set -euo pipefail

PROJECT=$SLURM_SUBMIT_DIR  # submit from the repo root
DATASETS=(pbmc brain bmmc bmmc_all)
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}
case $DATASET in
    brain) REFERENCE=human_motor_cortex ;;
    *) REFERENCE=human_pbmc ;;
esac
WORK=$PROJECT/data/scratch/processed/$DATASET

module load miniforge
mamba activate azimuth
export PYTHONNOUSERSITE=1
export R_HOME=$CONDA_PREFIX/lib/R
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

python "$PROJECT/scripts/run_azimuth.py" --input "$WORK/clustered.h5ad" \
    --reference "$PROJECT/data/reference/azimuth/$REFERENCE" --out "$WORK/azimuth.h5ad"
