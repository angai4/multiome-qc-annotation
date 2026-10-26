#!/bin/bash
#SBATCH -J atac_metrics
#SBATCH -o logs/%x.%A_%a.o
#SBATCH -e logs/%x.%A_%a.e
#SBATCH -p compute
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH -t 2:0:0
#SBATCH --array=2

# Per-cell ATAC QC metrics (scripts/atac_metrics.py) for the Cell Ranger-called cells of one dataset.
# Writes data/processed/<dataset>/atac_metrics.tsv. Array index: 0 = pbmc, 1 = brain, 2 = bmmc, 3 = bmmc_all.
set -euo pipefail

PROJECT=$SLURM_SUBMIT_DIR  # submit from the repo root
DATASETS=(pbmc brain bmmc bmmc_all)
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}
case $DATASET in
    pbmc) SAMPLE=pbmc_unsorted_10k ;;
    brain) SAMPLE=human_brain_3k ;;
    bmmc) SAMPLE=bmmc_s1d1_s1d3 ;;  # NeurIPS 2021 s1d1 + s1d3, aggregated by jobs/07_bmmc_aggr.sbatch
    bmmc_all) SAMPLE=bmmc_all ;;  # all 13 NeurIPS 2021 multiome samples (jobs/07_bmmc_aggr.sbatch)
esac
RAW=$PROJECT/data/raw/$DATASET
OUT=$PROJECT/data/processed/$DATASET
mkdir -p "$OUT"

module load miniforge
mamba activate multiome
export PYTHONNOUSERSITE=1

python "$PROJECT/scripts/atac_metrics.py" --h5 "$RAW/${SAMPLE}_filtered_feature_bc_matrix.h5" \
    --fragments "$RAW/${SAMPLE}_atac_fragments.tsv.gz" --out "$OUT/atac_metrics.tsv" --n-jobs "$SLURM_CPUS_PER_TASK"
