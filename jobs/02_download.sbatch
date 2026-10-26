#!/bin/bash
#SBATCH -J download_multiome
#SBATCH -o logs/%x.%A_%a.o
#SBATCH -e logs/%x.%A_%a.e
#SBATCH -p compute
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH -t 2:0:0
#SBATCH --array=0-1

# Downloads the 10x filtered matrix h5 and ATAC fragments for one dataset.
# Array index: 0 = pbmc, 1 = brain.
set -euo pipefail

PROJECT=$SLURM_SUBMIT_DIR  # submit from the repo root
DATASETS=(pbmc brain)
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}

case $DATASET in
    pbmc)
        # PBMC from a Healthy Donor - No Cell Sorting (10k), Cell Ranger ARC 2.0.0
        SAMPLE=pbmc_unsorted_10k
        BASE=https://cf.10xgenomics.com/samples/cell-arc/2.0.0/$SAMPLE ;;
    brain)
        # Frozen Human Healthy Brain Tissue (3k), Cell Ranger ARC 2.0.0
        SAMPLE=human_brain_3k
        BASE=https://cf.10xgenomics.com/samples/cell-arc/2.0.0/$SAMPLE ;;
esac

OUTDIR=$PROJECT/data/raw/$DATASET
mkdir -p "$OUTDIR"

echo "Downloading $DATASET ($SAMPLE) to $OUTDIR"
wget -c -q -P "$OUTDIR" \
    "$BASE/${SAMPLE}_filtered_feature_bc_matrix.h5" \
    "$BASE/${SAMPLE}_atac_fragments.tsv.gz" \
    "$BASE/${SAMPLE}_atac_fragments.tsv.gz.tbi"

ls -lh "$OUTDIR"
gzip -t "$OUTDIR/${SAMPLE}_atac_fragments.tsv.gz" && echo "fragments gzip OK"

echo "Barcode suffixes in fragments:"
# head closes the pipe early; the resulting SIGPIPE is expected.
(zcat "$OUTDIR/${SAMPLE}_atac_fragments.tsv.gz" | grep -v '^#' | head -n 2000000 \
    | awk '{print $4}' | sed 's/.*-//' | sort | uniq -c) || true
