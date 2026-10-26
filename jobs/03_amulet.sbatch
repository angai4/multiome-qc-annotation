#!/bin/bash
#SBATCH -J amulet
#SBATCH -o logs/%x.%A_%a.o
#SBATCH -e logs/%x.%A_%a.e
#SBATCH -p compute
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH -t 12:0:0
#SBATCH --array=2

# AMULET multiplet detection on the Cell Ranger-called cells of one dataset.
# Overlaps are counted once, then multiplets are called per GEM well (barcode suffix).
# Array index: 0 = pbmc, 1 = brain, 2 = bmmc, 3 = bmmc_all.
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
OUT=$PROJECT/data/processed/$DATASET/amulet
AMULET=$PROJECT/tools/AMULET
mkdir -p "$OUT"

module load miniforge
mamba activate multiome
export PYTHONNOUSERSITE=1

python "$PROJECT/scripts/cell_barcodes.py" --h5 "$RAW/${SAMPLE}_filtered_feature_bc_matrix.h5" --out "$OUT/singlecell.csv"

python "$AMULET/FragmentFileOverlapCounter.py" \
    "$RAW/${SAMPLE}_atac_fragments.tsv.gz" "$OUT/singlecell.csv" "$AMULET/human_autosomes.txt" "$OUT"

for B in $(tail -n +2 "$OUT/singlecell.csv" | sed -E 's/.*-([0-9]+),.*/\1/' | sort -u); do
    BDIR=$OUT/batch_$B
    mkdir -p "$BDIR"
    awk -F'\t' -v s="-$B" 'NR==1 || substr($4, length($4)-length(s)+1) == s' \
        "$OUT/Overlaps.txt" > "$BDIR/Overlaps.txt"
    awk -F'\t' -v s="-$B" 'NR==1 || substr($1, length($1)-length(s)+1) == s' \
        "$OUT/OverlapSummary.txt" > "$BDIR/OverlapSummary.txt"
    python "$AMULET/AMULET.py" \
        --rfilter "$AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed" \
        "$BDIR/Overlaps.txt" "$BDIR/OverlapSummary.txt" "$BDIR"
    echo "batch $B:"; cat "$BDIR/MultipletSummary.txt"
done
