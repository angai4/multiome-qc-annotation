#!/bin/bash
#SBATCH -J bmmc_aggr
#SBATCH -o logs/%x.o%j
#SBATCH -e logs/%x.e%j
#SBATCH -p compute
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH -t 12:0:0

# NeurIPS 2021 BMMC multiome: downloads each sample's ATAC fragments (S3) next to its rebuilt GEX matrix in
# data/raw/bmmc_all/<sample>/ (jobs/06_rebuild_gex.sbatch), checks their Cell Ranger header, then combines the samples
# aggr-style (suffix -1, -2, ... in table order) into data/raw/$DATASET/. Reads DATASET: bmmc (s1d1 + s1d3) or
# bmmc_all (default; every row of metadata/bmmc_samples.tsv).
set -euo pipefail

PROJECT=$SLURM_SUBMIT_DIR  # submit from the repo root
cd "$PROJECT"
module load miniforge
mamba activate multiome

S3=https://openproblems-bio.s3.amazonaws.com/public/post_competition/multiome
DATASET=${DATASET:-bmmc_all}
case $DATASET in
    bmmc) SAMPLES=(s1d1 s1d3); PREFIX=bmmc_s1d1_s1d3 ;;
    bmmc_all) SAMPLES=($(awk -F'\t' 'NR > 1 {print $1}' metadata/bmmc_samples.tsv)); PREFIX=bmmc_all ;;
    *) echo "unknown DATASET $DATASET" >&2; exit 1 ;;
esac
echo "$DATASET: ${#SAMPLES[@]} samples: ${SAMPLES[*]}"
ARGS=()
for S in "${SAMPLES[@]}"; do
    D=data/raw/bmmc_all/$S
    mkdir -p "$D"
    wget -c -q -P "$D" "$S3/$S/atac_fragments.tsv.gz" "$S3/$S/atac_fragments.tsv.gz.tbi"
    gzip -t "$D/atac_fragments.tsv.gz"
    zcat "$D/atac_fragments.tsv.gz" | head -40 | grep -E "^# (id|pipeline_version|reference_fasta_hash|reference_gtf_hash)=" || true
    [ -s "$D/${S}_filtered_feature_bc_matrix.h5" ] || { echo "missing $D/${S}_filtered_feature_bc_matrix.h5 (run 06 first)" >&2; exit 1; }
    ARGS+=(--sample "$S" "$D/${S}_filtered_feature_bc_matrix.h5" "$D/atac_fragments.tsv.gz")
done

mkdir -p "data/raw/$DATASET"
python scripts/aggregate_samples.py --out "data/raw/$DATASET/$PREFIX" "${ARGS[@]}"
ls -lh "data/raw/$DATASET/"
