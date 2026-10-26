#!/bin/bash
#SBATCH -J rebuild_gex
#SBATCH -o logs/%x.%A_%a.o
#SBATCH -e logs/%x.%A_%a.e
#SBATCH -p compute
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH -t 6:0:0

# GEX filtered_feature_bc_matrix.h5 rebuilt from a Cell Ranger ARC BAM (scripts/rebuild_gex_h5.py).
# Reads MODE (probe | genes | rebuild | compare | curated), SAMPLE (pbmc, or a BMMC sample in metadata/bmmc_samples.tsv;
# with --array, task i is the table's row i) and optional REGIONS (contigs joined by +, e.g. chr19+chr21+chr22).
# BMMC samples are written to data/raw/bmmc_all/<sample>/ and checked against their metrics after a rebuild.
set -euo pipefail

PROJECT=$SLURM_SUBMIT_DIR  # submit from the repo root
cd "$PROJECT"
REGIONS=${REGIONS:+${REGIONS//+/,}}  # --export splits on commas, so contigs arrive joined by +
module load miniforge
mamba activate multiome

PBMC=https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k
S3=https://openproblems-bio.s3.amazonaws.com/public/post_competition/multiome
ENA=https://ftp.sra.ebi.ac.uk/vol1/run/SRR176
GENES=data/reference/grch38_2020a_genes.tsv

if [ -z "${SAMPLE:-}" ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    SAMPLE=$(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR == i + 1 {print $1}' metadata/bmmc_samples.tsv)
fi
case ${SAMPLE:-pbmc} in
    pbmc)
        BAM=${PBMC}_gex_possorted_bam.bam
        METRICS_URL=${PBMC}_per_barcode_metrics.csv
        METRICS=data/raw/pbmc/pbmc_unsorted_10k_per_barcode_metrics.csv
        OUT=data/rebuilt/pbmc/pbmc_unsorted_10k_filtered_feature_bc_matrix${REGIONS:+_${REGIONS//,/_}}.h5
        TRUTH=data/raw/pbmc/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5 ;;
    *)
        BAM=$(awk -F'\t' -v s="$SAMPLE" '$1 == s {print $5}' metadata/bmmc_samples.tsv)
        [ -n "$BAM" ] || { echo "sample $SAMPLE is not in metadata/bmmc_samples.tsv" >&2; exit 1; }
        METRICS_URL=$S3/$SAMPLE/per_barcode_metrics.csv
        METRICS=data/raw/bmmc_all/$SAMPLE/per_barcode_metrics.csv
        OUT=data/raw/bmmc_all/$SAMPLE/${SAMPLE}_filtered_feature_bc_matrix.h5 ;;
esac
echo "MODE=$MODE SAMPLE=${SAMPLE:-pbmc} REGIONS=${REGIONS:-all} on $(hostname)"

case $MODE in
    probe)
        for url in "${PBMC}_gex_possorted_bam.bam" "$ENA/SRR17693278/site1_donor01_multiome_gex.possorted_genome_bam.bam.1" \
                   "$ENA/SRR17693276/site1_donor03_multiome_gex.possorted_bam.bam.1"; do
            echo "== $url"; python scripts/rebuild_gex_h5.py --bam "$url" --probe
        done
        for url in "${PBMC}_per_barcode_metrics.csv" "$S3/s1d1/per_barcode_metrics.csv" \
                   https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz; do
            echo "== $url: $(curl -s -I --max-time 30 "$url" | head -1)"
        done ;;
    genes)
        mkdir -p data/reference
        python scripts/reference_genes.py \
            --tarball https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz \
            --fragments data/raw/pbmc/pbmc_unsorted_10k_atac_fragments.tsv.gz --out "$GENES" \
            --check-h5 data/raw/pbmc/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5 ;;
    rebuild)
        mkdir -p "$(dirname "$METRICS")" "$(dirname "$OUT")"
        [ -s "$METRICS" ] || wget -q -O "$METRICS" "$METRICS_URL"
        python scripts/rebuild_gex_h5.py --bam "$BAM" --metrics "$METRICS" --genes "$GENES" --out "$OUT" \
            --library-id "${SAMPLE:-pbmc}" ${REGIONS:+--regions "$REGIONS"}
        if [ "${SAMPLE:-pbmc}" != pbmc ]; then
            python scripts/compare_h5.py --rebuilt "$OUT" --metrics "$METRICS"
        fi ;;
    curated)
        python scripts/compare_curated.py --rebuilt "$OUT" --sample "$SAMPLE" ;;
    compare)
        if [ -n "${TRUTH:-}" ]; then
            python scripts/compare_h5.py --rebuilt "$OUT" --truth "$TRUTH" ${REGIONS:+--regions "$REGIONS"}
        fi
        [ -n "${REGIONS:-}" ] || python scripts/compare_h5.py --rebuilt "$OUT" --metrics "$METRICS" ;;
esac
