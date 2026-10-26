#!/bin/bash
# Builds the bmmc_all dataset (all 13 NeurIPS 2021 multiome samples): GEX rebuilds (06, array over
# metadata/bmmc_samples.tsv, 3 at a time) -> fragments + aggregation (07) -> AMULET (03) and ATAC metrics (04) with
# more memory than their headers ask for. Dry run by default (prints the sbatch lines); --submit queues them.
set -euo pipefail
cd "$(dirname "$0")/.."

SUBMIT=0
[[ "${1:-}" == "--submit" ]] && SUBMIT=1
N=$(($(wc -l < metadata/bmmc_samples.tsv) - 1))

for f in metadata/bmmc_samples.tsv data/reference/grch38_2020a_genes.tsv; do
    [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
done
[ -d data/raw/bmmc_all/ ] || { echo "data/raw/bmmc_all (scratch symlink) is missing" >&2; exit 1; }
echo "samples ($N): $(awk -F'\t' 'NR > 1 {printf "%s ", $1}' metadata/bmmc_samples.tsv)"
echo "raw -> $(readlink -f data/raw/bmmc_all); processed -> data/processed/bmmc_all"
echo

run() {  # prints the sbatch line; with --submit also queues it and returns the job id
    echo "+ sbatch $*" >&2
    if (( SUBMIT )); then sbatch --parsable "$@"; else echo JOB; fi
}
R=$(run --array=1-"$N"%3 --export=ALL,MODE=rebuild jobs/06_rebuild_gex.sbatch)
A=$(run --dependency=afterok:"$R" --export=ALL,DATASET=bmmc_all jobs/07_bmmc_aggr.sbatch)
run --dependency=afterok:"$A" --array=3 --mem-per-cpu=96G -t 24:0:0 jobs/03_amulet.sbatch >/dev/null
run --dependency=afterok:"$A" --array=3 --mem-per-cpu=25G -t 6:0:0 jobs/04_atac_metrics.sbatch >/dev/null
(( SUBMIT )) || echo -e "\nDRY RUN. Re-run with --submit to queue."
