#!/bin/bash
# Re-runs the v3 notebooks as batch jobs (jobs/09_run_notebook.sbatch): pbmc, brain and bmmc as one array, bmmc_all
# on its own with more memory. Dry run by default (prints the sbatch lines); --submit queues them.
set -euo pipefail
cd "$(dirname "$0")/.."

SUBMIT=0
[[ "${1:-}" == "--submit" ]] && SUBMIT=1
# bmmc_all joins the annotation from an earlier jobs/08_azimuth.sbatch run
AZIMUTH=data/scratch/processed/bmmc_all/azimuth.h5ad
[ -s "$AZIMUTH" ] || { echo "missing $AZIMUTH (run jobs/08_azimuth.sbatch first)" >&2; exit 1; }
mkdir -p logs

run() {  # prints the sbatch line; with --submit also queues it
    echo "+ sbatch $*"
    if (( SUBMIT )); then sbatch "$@"; fi
}
run --array=0-2 jobs/09_run_notebook.sbatch
run --array=3 --mem=32G -t 12:0:0 jobs/09_run_notebook.sbatch
(( SUBMIT )) || echo -e "\nDRY RUN. Re-run with --submit to queue."
