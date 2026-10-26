#!/bin/bash
#SBATCH -J create_env_azimuth
#SBATCH -o logs/%x.o%j
#SBATCH -e logs/%x.e%j
#SBATCH -p compute
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH -t 3:0:0

# Builds the "azimuth" env with azimuthpy (pinned GitHub commit), registers its Jupyter kernel and
# downloads the Azimuth human PBMC and human motor cortex references into data/reference/azimuth/.
set -euo pipefail

PROJECT=$SLURM_SUBMIT_DIR  # submit from the repo root
export PYTHONNOUSERSITE=1

module load miniforge

# A failed create leaves a partial prefix that blocks the next attempt.
mamba env remove -n azimuth -y || true
rm -rf "$(conda config --show envs_dirs | awk 'NR == 2 {print $2}')/azimuth"

mamba env create -f "$PROJECT/envs/azimuth.yml"
mamba activate azimuth
pip install --no-deps "git+https://github.com/KalinNonchev/azimuthpy.git@2465452"  # 0.1.1

REF=$PROJECT/data/reference/azimuth
for ENTRY in human_pbmc:4546839 human_motor_cortex:4546932; do
    NAME=${ENTRY%%:*}
    RECORD=${ENTRY##*:}
    mkdir -p "$REF/$NAME"
    for F in ref.Rds idx.annoy; do
        wget -q -O "$REF/$NAME/$F" "https://zenodo.org/api/records/$RECORD/files/$F/content"
    done
done
ls -lh "$REF"/*

python - <<'PY'
import numpy, pandas, scanpy, anndata, rpy2
import azimuthpy  # loads Seurat and Azimuth through rpy2
import rpy2.robjects as ro
for m in (numpy, pandas, scanpy, anndata, rpy2):
    print(m.__name__, m.__version__)
print(ro.r('paste("R", getRversion(), "| Seurat", packageVersion("Seurat"), "| Azimuth", packageVersion("Azimuth"))')[0])
PY

# Annotation levels in each reference; the notebooks use celltype.l1 (PBMC) and subclass (motor cortex).
Rscript -e 'suppressMessages(library(Seurat)); for (d in commandArgs(TRUE)) { ref <- Azimuth::LoadReference(d)$map; cat(basename(d), ":", colnames(ref[[]]), "\n") }' \
    "$REF/human_pbmc" "$REF/human_motor_cortex"

python -m ipykernel install --user --name azimuth --display-name "Python (azimuth)"
mamba env export -n azimuth | grep -v '^prefix:' > "$PROJECT/envs/azimuth.lock.yml"
echo "azimuth env created."
