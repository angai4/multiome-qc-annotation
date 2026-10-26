#!/bin/bash
#SBATCH -J create_env_multiome
#SBATCH -o logs/%x.o%j
#SBATCH -e logs/%x.e%j
#SBATCH -p compute
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH -t 2:0:0

# Builds the "multiome" env and installs AMULET v1.1 into tools/AMULET.
set -euo pipefail

PROJECT=$SLURM_SUBMIT_DIR  # submit from the repo root
export PYTHONNOUSERSITE=1

module load miniforge

# A failed create leaves a partial prefix that blocks the next attempt.
mamba env remove -n multiome -y || true
rm -rf "$(conda config --show envs_dirs | awk 'NR == 2 {print $2}')/multiome"

mamba env create -f "$PROJECT/envs/multiome.yml"
mamba activate multiome

AMULET_DIR=$PROJECT/tools/AMULET
rm -rf "$AMULET_DIR"
mkdir -p "$AMULET_DIR"
cd "$AMULET_DIR"
wget -q https://github.com/UcarLab/AMULET/releases/download/v1.1/AMULET-v1.1.zip
wget -q https://github.com/UcarLab/AMULET/releases/download/v1.0/RestrictionRepeatLists.zip
unzip -oq AMULET-v1.1.zip
unzip -oq RestrictionRepeatLists.zip
rm -rf __MACOSX ./*.zip
# np.object was removed in numpy 1.24; it was an alias of the builtin object.
sed -i 's/np\.object\b/object/g' AMULET.py peakoverlap.py FragmentFileOverlapCounter.py

python - <<'PY'
import scanpy, anndata, snapatac2, numpy, pandas, scipy, statsmodels, igraph, leidenalg
for m in (scanpy, anndata, snapatac2, numpy, pandas, scipy, statsmodels, igraph):
    print(m.__name__, m.__version__)
print("leidenalg", leidenalg.version)
print("hg38 annotation:", snapatac2.genome.hg38.annotation)
PY

python -m ipykernel install --user --name multiome --display-name "Python (multiome)"
mamba env export -n multiome | grep -v '^prefix:' > "$PROJECT/envs/multiome.lock.yml"
echo "multiome env created."
