# multiome-qc-annotation

QC and cell-type annotation of 10x multiome (RNA + ATAC) data. Each dataset is QC'd on its own (RNA and ATAC
metrics, Scrublet and AMULET doublets), clustered with scanpy and annotated with Azimuth.

| Dataset | Data | Notebook |
| --- | --- | --- |
| pbmc | 10x PBMC 10k (no sorting) | `notebooks/pbmc_v3.ipynb` |
| brain | 10x human brain 3k | `notebooks/brain_v3.ipynb` |
| bmmc | NeurIPS 2021 bone marrow, samples s1d1 + s1d3 | `notebooks/bmmc_s1d1_s1d3_v3.ipynb` |
| bmmc_all | NeurIPS 2021 bone marrow, all 13 samples | `notebooks/bmmc_all_v3.ipynb` |

## Layout

| Path | Contents |
| --- | --- |
| `jobs/` | Slurm jobs, numbered in run order |
| `scripts/` | Python run by the jobs |
| `notebooks/` | One notebook per dataset |
| `metadata/` | `bmmc_samples.tsv`: the 13 BMMC samples, in the order they are combined (don't reorder) |
| `envs/` | Conda envs: `<name>.yml` (packages) and `<name>.lock.yml` (exact versions used) |

`data/`, `figures/`, `tools/` and `logs/` are created by the jobs and notebooks and are not in git.

## Running it

Built for QMUL's Apocrita (Slurm, conda/mamba); on another cluster, change the `#SBATCH -p` partitions. Submit jobs
from the repo root. Large files go to scratch through two links:

```bash
mkdir -p data/raw /gpfs/scratch/$USER/multiome_data/raw/bmmc_all
ln -s /gpfs/scratch/$USER/multiome_data/raw/bmmc_all data/raw/bmmc_all
ln -s /gpfs/scratch/$USER/multiome_data data/scratch
```

| Job | Does | Datasets |
| --- | --- | --- |
| `01_create_env` | `multiome` env (for the jobs) and AMULET | all |
| `02_download` | 10x matrices and ATAC fragments | pbmc, brain |
| `03_amulet` | ATAC doublets (AMULET) | all |
| `04_atac_metrics` | Fragments, TSS enrichment, nucleosome signal per cell | all |
| `05_create_env_azimuth` | `azimuth` env (for the notebooks and job 08) and Azimuth references | all |
| `06_rebuild_gex` | GEX matrix rebuilt from each sample's BAM (none is public) | bmmc |
| `07_bmmc_aggr` | Downloads fragments, combines the samples | bmmc |
| `08_azimuth` | Azimuth as a batch job, too large for the notebook | bmmc_all |
| `09_run_notebook` | Runs a notebook as a batch job | all |

Jobs that take an array index use 0 = pbmc, 1 = brain, 2 = bmmc, 3 = bmmc_all. `jobs/run_bmmc_all.sh` queues
06 → 07 → 03 + 04 for bmmc_all, and `jobs/run_notebooks.sh` re-runs all four notebooks (both are dry runs unless
given `--submit`).

Run the notebooks from `notebooks/` with the "Python (azimuth)" kernel. `bmmc_all_v3` stops before annotation for
job 08, then continues from its output. Each notebook also saves every figure as an SVG in `figures/<notebook>/`.

To reproduce the exact package versions, build the envs from `envs/*.lock.yml` instead of `envs/*.yml`.

## Outputs

In `data/processed/<dataset>/`:

- `rna_annotated.h5ad`: final cells × genes with QC metrics, doublet calls, clusters, Azimuth labels, HVG flags
  and UMAP.
- `atac_annotated.h5ad`: ATAC peak counts for the same cells (pbmc and brain; BMMC has no public peak matrix).
- `atac_metrics.tsv` and `amulet/`: the per-cell ATAC QC and doublet calls from jobs 03 and 04.

## Data sources

- **pbmc, brain:** 10x Genomics datasets *PBMC from a Healthy Donor - No Cell Sorting (10k)* and *Frozen Human
  Healthy Brain Tissue (3k)*, Cell Ranger ARC 2.0.0.
- **bmmc:** NeurIPS 2021 multimodal single-cell competition (Luecken et al. 2021; GEO GSE194122). GEX BAMs from ENA
  study PRJNA799242 (the 13 "Multiome (GEX)" runs, listed with URLs and md5s in `metadata/bmmc_samples.tsv`); ATAC
  fragments from the Open Problems S3 bucket.
- **References:** 10x GRCh38 2020-A genes; Azimuth human PBMC (Zenodo 4546839; Hao et al. 2021) and human motor
  cortex (Zenodo 4546932; Bakken et al. 2021).
- **Methods:** [scanpy](https://scanpy.readthedocs.io) (Wolf et al. 2018), Scrublet (Wolock et al. 2019),
  [SnapATAC2](https://github.com/kaizhang/SnapATAC2) (Zhang et al. 2024), [AMULET](https://github.com/UcarLab/AMULET)
  (Thibodeau et al. 2021), [Azimuth](https://azimuth.hubmapconsortium.org/) / Seurat (Hao et al. 2021).

Exact download URLs are in jobs 02, 05, 06 and 07.

## Notes

- The `azimuth` env pins R 4.3 and rpy2 < 3.6, because azimuthpy 0.1.1 calls `pandas2ri.activate()`, which rpy2 3.6
  removed.
- AMULET is v1.1 with `np.object` replaced by `object` (removed in numpy 1.24).

## License

[MIT](LICENSE)
