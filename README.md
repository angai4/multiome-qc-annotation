# multiome-qc-annotation

10x multiome (RNA + ATAC) datasets, each QC'd on its own and then annotated with Azimuth: **pbmc**, **brain**
(10x public datasets) and **bmmc** (NeurIPS 2021 BMMC multiome; two samples, then all 13).

Built for QMUL's Apocrita (Slurm, `module load miniforge`); the jobs should run on any Slurm cluster with
conda/mamba after editing the `#SBATCH -p` partition and pointing `data/scratch` at a large disk (see Setup).

## Layout

| Path | Contents |
| --- | --- |
| `jobs/` | Slurm jobs, numbered in the order they run; `run_*.sh` queue several as a chain |
| `scripts/` | Python run by the jobs |
| `metadata/` | `bmmc_samples.tsv`: the 13 BMMC samples (site, donor, ENA run, BAM URL and md5), in aggregation order |
| `notebooks/` | One notebook per dataset: QC, clustering and Azimuth annotation |
| `envs/` | Conda envs: `<name>.yml` (the packages asked for) and `<name>.lock.yml` (exact versions installed) |
| `data/`, `figures/`, `tools/`, `logs/` | Not in git: downloads, references, outputs, figure SVGs, AMULET and job logs, all made by the jobs and notebooks |

## Setup

```bash
git clone https://github.com/<you>/multiome-qc-annotation.git && cd multiome-qc-annotation

# data/ is not in git. Two links put the large files on scratch, outside the 100 GB home quota:
# data/raw/bmmc_all (the 13 BMMC samples and their aggregate) and data/scratch (bmmc_all's clustered.h5ad and
# azimuth.h5ad). Everything else under data/ is created by the jobs.
mkdir -p data/raw /gpfs/scratch/$USER/multiome_data/raw/bmmc_all
ln -s /gpfs/scratch/$USER/multiome_data/raw/bmmc_all data/raw/bmmc_all
ln -s /gpfs/scratch/$USER/multiome_data data/scratch
```

Submit every job **from the repo root** (`sbatch jobs/...`): jobs find the project through `$SLURM_SUBMIT_DIR`
and write their logs to `logs/`.

**Envs.** `envs/<name>.yml` lists the packages, mostly unpinned; `envs/<name>.lock.yml` records the exact version
and build of every package that was installed. Jobs 01 and 05 build from the `.yml` (so a rebuild gets current
versions) and then overwrite the lock file. To reproduce the versions used here, create the env from the lock file
instead (`mamba env create -n multiome -f envs/multiome.lock.yml`) and run the rest of the job by hand.

## Steps

Array index, wherever a job takes one: 0 = pbmc, 1 = brain, 2 = bmmc, 3 = bmmc_all.

### pbmc and brain

- `sbatch jobs/01_create_env.sbatch`: builds the `multiome` env, registers the "Python (multiome)" kernel and
  installs AMULET into `tools/`.
- `sbatch jobs/02_download.sbatch`: 10x h5 and ATAC fragments into `data/raw/<dataset>/`.
- `sbatch jobs/03_amulet.sbatch`: AMULET multiplets into `data/processed/<dataset>/amulet/`, called per GEM well.
- `sbatch jobs/04_atac_metrics.sbatch`: per-cell ATAC QC (unique fragments, TSS enrichment, nucleosome signal) into
  `data/processed/<dataset>/atac_metrics.tsv`.
- `sbatch jobs/05_create_env_azimuth.sbatch`: builds the `azimuth` env (R + Seurat + Azimuth +
  [azimuthpy](https://github.com/KalinNonchev/azimuthpy) at a pinned commit), registers the "Python (azimuth)" kernel
  and downloads the Azimuth human PBMC and human motor cortex references into `data/reference/azimuth/`.
- Notebooks `notebooks/{pbmc,brain}_v3.ipynb` (kernel "Python (azimuth)", 16 GB RAM; run from `notebooks/`): follow
  the scanpy [clustering tutorial](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html) with
  ATAC QC added. RNA QC (genes, counts, % mito, min cells) and ATAC QC (fragments, TSSE, nucleosome signal); keep
  cells passing both; Scrublet and AMULET doublets; normalise, log1p, HVGs, PCA, UMAP, Leiden; re-assess QC; Azimuth
  annotation (PBMC `celltype.l1`, motor cortex `subclass`).

### BMMC (NeurIPS 2021)

No multiome `filtered_feature_bc_matrix.h5` is public for these samples, so the GEX matrix is rebuilt from each
sample's Cell Ranger ARC BAM. `jobs/run_bmmc_all.sh [--submit]` queues 06 → 07 → 03 + 04 for `bmmc_all` as a
dependency chain (dry run unless `--submit`).

- `jobs/06_rebuild_gex.sbatch`: rebuilds a sample's GEX matrix from its BAM (`scripts/rebuild_gex_h5.py`). It sums
  the reads Cell Ranger marked as counted UMIs (`xf` bit 8), per `CB` and `GX`, for its called cells (`is_cell`).
  This is not a Cell Ranger rerun; it reproduces our PBMC `.h5` exactly. Run with
  `--export=ALL,MODE=...,SAMPLE=...[,REGIONS=chr19+chr21+chr22]`, or as an array over `metadata/bmmc_samples.tsv`.
  - `MODE=genes`: the gene table from the 10x reference tarball (once, first).
  - `MODE=rebuild`: rebuild, then check per-cell totals against Cell Ranger's metrics.
  - `MODE=compare` / `MODE=curated`: compare a rebuild with a real `.h5` (PBMC) or the curated NeurIPS object.
- `jobs/07_bmmc_aggr.sbatch`: downloads each sample's ATAC fragments next to its rebuilt matrix
  (`data/raw/bmmc_all/<sample>/`) and combines the samples `cellranger-arc aggr` style (barcode suffix `-1`, `-2`, …
  in table order; no depth normalisation) into `data/raw/<dataset>/`. `--export=ALL,DATASET=bmmc` for the
  two-sample test (s1d1 + s1d3), default `bmmc_all`.
- Then `03_amulet` and `04_atac_metrics` with array index 2 or 3.
- Notebooks `notebooks/bmmc_s1d1_s1d3_v3.ipynb` (annotates in the notebook) and `notebooks/bmmc_all_v3.ipynb`
  (peaks at about 15 GB RAM). `bmmc_all_v3` saves its final cells to `data/scratch/processed/bmmc_all/clustered.h5ad` and
  stops for the next job; once that has run, restart from the first cell's imports and continue.
- `jobs/08_azimuth.sbatch`: Azimuth annotation as a batch job (`scripts/run_azimuth.py`), for datasets too large to
  annotate in the notebook. Writes `azimuth.h5ad` (Azimuth's columns and `ref.umap`) next to `clustered.h5ad`, and
  the notebook joins it back on. 1 core, 160 GB.

### Re-running the notebooks

`jobs/run_notebooks.sh [--submit]` re-runs all four notebooks as batch jobs (`jobs/09_run_notebook.sbatch`, in place;
`bmmc_all` separately with a longer walltime, using the existing `azimuth.h5ad`). Every figure a notebook shows is also saved as an SVG to
`figures/<notebook>/` (not in git), named after its `In[]` number and title, with point clouds rasterised at 300 dpi.

## Outputs

Per dataset, in `data/processed/<dataset>/`:

| File | Written by | Contents |
| --- | --- | --- |
| `amulet/batch_<well>/MultipletProbabilities.txt` | `03_amulet` | AMULET p- and q-values per cell, per GEM well |
| `atac_metrics.tsv` | `04_atac_metrics` | Unique fragments, TSS enrichment and nucleosome signal per called cell |
| `rna_annotated.h5ad` | notebook | Final cells × genes: QC metrics, doublet calls, clusters, Azimuth predictions, HVG flags, UMAP |
| `atac_annotated.h5ad` | notebook | ATAC peak counts for the same cells, with the same `obs` |

`bmmc_all` also keeps `clustered.h5ad` and `azimuth.h5ad` in `data/scratch/processed/bmmc_all/`. Figures are in the
notebooks' outputs, and as SVGs in `figures/<notebook>/` (local only).

## Data sources

| Dataset | Source |
| --- | --- |
| pbmc | 10x Genomics, *PBMC from a Healthy Donor - No Cell Sorting (10k)*, Cell Ranger ARC 2.0.0 (`pbmc_unsorted_10k`) |
| brain | 10x Genomics, *Frozen Human Healthy Brain Tissue (3k)*, Cell Ranger ARC 2.0.0 (`human_brain_3k`) |
| bmmc | NeurIPS 2021 multimodal single-cell competition, bone marrow multiome (Luecken et al. 2021; GEO GSE194122). GEX BAMs from ENA (runs in `metadata/bmmc_samples.tsv`); ATAC fragments and per-barcode metrics from the Open Problems S3 bucket (`openproblems-bio/public/post_competition/multiome`) |

| Reference | Source |
| --- | --- |
| Genes (GRCh38) | 10x `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` |
| Azimuth human PBMC | Zenodo record 4546839 (Hao et al. 2021) |
| Azimuth human motor cortex | Zenodo record 4546932 (Bakken et al. 2021) |
| AMULET restriction list | AMULET v1.0 release, `RestrictionRepeatLists` (hg38) |

Exact URLs are in `jobs/02_download.sbatch`, `jobs/05_create_env_azimuth.sbatch`, `jobs/06_rebuild_gex.sbatch` and
`jobs/07_bmmc_aggr.sbatch`. `metadata/bmmc_samples.tsv` lists the 13 "Site *N* Donor *M* Multiome (GEX)" runs of ENA
study PRJNA799242 (SRA SRP356158), with ENA's submitted-file URL and md5 for each BAM; its row order sets the
barcode suffixes (`-1` … `-13`), so don't reorder it.

Methods: [scanpy](https://scanpy.readthedocs.io) (Wolf et al. 2018), Scrublet (Wolock et al. 2019),
[SnapATAC2](https://github.com/kaizhang/SnapATAC2) (Zhang et al. 2024), [AMULET](https://github.com/UcarLab/AMULET)
(Thibodeau et al. 2021), [Azimuth](https://azimuth.hubmapconsortium.org/) / Seurat (Hao et al. 2021).

## Notes

The `azimuth` env pins R 4.3 and rpy2 < 3.6: rpy2 3.6 removed `pandas2ri.activate()`, which azimuthpy
0.1.1 calls on import.

AMULET is the v1.1 release with `np.object` replaced by `object` (removed in numpy 1.24).
