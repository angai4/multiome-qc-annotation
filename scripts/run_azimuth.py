"""Azimuth cell-type annotation of a notebook's final cells, outside the notebook.

Reads layers["counts"] of --input (clustered.h5ad, written by the notebook before annotation) and writes --out
(azimuth.h5ad): Azimuth's obs columns and obsm["ref.umap"], same cells in the same order. Run in the azimuth env with
R_HOME set (jobs/08_azimuth.sbatch).

  python scripts/run_azimuth.py --input clustered.h5ad --reference data/reference/azimuth/human_pbmc --out azimuth.h5ad
"""

import argparse
import time

import anndata as ad
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from azimuthpy.python_wrapper import adata_to_seurat, extract_meta_data


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    t0 = time.time()
    full = ad.read_h5ad(args.input)
    # only the counts go to Azimuth; azimuthpy reads layers["raw_counts"] whatever use_layer says
    query = ad.AnnData(obs=pd.DataFrame(index=full.obs_names), var=pd.DataFrame(index=full.var_names),
                       layers={"raw_counts": full.layers["counts"]})
    del full
    print(f"{query.n_obs:,} cells x {query.n_vars:,} genes from {args.input}; reference {args.reference}", flush=True)

    # azimuthpy's annotate_adata_azimuth calls RunAzimuth from Python, which puts the whole query object inside the
    # top-level R call; Seurat's LogSeuratCommand then deparses it, which never finishes at this size. Calling
    # RunAzimuth from R on the seurat_query variable (created by adata_to_seurat) keeps that call small.
    adata_to_seurat(query)
    ro.r(f'mapped <- Azimuth::RunAzimuth(query = seurat_query, reference = "{args.reference}")')
    meta, umap = extract_meta_data(ro.r("mapped"), "ref.umap")
    umap_cells = np.asarray(ro.r('rownames(Seurat::Embeddings(mapped[["ref.umap"]]))'))
    assert (umap_cells == query.obs_names).all()
    query.obs = query.obs.join(meta)

    out = ad.AnnData(obs=query.obs, obsm={"ref.umap": np.asarray(umap)})
    assert out.obs_names.equals(query.obs_names) and out.obsm["ref.umap"].shape == (out.n_obs, 2)
    assert out.obs.filter(like="predicted.").notna().all().all() and out.obs["mapping.score"].notna().all()
    out.write_h5ad(args.out)
    print(f"wrote {args.out}: columns {list(out.obs.columns)}; {(time.time() - t0) / 60:.0f} min")


if __name__ == "__main__":
    main()
