"""Per-cell ATAC QC metrics for the Cell Ranger-called cells of a multiome sample.

Unique fragments and TSS enrichment (SnapATAC2), and nucleosome signal (Signac's definition: mononucleosomal
fragments, 147-294 bp, over nucleosome-free ones, under 147 bp). Writes a TSV with one row per called cell, in the
h5's barcode order.

  python scripts/atac_metrics.py --h5 x_filtered_feature_bc_matrix.h5 --fragments x_atac_fragments.tsv.gz \
      --out atac_metrics.tsv [--n-jobs 8]
"""

import argparse

import h5py
import numpy as np
import snapatac2 as snap


def count_fragments(lengths, low, high):
    """Fragments per cell with low <= length < high, counted without copying the matrix."""
    in_range = (lengths.data >= low) & (lengths.data < high)
    counts = np.zeros(lengths.shape[0], dtype=np.int64)
    has_fragments = np.diff(lengths.indptr) > 0
    counts[has_fragments] = np.add.reduceat(in_range, lengths.indptr[:-1][has_fragments], dtype=np.int64)
    return counts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5", required=True)
    p.add_argument("--fragments", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--n-jobs", type=int, default=8)
    args = p.parse_args()

    with h5py.File(args.h5) as f:
        barcodes = [b.decode() for b in f["matrix/barcodes"][:]]

    frags = snap.pp.import_fragments(
        args.fragments,
        chrom_sizes=snap.genome.hg38,
        whitelist=barcodes,
        sorted_by_barcode=False,
        min_num_fragments=0,
        n_jobs=args.n_jobs,
    )
    snap.metrics.tsse(frags, snap.genome.hg38, n_jobs=args.n_jobs)

    # one row per cell, values are fragment lengths
    lengths = frags.obsm["fragment_paired"]
    mononucleosome = count_fragments(lengths, 147, 294)
    nucleosome_free = count_fragments(lengths, 0, 147)
    frags.obs["nucleosome_signal"] = mononucleosome / np.maximum(nucleosome_free, 1)

    metrics = frags.obs[["n_fragment", "tsse", "nucleosome_signal"]].reindex(barcodes)
    metrics.index.name = "barcode"
    metrics.to_csv(args.out, sep="\t")
    print(f"{len(barcodes)} cells, {metrics['n_fragment'].isna().sum()} without fragments -> {args.out}")
    print(metrics.describe().T.to_string())


if __name__ == "__main__":
    main()
