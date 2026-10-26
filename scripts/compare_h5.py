"""Compare a rebuilt GEX matrix with a Cell Ranger filtered_feature_bc_matrix.h5 (its Gene Expression block).

Barcodes must be the same set and genes the same, in the same order. Every count must be identical. With
--regions, only genes whose interval lies on those contigs are compared (for a rebuild that read only them).

--metrics (no --truth) instead checks every cell's total UMIs against Cell Ranger's per_barcode_metrics.csv:
gex_exonic_umis + gex_intronic_umis, which equals the real matrix's per-cell total for all PBMC cells
(gex_umis_count is 0 in Cell Ranger ARC 2.0.0 output).

  python scripts/compare_h5.py --rebuilt x.h5 --truth y.h5 [--regions chr19,chr21,chr22]
  python scripts/compare_h5.py --rebuilt x.h5 --metrics per_barcode_metrics.csv
"""

import argparse

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_gex(path):
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        x = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:]))
        ft = g["features"]
        gex = ft["feature_type"][:] == b"Gene Expression"
        return (x[np.flatnonzero(gex)].tocsc(), g["barcodes"][:].astype(str),
                ft["id"][:][gex].astype(str), ft["interval"][:][gex].astype(str))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--rebuilt", required=True)
    p.add_argument("--truth")
    p.add_argument("--metrics")
    p.add_argument("--regions")
    args = p.parse_args()

    x, bc_x, ids_x, _ = read_gex(args.rebuilt)
    if args.metrics:
        m = pd.read_csv(args.metrics).set_index("barcode")
        assert set(m.index[m["is_cell"] == 1]) == set(bc_x), "rebuilt barcodes are not the called cells"
        expected = (m["gex_exonic_umis"] + m["gex_intronic_umis"]).loc[bc_x].to_numpy()
        total = np.asarray(x.sum(axis=0)).ravel()
        bad = np.flatnonzero(total != expected)
        print(f"{len(bc_x)} cells; UMIs rebuilt {total.sum():,}, metrics {expected.sum():,}; cells differing: {len(bad)}")
        for i in bad[:20]:
            print(f"  {bc_x[i]}: rebuilt {total[i]} metrics {expected[i]}")
        print("PER-CELL TOTALS EXACT" if len(bad) == 0 else "MISMATCH")
        raise SystemExit(0 if len(bad) == 0 else 1)

    y, bc_y, ids_y, interval = read_gex(args.truth)
    print(f"rebuilt {x.shape}, truth {y.shape}")
    assert set(bc_x) == set(bc_y), f"barcodes differ: {len(set(bc_x) - set(bc_y))} only rebuilt, {len(set(bc_y) - set(bc_x))} only truth"
    assert (ids_x == ids_y).all(), "gene ids or order differ"
    y = y[:, pd.Index(bc_y).get_indexer(bc_x)]  # truth columns in the rebuilt barcode order
    rows = np.arange(len(ids_x))
    if args.regions:
        rows = np.flatnonzero(np.isin([iv.split(":")[0] for iv in interval], args.regions.split(",")))
        print(f"comparing {len(rows)} genes on {args.regions}")
    xs, ys = x[rows].tocsr(), y[rows].tocsr()
    diff = (xs - ys).tocoo()
    diff.eliminate_zeros()
    print(f"UMIs: rebuilt {int(xs.sum()):,}, truth {int(ys.sum()):,}; differing entries: {diff.nnz:,}")
    for r, c, v in list(zip(diff.row, diff.col, diff.data))[:20]:
        print(f"  gene {ids_x[rows[r]]} cell {bc_x[c]}: rebuilt {xs[r, c]} truth {ys[r, c]} (diff {v})")
    print("EXACT MATCH" if diff.nnz == 0 else "MISMATCH")
    raise SystemExit(0 if diff.nnz == 0 else 1)


if __name__ == "__main__":
    main()
