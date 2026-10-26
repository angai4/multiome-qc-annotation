"""Compare a rebuilt GEX matrix with the NeurIPS 2021 curated MultiVI h5ad, read remotely (range requests).

The curated file keeps the authors' QC'd cells and 25,629 genes, with raw `counts`. For one sample, every
(cell, gene) count it holds must equal the rebuilt matrix's. Its cells must all be Cell Ranger called cells.
Only the sample's block of rows is downloaded.

  python scripts/compare_curated.py --rebuilt s1d1_filtered_feature_bc_matrix.h5 --sample s1d1 [--curated URL]
"""

import argparse
import io
import urllib.request

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

CURATED = ("https://openproblems-bio.s3.amazonaws.com/public/post_competition/multiome/"
           "bmmc_multiome_multivi_neurips21_curated.h5ad")


class RangeFile(io.RawIOBase):
    """A read-only, seekable view of a remote file, fetched with HTTP range requests."""

    def __init__(self, url):
        self.url, self.pos = url, 0
        head = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"))
        self.size = int(head.headers["Content-Length"])

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def readinto(self, buffer):
        if self.pos >= self.size:
            return 0
        end = min(self.pos + len(buffer), self.size) - 1
        request = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        data = urllib.request.urlopen(request).read()
        buffer[:len(data)] = data
        self.pos += len(data)
        return len(data)


def categorical(group):
    return np.asarray(group["categories"][:]).astype(str)[group["codes"][:]]


def read_rebuilt(path):
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        x = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:]))
        return x.T.tocsr(), g["barcodes"][:].astype(str), g["features"]["id"][:].astype(str)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--rebuilt", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--curated", default=CURATED)
    args = p.parse_args()

    rebuilt, barcodes, gene_ids = read_rebuilt(args.rebuilt)
    f = h5py.File(io.BufferedReader(RangeFile(args.curated), buffer_size=1 << 24), "r")
    rows = np.flatnonzero(categorical(f["obs"]["batch"]) == args.sample)
    assert len(rows) and rows[-1] - rows[0] + 1 == len(rows), "sample rows are not one contiguous block"
    names = np.asarray(f["obs"][f["obs"].attrs["_index"]][rows[0]:rows[-1] + 1]).astype(str)
    cells = np.array([n.split(":", 1)[1] + "-1" for n in names])  # s1d1_multiome:ACGT... -> ACGT...-1
    features = categorical(f["var"]["feature_types"])
    ids = f["var"]["gene_ids"]
    ids = categorical(ids) if isinstance(ids, h5py.Group) else np.asarray(ids[:]).astype(str)
    gex = np.flatnonzero(features == "GEX")
    print(f"curated {args.sample}: {len(rows)} cells, {len(gex)} genes; rebuilt: {len(barcodes)} cells, {len(gene_ids)} genes")

    counts = f["layers"]["counts"]
    indptr = counts["indptr"][rows[0]:rows[-1] + 2]
    lo, hi = int(indptr[0]), int(indptr[-1])
    block = sp.csr_matrix((counts["data"][lo:hi], counts["indices"][lo:hi], indptr - lo),
                          shape=(len(rows), counts.attrs["shape"][1]))
    curated = block[:, gex]
    assert np.array_equal(curated.data, np.round(curated.data)), "curated counts are not integers"

    cell_pos = pd.Index(barcodes).get_indexer(cells)
    gene_pos = pd.Index(gene_ids).get_indexer(ids[gex])
    print(f"curated cells among Cell Ranger called cells: {(cell_pos >= 0).sum()} / {len(cells)}; "
          f"curated genes in the gene table: {(gene_pos >= 0).sum()} / {len(gex)}")
    assert (cell_pos >= 0).all() and (gene_pos >= 0).all()
    ours = rebuilt[cell_pos][:, gene_pos]
    diff = (ours - curated.astype(ours.dtype)).tocoo()
    diff.eliminate_zeros()
    print(f"UMIs on the curated cells x genes: rebuilt {int(ours.sum()):,}, curated {int(curated.sum()):,}; "
          f"differing entries: {diff.nnz:,}")
    for r, c, v in list(zip(diff.row, diff.col, diff.data))[:20]:
        print(f"  cell {cells[r]} gene {ids[gex][c]}: rebuilt {ours[r, c]} curated {curated[r, c]}")
    print("IDENTICAL TO CURATED" if diff.nnz == 0 else "MISMATCH")
    raise SystemExit(0 if diff.nnz == 0 else 1)


if __name__ == "__main__":
    main()
