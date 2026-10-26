"""Rebuild the GEX part of a Cell Ranger ARC filtered_feature_bc_matrix.h5 from its GEX BAM.

This is not a Cell Ranger rerun. The BAM already records, for every read, Cell Ranger's corrected barcode (CB),
its gene assignment (GX) and whether the read is the one representative read counted as a UMI (xf bit 8).
Summing those reads per (cell, gene) for Cell Ranger's own called cells (is_cell in per_barcode_metrics.csv)
replays only its last step, which turns molecules into the matrix.

Cells: is_cell == 1 barcodes from --metrics, sorted. Genes: --genes (scripts/reference_genes.py), in reference
order. The BAM can be a local path or an http(s)/ftp URL. It is read one contig at a time through its index
(retrying a contig whose stream breaks); --regions reads only those contigs.

  python scripts/rebuild_gex_h5.py --bam URL --metrics per_barcode_metrics.csv --genes genes.tsv --out x.h5 \
      [--regions chr19,chr21,chr22] [--library-id name]
  python scripts/rebuild_gex_h5.py --bam URL --probe      (header and the tags of the first 1,000 reads)
"""

import argparse
import collections
import time
from array import array

import h5py
import numpy as np
import pandas as pd
import pysam
import scipy.sparse as sp

COUNTED = 8  # xf bit 8: "representative for a molecule and is treated as a UMI count"


def probe(bam_url: str):
    with pysam.AlignmentFile(bam_url, "rb") as bam:
        header = bam.header.to_dict()
        for pg in header.get("PG", []):
            print("@PG", {k: pg[k] for k in ("ID", "PN", "VN") if k in pg})
        for co in header.get("CO", [])[:10]:
            print("@CO", co[:150])
        print("contigs:", len(header.get("SQ", [])), [sq["SN"] for sq in header["SQ"][:3]])
        tags, xf = collections.Counter(), collections.Counter()
        for i, read in enumerate(bam.fetch(until_eof=True)):
            if i == 1000:
                break
            tags.update(t for t, _ in read.get_tags())
            if read.has_tag("xf"):
                xf[read.get_tag("xf")] += 1
        print("tags on the first 1,000 reads:", dict(tags))
        print("xf values:", dict(xf))


MAX_TRIES = 5


def count_contig(bam_url, contig, cell_index, gene_index):
    """Counted reads of one contig, as (cell indices, gene indices, tally of skipped reads)."""
    cells, genes, tally = array("i"), array("i"), collections.Counter()
    with pysam.AlignmentFile(bam_url, "rb") as bam:
        for read in bam.fetch(contig):
            tally["reads seen"] += 1
            if not read.has_tag("xf"):
                tally["no xf"] += 1
                continue
            if not read.get_tag("xf") & COUNTED:
                tally["not counted (xf bit 8 unset)"] += 1
                continue
            if not read.has_tag("CB"):
                tally["counted, no CB"] += 1
                continue
            cell = cell_index.get(read.get_tag("CB"))
            if cell is None:
                tally["counted, barcode not a called cell"] += 1
                continue
            gx = read.get_tag("GX") if read.has_tag("GX") else ""
            if ";" in gx or gx == "":
                tally["counted, GX empty or multi-gene"] += 1
                continue
            gene = gene_index.get(gx)
            if gene is None:
                tally["counted, GX not in gene table"] += 1
                continue
            cells.append(cell)
            genes.append(gene)
    return cells, genes, tally


def rebuild(bam_url, metrics, genes_tsv, out, regions, library_id):
    m = pd.read_csv(metrics, usecols=["barcode", "is_cell"])
    cells = np.sort(m.loc[m["is_cell"] == 1, "barcode"].to_numpy(str))
    cell_index = {b: i for i, b in enumerate(cells)}
    genes = pd.read_csv(genes_tsv, sep="\t")
    gene_index = {g: i for i, g in enumerate(genes["id"])}
    print(f"{len(cells)} called cells, {len(genes)} genes; BAM {bam_url}" + (f", regions {regions}" if regions else ""))

    cell_hits, gene_hits = array("i"), array("i")
    tally = collections.Counter()
    start = time.time()
    with pysam.AlignmentFile(bam_url, "rb") as bam:
        contigs = regions or [sq["SN"] for sq in bam.header.to_dict()["SQ"]]
    # One contig at a time through the index, so a dropped connection costs one contig, not the whole BAM.
    # Counted reads always carry a gene, so they are always on a contig; unmapped reads are never counted.
    for contig in contigs:
        for attempt in range(1, MAX_TRIES + 1):
            try:
                cells_c, genes_c, tally_c = count_contig(bam_url, contig, cell_index, gene_index)
                break
            except OSError as error:
                print(f"  {contig}: attempt {attempt} failed ({error}); retrying", flush=True)
                time.sleep(30 * attempt)
        else:
            raise RuntimeError(f"{contig}: failed {MAX_TRIES} times")
        cell_hits.extend(cells_c)
        gene_hits.extend(genes_c)
        tally.update(tally_c)
        if tally_c["reads seen"]:
            print(f"  {contig}: {tally_c['reads seen']:,} reads, {len(cells_c):,} counted; total {len(cell_hits):,}, "
                  f"{time.time() - start:.0f} s", flush=True)
    tally["counted into the matrix"] = len(cell_hits)
    for k, v in tally.most_common():
        print(f"  {k}: {v:,}")

    counts = sp.csc_matrix((np.ones(len(cell_hits), dtype=np.int32),
                            (np.frombuffer(gene_hits, dtype=np.int32), np.frombuffer(cell_hits, dtype=np.int32))),
                           shape=(len(genes), len(cells)))
    counts.sum_duplicates()
    counts.sort_indices()
    with h5py.File(out, "w") as f:
        f.attrs["filetype"], f.attrs["version"] = "matrix", 2
        f.attrs["software_version"] = "rebuilt from Cell Ranger BAM (multiome/scripts/rebuild_gex_h5.py)"
        f.attrs["library_ids"] = np.array([library_id.encode()])
        g = f.create_group("matrix")
        g["barcodes"] = cells.astype("S")
        g["data"], g["indices"], g["indptr"] = counts.data, counts.indices.astype(np.int64), counts.indptr.astype(np.int64)
        g["shape"] = np.array(counts.shape, dtype=np.int32)
        ft = g.create_group("features")
        for k in ["id", "name", "feature_type", "genome", "interval"]:
            ft[k] = genes[k].to_numpy(str).astype("S")
        ft["_all_tag_keys"] = np.array([b"genome", b"interval"])
    print(f"wrote {out}: {counts.shape[0]} genes x {counts.shape[1]} cells, {counts.sum():,} UMIs, "
          f"{time.time() - start:.0f} s")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--probe", action="store_true")
    p.add_argument("--metrics")
    p.add_argument("--genes")
    p.add_argument("--out")
    p.add_argument("--regions", help="comma-separated contigs (uses the BAM index)")
    p.add_argument("--library-id", default="")
    args = p.parse_args()
    if args.probe:
        probe(args.bam)
    else:
        rebuild(args.bam, args.metrics, args.genes, args.out,
                args.regions.split(",") if args.regions else None, args.library_id)


if __name__ == "__main__":
    main()
