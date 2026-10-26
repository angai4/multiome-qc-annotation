"""Combine per-sample 10x multiome outputs into one dataset, as `cellranger-arc aggr` lays them out (no depth
normalisation): sample k's barcodes get the GEM-well suffix -k, in the order given.

Writes <out>_filtered_feature_bc_matrix.h5 (the samples' GEX matrices side by side; genes must match) and
<out>_atac_fragments.tsv.gz (+ .tbi): all fragments with renamed barcodes, sorted by position, bgzipped and indexed.
The samples' fragment files must be sorted and tabix-indexed (Cell Ranger's are). They are merged one contig at a
time through their indexes, so no temporary files are written.

  python scripts/aggregate_samples.py --out data/raw/bmmc/bmmc_s1d1_s1d3 \
      --sample s1d1 h5 fragments --sample s1d3 h5 fragments
"""

import argparse
import heapq

import h5py
import numpy as np
import pysam
import scipy.sparse as sp


def rename(barcodes, well):
    assert all(b.endswith("-1") for b in barcodes), "expected single-library barcodes ending in -1"
    return np.array([b[:-2] + f"-{well}" for b in barcodes])


def aggregate_matrices(samples, out):
    blocks, barcodes, features = [], [], None
    for well, (name, h5, _) in enumerate(samples, 1):
        with h5py.File(h5, "r") as f:
            g = f["matrix"]
            blocks.append(sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:])))
            barcodes.append(rename(g["barcodes"][:].astype(str), well))
            ft = {k: g["features"][k][:] for k in g["features"]}
        if features is None:
            features = ft
        assert all(np.array_equal(features[k], ft[k]) for k in features), f"{name}: features differ"
        print(f"  {name}: {blocks[-1].shape[1]} cells -> barcode suffix -{well}")
    x = sp.hstack(blocks).tocsc()
    with h5py.File(f"{out}_filtered_feature_bc_matrix.h5", "w") as f:
        f.attrs["filetype"], f.attrs["version"] = "matrix", 2
        f.attrs["software_version"] = "aggregated (multiome/scripts/aggregate_samples.py)"
        f.attrs["library_ids"] = np.array([name.encode() for name, _, _ in samples])
        g = f.create_group("matrix")
        g["barcodes"] = np.concatenate(barcodes).astype("S")
        g["data"], g["indices"], g["indptr"] = x.data, x.indices.astype(np.int64), x.indptr.astype(np.int64)
        g["shape"] = np.array(x.shape, dtype=np.int32)
        ft = g.create_group("features")
        for k, v in features.items():
            ft[k] = v
    print(f"wrote {out}_filtered_feature_bc_matrix.h5: {x.shape[0]} genes x {x.shape[1]} cells")


def aggregate_fragments(samples, out):
    tabix = [pysam.TabixFile(fragments) for _, _, fragments in samples]
    contigs = list(tabix[0].contigs) + [c for t in tabix[1:] for c in t.contigs if c not in tabix[0].contigs]
    contigs = list(dict.fromkeys(contigs))
    target = f"{out}_atac_fragments.tsv.gz"
    counts = [0] * len(samples)
    with pysam.BGZFile(target, "w") as w:
        w.write("".join(f"{line}\n" for line in tabix[0].header).encode())
        for contig in contigs:
            def renamed(well, t):
                for line in t.fetch(contig):
                    chrom, start, end, barcode, count = line.split("\t")
                    assert barcode.endswith("-1")
                    counts[well - 1] += 1
                    yield int(start), int(end), f"{chrom}\t{start}\t{end}\t{barcode[:-2]}-{well}\t{count}\n"
            streams = [renamed(well, t) for well, t in enumerate(tabix, 1) if contig in t.contigs]
            batch = []
            for _, _, line in heapq.merge(*streams):
                batch.append(line)
                if len(batch) == 100_000:  # write in batches, never a whole contig at once
                    w.write("".join(batch).encode())
                    batch = []
            w.write("".join(batch).encode())
            print(f"  {contig} done; fragments so far {sum(counts):,}", flush=True)
    for (name, _, _), n in zip(samples, counts):
        print(f"  {name}: {n:,} fragments")
    pysam.tabix_index(target, preset="bed", force=True)
    print(f"wrote {target} (+ .tbi): {sum(counts):,} fragments")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--sample", nargs=3, action="append", metavar=("NAME", "H5", "FRAGMENTS"), required=True)
    args = p.parse_args()
    aggregate_matrices(args.sample, args.out)
    aggregate_fragments(args.sample, args.out)


if __name__ == "__main__":
    main()
