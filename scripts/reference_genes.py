"""Gene table of a 10x reference, in Cell Ranger's feature order, from the reference tarball itself.

Streams the reference tarball once and keeps only reference.json, genes/genes.gtf.gz and regions/{tss,transcripts}.bed
(saved next to --out, about 60 MB). It checks that the
reference's GTF hash matches the hash recorded in a fragments file header, i.e. that this is the reference the
data were processed with. Writes one row per gene, in GTF order:
  id, name, feature_type, genome, interval
interval is Cell Ranger ARC's: chrom:start-end spanning the gene's entries in the reference's regions/tss.bed
(0-based start), or NA for genes with none (the chrM genes). The members are reused if already saved.

  python scripts/reference_genes.py --tarball URL --fragments <atac_fragments.tsv.gz> --out genes.tsv \
      [--check-h5 <a filtered_feature_bc_matrix.h5 made with the same reference>]
"""

import argparse
import gzip
import hashlib
import io
import json
import pathlib
import re
import tarfile
import urllib.request

import h5py
import numpy as np
import pandas as pd


def fragment_header(path):
    """The '# key=value' header lines of a 10x fragments file."""
    header = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                break
            key, _, value = line[2:].strip().partition("=")
            header.setdefault(key, value)
    return header


KEEP = ["reference.json", "genes/genes.gtf.gz", "regions/tss.bed", "regions/transcripts.bed"]


def read_tarball(url, keep_dir):
    """The KEEP members, streamed out of the tarball and saved under keep_dir."""
    saved = {name: next(pathlib.Path(keep_dir).glob(f"refdata-*/{name}"), None) for name in KEEP}
    if all(saved.values()):
        print(f"  reusing saved members under {keep_dir}")
        reference = json.loads(saved["reference.json"].read_text())
        return reference, saved["genes/genes.gtf.gz"].read_bytes(), saved["regions/tss.bed"].read_text()
    wanted = dict.fromkeys(KEEP)
    request = urllib.request.Request(url, headers={"User-Agent": "curl/8.0"})  # the 10x CDN rejects urllib's default
    with urllib.request.urlopen(request) as response, tarfile.open(fileobj=response, mode="r|gz") as tar:
        for member in tar:
            name = member.name.split("/", 1)[-1]
            if name in wanted:
                wanted[name] = tar.extractfile(member).read()
                path = pathlib.Path(keep_dir) / member.name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(wanted[name])
                print(f"  extracted {member.name} ({member.size:,} bytes) to {path}", flush=True)
            if all(v is not None for v in wanted.values()):
                break
    assert all(v is not None for v in wanted.values()), f"missing from tarball: {wanted.keys()}"
    return json.loads(wanted["reference.json"]), wanted["genes/genes.gtf.gz"], wanted["regions/tss.bed"].decode()


def gene_table(gtf_gz: bytes, tss_bed: str) -> pd.DataFrame:
    """One row per gene in GTF order; interval spans the gene's entries in regions/tss.bed (NA if none)."""
    genes = {}
    attr = re.compile(r'(\w+) "([^"]*)"')
    with gzip.open(io.BytesIO(gtf_gz), "rt") as f:
        for line in f:
            fields = line.split("\t")
            if line.startswith("#") or fields[2] != "gene":
                continue
            a = dict(attr.findall(fields[8]))
            genes[a["gene_id"]] = a.get("gene_name", a["gene_id"])
    tss = pd.read_csv(io.StringIO(tss_bed), sep="\t", header=None, names=["chrom", "start", "end", "id", "score", "strand"])
    span = tss.groupby("id").agg(chrom=("chrom", "first"), start=("start", "min"), end=("end", "max"))
    interval = (span["chrom"] + ":" + span["start"].astype(str) + "-" + span["end"].astype(str)).to_dict()
    rows = [(gene_id, name, "Gene Expression", "GRCh38", interval.get(gene_id, "NA")) for gene_id, name in genes.items()]
    return pd.DataFrame(rows, columns=["id", "name", "feature_type", "genome", "interval"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tarball", required=True)
    p.add_argument("--fragments", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--check-h5")
    args = p.parse_args()

    header = fragment_header(args.fragments)
    print(f"fragments header: {header.get('pipeline_version')}, reference {header.get('reference_version')}, "
          f"gtf hash {header.get('reference_gtf_hash')}")
    reference, gtf_gz, tss_bed = read_tarball(args.tarball, pathlib.Path(args.out).parent)
    print("reference.json:", {k: v for k, v in reference.items() if not isinstance(v, (list, dict)) or k == "genomes"})
    assert reference["fasta_hash"] == header["reference_fasta_hash"], "FASTA hash differs from the fragments header"
    # The GTF hash is recorded under a key naming its format (e.g. "gtf_hash.gz"); also hash the extracted file.
    gtf_hashes = {k: v for k, v in reference.items() if k.startswith("gtf_hash")}
    gtf_hashes["sha1(genes.gtf.gz)"] = hashlib.sha1(gtf_gz).hexdigest()
    gtf_hashes["sha1(genes.gtf)"] = hashlib.sha1(gzip.decompress(gtf_gz)).hexdigest()
    print("GTF hashes:", gtf_hashes)
    assert header["reference_gtf_hash"] in gtf_hashes.values(), "GTF hash differs from the fragments header"
    print("reference hashes match the fragments header")

    genes = gene_table(gtf_gz, tss_bed)
    genes.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}: {len(genes)} genes")

    if args.check_h5:
        with h5py.File(args.check_h5, "r") as f:
            ft = f["matrix"]["features"]
            gex = ft["feature_type"][:] == b"Gene Expression"
            h5 = pd.DataFrame({k: ft[k][:][gex].astype(str) for k in ["id", "name", "feature_type", "genome", "interval"]})
        same = {k: bool((h5[k].to_numpy() == genes[k].to_numpy()).all()) if len(h5) == len(genes) else False
                for k in h5.columns}
        print(f"check against {args.check_h5}: {len(h5)} vs {len(genes)} genes; columns equal: {same}")
        if not all(same.values()) and len(h5) == len(genes):
            for k in [k for k, ok in same.items() if not ok]:
                bad = np.flatnonzero(h5[k].to_numpy() != genes[k].to_numpy())
                print(f"  {k}: {len(bad)} differ, e.g.", list(zip(h5[k].iloc[bad[:5]], genes[k].iloc[bad[:5]])))
        assert all(same.values()), "gene table does not match the h5 features"
        print("gene table matches the h5 features exactly")


if __name__ == "__main__":
    main()
