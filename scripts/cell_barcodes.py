"""Write the Cell Ranger-called barcodes of a filtered_feature_bc_matrix.h5 as AMULET's single-cell CSV.

  python scripts/cell_barcodes.py --h5 x_filtered_feature_bc_matrix.h5 --out singlecell.csv
"""

import argparse

import h5py


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    with h5py.File(args.h5) as f:
        barcodes = [b.decode() for b in f["matrix/barcodes"][:]]
    with open(args.out, "w") as out:
        out.write("barcode,is__cell_barcode\n")
        out.writelines(f"{b},1\n" for b in barcodes)


if __name__ == "__main__":
    main()
