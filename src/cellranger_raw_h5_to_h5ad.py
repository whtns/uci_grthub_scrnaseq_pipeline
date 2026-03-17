#!/usr/bin/env python
import argparse
import os
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert CellRanger raw_feature_bc_matrix.h5 to h5ad for CellSweep"
    )
    parser.add_argument("--input", required=True, help="Path to raw_feature_bc_matrix.h5")
    parser.add_argument("--output", required=True, help="Path to output .h5ad")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    adata = sc.read_10x_h5(args.input)
    adata.var_names_make_unique()

    # CellSweep expects a celltype column in obs.
    adata.obs["celltype"] = args.sample

    outdir = os.path.dirname(args.output)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    adata.write(args.output)


if __name__ == "__main__":
    main()
