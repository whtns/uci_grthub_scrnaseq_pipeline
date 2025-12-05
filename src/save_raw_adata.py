#!/usr/bin/env python3

import argparse
from pathlib import Path
import scanpy as sc
import sys
import grthub_tools as gt


def main():
	parser = argparse.ArgumentParser(description="Read an AnnData file and save it (works with Snakemake fallback)")
	parser.add_argument("input", nargs="?", help="Path to input .h5ad file")
	parser.add_argument("raw_input", nargs="?", help="Path to raw input .h5ad file")
	parser.add_argument("data_output", nargs="?", help="Path to output .h5ad file")
	parser.add_argument("scaled_data_output", nargs="?", help="Path to output .h5ad file")
	args = parser.parse_args()

	in_path = args.input
	raw_in_path = args.input
	out_path = args.data_output

	# If running under Snakemake, allow using snakemake.input/output when args omitted
	if (in_path is None) or (out_path is None):
		try:
			# 'snakemake' is injected when executing via Snakemake
			snakemake  # type: ignore
			in_path = in_path or snakemake.input[0]
			raw_in_path = raw_in_path or snakemake.input[1]
			out_path = out_path or snakemake.output[0]
		except NameError:
			pass

	if in_path is None or raw_in_path is None or out_path is None:
		parser.error("input and output paths required (or run via Snakemake with input/output)")

	in_path = Path(in_path)
	raw_in_path = Path(raw_in_path)
	out_path = Path(out_path)
	if not in_path.exists():
		print(f"ERROR: input file not found: {in_path}", file=sys.stderr)
		sys.exit(2)
	adata = sc.read_h5ad(str(in_path))
	adata.var_names_make_unique()
	adata.obs_names_make_unique()
	raw_adata = sc.read_h5ad(str(raw_in_path))
	raw_adata = raw_adata.raw.copy().to_adata()
	gt.write_adata_to_csv(raw_adata, in_path.stem)

	raw_adata.var_names_make_unique()
	raw_adata.obs_names_make_unique()

	# sc.pp.filter_genes(raw_adata, min_cells=round(raw_adata.shape[0]*0.0005))

	raw_adata.layers['counts'] = raw_adata.copy()
	raw_adata = gt.preprocess_adata(raw_adata, batch = True)
	raw_adata.obsm = adata.obsm.copy()
	raw_adata.uns = adata.uns.copy()
	raw_adata.obs = adata.obs.copy()
	raw_adata.layers['data'] = adata.X.copy()

	raw_adata.write_h5ad(str(args.data_output))
	print(f"Wrote AnnData to: {args.data_output}")
	sc.pp.scale(raw_adata, max_value=10)
	raw_adata.layers['scale.data'] = adata.X.copy()
	raw_adata.write_h5ad(str(args.scaled_data_output))


if __name__ == "__main__":
	main()
