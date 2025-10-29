#!/usr/bin/env python
import argparse
import scanpy as sc
import os
import sys
import warnings
warnings.filterwarnings("ignore")


def main():
	parser = argparse.ArgumentParser(description="Load CellRanger filtered HDF5 and add QC gene flags to adata.var")
	parser.add_argument("--adata_path", default="filtered_feature_bc_matrix.h5",
						help="Path to the 10x HDF5 file (default: filtered_feature_bc_matrix.h5)")
	parser.add_argument("--output", default=None,
						help="Optional path to write the modified AnnData (h5ad). If omitted, the script will not save.")
	parser.add_argument("--min_genes", type=int, default=200,
					help="Minimum number of genes required per cell (default: 200)")
	parser.add_argument("--min_cells", type=int, default=3,
					help="Minimum number of cells required per gene (default: 3)")
	parser.add_argument("--organism", choices=["human", "mouse", "auto"], default="auto",
					help="Organism to determine mitochondrial gene prefix: 'human' -> MT-, 'mouse' -> Mt-, or 'auto' to guess from gene names (default: auto)")
	parser.add_argument("--mt_thresh", type=float, default=5.0,
					help="Maximum percent mitochondrial counts allowed per cell (default: 5.0)")
	args = parser.parse_args()

	if not os.path.exists(args.adata_path):
		print(f"Error: Specified adata_path not found: {args.adata_path}", file=sys.stderr)
		sys.exit(2)

	try:
		adata = sc.read_10x_h5(args.adata_path)
		adata.var_names_make_unique()
	except Exception as e:
		print(f"Error reading adata_path {args.adata_path}: {e}", file=sys.stderr)
		sys.exit(3)

	# mitochondrial genes detection
	if args.organism == 'human':
		mt_prefix = 'MT-'
	elif args.organism == 'mouse':
		mt_prefix = 'Mt-'
	else:
		# auto-detect: if any gene starts with 'MT-' assume human, else if 'Mt-' assume mouse, else default to 'MT-'
		if any(gn.startswith('MT-') for gn in adata.var_names):
			mt_prefix = 'MT-'
		elif any(gn.startswith('Mt-') for gn in adata.var_names):
			mt_prefix = 'Mt-'
		else:
			mt_prefix = 'MT-'

	adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
	# ribosomal genes
	adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
	# hemoglobin genes
	adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

	# Calculate common QC metrics (adds e.g. 'pct_counts_mt' to adata.obs)
	sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True)

	# Apply configurable filters
	orig_n_cells, orig_n_genes = adata.n_obs, adata.n_vars
	sc.pp.filter_cells(adata, min_genes=args.min_genes)
	sc.pp.filter_genes(adata, min_cells=args.min_cells)

	# Optional mitochondrial filtering (hard-coded threshold here); keep rows for later use
	if 'pct_counts_mt' not in adata.obs:
		# calculate_qc_metrics should have added this, but guard just in case
		sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

	# Apply mitochondrial pct filter using provided threshold
	adata = adata[adata.obs.pct_counts_mt < args.mt_thresh, :]

	new_n_cells, new_n_genes = adata.n_obs, adata.n_vars
	# print(f"Filtered cells: {orig_n_cells} -> {new_n_cells} (min_genes={args.min_genes}, mt<{args.mt_thresh}%)")
	# print(f"Filtered genes: {orig_n_genes} -> {new_n_genes} (min_cells={args.min_cells})")

	if args.output:
		outdir = os.path.dirname(args.output)
		if outdir:
			os.makedirs(outdir, exist_ok=True)
		try:
			adata.write(args.output)
			print(f"Wrote annotated/filtered AnnData to {args.output}", file=sys.stderr)
		except Exception as e:
			print(f"Warning: failed to write output {args.output}: {e}", file=sys.stderr)

	# Return number of cells remaining (callers should print this to stdout if desired)
	return new_n_cells


if __name__ == '__main__':
	# Print only the integer result to stdout. All diagnostics go to stderr.
	result = main()
	out_int = None
	# Try direct int conversion
	try:
		out_int = int(result)
	except Exception:
		# If result is a string that contains digits, try to extract
		try:
			s = str(result).strip()
			# Remove any non-digit (and non-dot) characters, attempt float->int
			import re
			m = re.search(r"(-?\d+)", s)
			if m:
				out_int = int(m.group(1))
			else:
				# Try float conversion
				out_int = int(float(s))
		except Exception:
			out_int = int(result)
			print(f"Warning: unexpected result from main(), falling back to 0. original: {result}", file=sys.stderr)

	print(out_int)