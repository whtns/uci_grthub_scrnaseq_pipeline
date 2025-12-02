#!/usr/bin/env python3
"""
Plot number of genes detected per cell grouped by batch at key filtering steps.

Stages shown:
 - raw: original AnnData
 - after_filter_cells: after sc.pp.filter_cells(min_genes=...)
 - after_filter_cells_and_genes: after sc.pp.filter_genes(min_cells=...) applied after filter_cells

Usage:
  python src/plot_filtering_timeline.py --input path/to/adata.h5ad --batch-key batch --min-genes 200 --min-cells 3 --out plot.png

The script saves a PNG (or other matplotlib-supported format) with a boxplot/violin showing n_genes per cell by batch for each stage.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse


def n_genes_per_cell(adata: ad.AnnData) -> pd.Series:
    """Return a Series indexed by adata.obs_names giving the number of genes detected (>0) per cell."""
    X = adata.X
    if X is None:
        # no matrix
        return pd.Series(index=adata.obs_names, dtype=int)
    if sparse.issparse(X):
        # (n_obs, n_vars) sparse matrix
        vals = (X > 0).sum(axis=1)
        # vals is matrix-like; convert
        try:
            arr = np.asarray(vals).reshape(-1)
        except Exception:
            arr = np.array(vals).ravel()
    else:
        arr = (X > 0).sum(axis=1)
        arr = np.asarray(arr).reshape(-1)
    return pd.Series(arr, index=adata.obs_names, name="n_genes")


def build_stage_df(adata: ad.AnnData, batch_key: str, stage_label: str) -> pd.DataFrame:
    s = n_genes_per_cell(adata)
    # get batch values; if missing, use "unknown"
    if batch_key in adata.obs:
        # Align to the cell index returned by n_genes_per_cell. Reindex first
        # (to get the same index), then coerce to a non-categorical/object
        # dtype before filling missing values with 'unknown'. Doing astype
        # before fillna avoids errors when the column is categorical.
        batch = adata.obs[batch_key].reindex(s.index)
        batch = batch.astype(object).fillna('unknown').astype(str)
    else:
        batch = pd.Series(["unknown"] * len(s), index=s.index)
    # construct DataFrame using the aligned index to ensure both arrays have
    # identical length and order
    df = pd.DataFrame({"batch": batch.values, "n_genes": s.values}, index=s.index)
    df["stage"] = stage_label
    return df


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Plot genes detected (n_genes) per cell by batch at filtering steps"
    )
    p.add_argument("--input", "-i", required=True, nargs='+', help="One or more paths to input AnnData (.h5ad) or 10x h5 files")
    p.add_argument("--batch-key", default="batch", help="Column in adata.obs to use as batch (default: 'batch')")
    p.add_argument("--batch-value", default=None,
                   help="If the batch column is missing, create it with this constant value for all cells (optional).")
    p.add_argument("--min-genes", type=int, default=200, help="min_genes for sc.pp.filter_cells (default: 200)")
    p.add_argument("--min-cells", type=int, default=3, help="min_cells for sc.pp.filter_genes (default: 3)")
    p.add_argument("--max-genes", type=int, default=8000, help="Maximum genes detected to retain a cell (default: 8000)")
    p.add_argument("--max-pct-mito", type=float, default=5.0, help="Maximum percent mitochondrial counts to retain a cell (default: 5.0)")
    p.add_argument("--out", "-o", default="output/qc/filtering_timeline.png", help="Output figure path (png/svg/pdf)")
    # use violin plots by default; allow opting out with --no-violin
    p.add_argument("--no-violin", action="store_true", help="Use boxplot instead of violin (default: violin)")
    p.add_argument("--plot-doublets", action="store_true", help="Also produce a doublet-summary plot (requires doublet data in adata.obs or --doublet-csv)")
    p.add_argument("--doublet-csv", nargs='*', default=None, help="One or more Scrublet CSV files (barcode,doublet_score,predicted_doublet) to merge into adata.obs before plotting")
    args = p.parse_args(argv)

    input_paths = [Path(p) for p in args.input]
    for pth in input_paths:
        if not pth.exists():
            p.error(f"Input file not found: {pth}")

    def _read_one(path: Path) -> ad.AnnData:
        # Read AnnData. Support both AnnData (.h5ad) and 10x matrix HDF5 (.h5)
        # If a matrix .h5 (CellRanger filtered_feature_bc_matrix.h5) is provided,
        # use scanpy.read_10x_h5 which returns an AnnData. Otherwise try scanpy/anndata
        if path.suffix == ".h5":
            try:
                return sc.read_10x_h5(str(path))
            except Exception:
                try:
                    return sc.read(path)
                except Exception:
                    return ad.read_h5ad(path)
        else:
            try:
                return sc.read(path)
            except Exception:
                return ad.read_h5ad(path)

    adatas = [_read_one(p) for p in input_paths]

    for adata in adatas:
        adata.var_names_make_unique()

    if len(adatas) == 1:
        adata = adatas[0]
    else:
        # Concatenate multiple AnnData objects into one, creating the batch
        # observation column named by args.batch_key. Use the input filenames
        # (stem) as keys for the batch labels. If the user provided
        # --batch-value, ignore it in multi-input mode (the per-file sample
        # names are usually what the user wants); warn the user.
        keys = [p.parent.parent.stem for p in input_paths]
        if args.batch_value is not None:
            print(f"Warning: --batch-value ignored when concatenating multiple inputs; using file-based sample keys: {keys}")
        # anndata.concat (aliased as ad.concat) will add a column named
        # according to `label` with the provided keys
        adata = ad.concat(adatas, join='outer', label=args.batch_key, keys=keys)

    # Ensure obs index is set
    if adata.obs_names is None or len(adata.obs_names) == 0:
        raise ValueError("AnnData contains no observations (cells)")

    # If the requested batch key is missing and a batch value was provided, add it
    if args.batch_key not in adata.obs:
        if args.batch_value is not None:
            # assign the provided batch value (string) to all cells
            adata.obs[args.batch_key] = str(args.batch_value)
        else:
            # leave it missing; build_stage_df will use 'unknown'
            pass

    # Compute mitochondrial percent and n_genes per cell
    var_names = adata.var_names.astype(str)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Stage 1: raw
    df_raw = build_stage_df(adata, args.batch_key, "raw")

    # Stage 2: apply QC filters requested by the user: retain cells with
    # fewer than `--max-genes` genes and mitochondrial percent less than
    # `--max-pct-mito`.
    s_raw = n_genes_per_cell(adata)
    if 'pct_counts_mt' not in adata.obs:
        # Ensure the column exists to avoid KeyError; default to 0
        adata.obs['pct_counts_mt'] = 0.0

    # Align pct_counts_mt to the same index/order as s_raw to avoid boolean
    # operation dimension/index mismatch. Fill missing values conservatively
    # with 0.0 (no mitochondrial reads) and ensure numeric dtype.
    pct_mt = adata.obs.get('pct_counts_mt')
    if pct_mt is None:
        pct_mt = pd.Series(0.0, index=adata.obs_names)
    else:
        pct_mt = pct_mt.reindex(s_raw.index).fillna(0.0).astype(float)

    mask_keep = (s_raw.reindex(pct_mt.index).fillna(0) < args.max_genes) & (pct_mt < float(args.max_pct_mito))
    # Subset AnnData with a boolean mask: use adata[mask] (AnnData does not
    # implement .loc like pandas.DataFrame)
    ad_qc = adata[mask_keep.values].copy()
    df_qc = build_stage_df(ad_qc, args.batch_key, "after_qc_filters")

    # Stage 3: after filter_genes applied to the filtered cells
    ad_fcg = ad_qc.copy()
    sc.pp.filter_genes(ad_fcg, min_cells=args.min_cells)
    df_fcg = build_stage_df(ad_fcg, args.batch_key, "after_filter_cells_and_genes")

    df_all = pd.concat([df_raw, df_qc, df_fcg], axis=0)

    # Order stages for plotting
    stage_order = ["raw", "after_qc_filters", "after_filter_cells_and_genes"]
    df_all["stage"] = pd.Categorical(df_all["stage"], categories=stage_order, ordered=True)

    # Prepare the plot: we want genes detected (y) and batch on x, with separate facet/hue per stage.
    sns.set(style="whitegrid")
    n_stages = len(stage_order)
    fig, axes = plt.subplots(1, n_stages, figsize=(5 * n_stages, 5), sharey=True)

    for ax, stage in zip(axes, stage_order):
        d = df_all[df_all.stage == stage]
        if not args.no_violin:
            sns.violinplot(x="batch", y="n_genes", data=d, inner="quartile", ax=ax)
        else:
            sns.boxplot(x="batch", y="n_genes", data=d, ax=ax)
            sns.stripplot(x="batch", y="n_genes", data=d, color="0.3", size=2, jitter=True, ax=ax)

        # annotate number of cells per batch above each violin/box
        # get tick positions and labels (they align)
        xticks = ax.get_xticks()
        xticklabels = [t.get_text() for t in ax.get_xticklabels()]
        # compute y-range and offset
        y_min, y_max = ax.get_ylim()
        if y_max > y_min:
            offset = (y_max - y_min) * 0.06
        else:
            offset = max(1.0, y_max * 0.05)

        for x_pos, lbl in zip(xticks, xticklabels):
            # if label is empty (rare), skip
            if lbl == "":
                continue
            cnt = int((d["batch"] == lbl).sum())
            y_text = y_max + offset
            ax.text(x_pos, y_text, f"n={cnt}", ha="center", va="bottom", fontsize=9)

        # extend y-limits so annotations are visible
        ax.set_ylim(y_min, y_max + 2 * offset)

        ax.set_title(stage)
        ax.set_xlabel(args.batch_key)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    axes[0].set_ylabel("n_genes (genes detected per cell)")
    # Add text showing the filtering thresholds used
    info_text = (
        f"min_genes={args.min_genes}, min_cells={args.min_cells}, "
        f"max_genes={args.max_genes}, max_pct_mito={args.max_pct_mito}"
    )
    # Put as a suptitle and make room for it
    fig.suptitle(info_text, fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    out_path = Path(args.out)
    plt.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")

    # Optionally plot doublet summaries
    if args.plot_doublets:
        try:
            out_doublet = out_path.with_name(out_path.stem + "_doublets" + out_path.suffix)
        except Exception:
            out_doublet = Path(str(out_path) + "_doublets.png")
        plot_doublets(adata, args.batch_key, doublet_csvs=args.doublet_csv, out_path=out_doublet)


def plot_doublets(adata: ad.AnnData, batch_key: str, doublet_csvs: list | None = None, out_path: Path | str | None = None):
    """Plot doublet score distributions per batch and percent predicted doublets per batch.

    The function looks for columns `doublet_score` and `predicted_doublet` in `adata.obs`.
    If `doublet_csvs` is provided, it will attempt to read one or more Scrublet CSV files
    (with columns `barcode`, `doublet_score`, `predicted_doublet`) and merge by barcode.
    """
    # Prepare a DataFrame with batch, doublet_score, predicted_doublet
    obs = adata.obs.copy()

    # If CSVs provided, read and merge into obs (by barcode)
    if doublet_csvs:
        csvs = []
        for p in doublet_csvs:
            try:
                df = pd.read_csv(p)
            except Exception:
                # try Path conversion
                df = pd.read_csv(str(p))
            # ensure expected columns
            if 'barcode' not in df.columns:
                continue
            # keep only relevant columns
            keep_cols = [c for c in ['barcode', 'doublet_score', 'predicted_doublet'] if c in df.columns]
            df = df[keep_cols].drop_duplicates(subset=['barcode'])
            csvs.append(df.set_index('barcode'))
        if len(csvs) > 0:
            df_csv = pd.concat(csvs, axis=0)
            # join to obs by index (barcode). Many AnnData objects store barcodes as obs_names
            # If obs index contains full barcode strings, merge directly; otherwise try to match prefix/suffix.
            try:
                merged = obs.join(df_csv, how='left')
            except Exception:
                # fallback: align by barcode column if exists
                if 'barcode' in obs.columns:
                    merged = obs.set_index('barcode', drop=False).join(df_csv, how='left')
                else:
                    merged = obs.join(df_csv, how='left')
            obs = merged

    # Check columns exist
    if 'doublet_score' not in obs.columns and 'predicted_doublet' not in obs.columns:
        print('No doublet information found in AnnData.obs or provided CSVs. Skipping doublet plot.')
        return

    # Ensure batch column exists
    if batch_key not in obs.columns:
        obs[batch_key] = 'unknown'

    # Normalize predicted_doublet to boolean if present
    if 'predicted_doublet' in obs.columns:
        try:
            pred = obs['predicted_doublet'].astype(bool)
        except Exception:
            pred = obs['predicted_doublet'].fillna(False).astype(bool)
        obs['predicted_doublet'] = pred

    # Build DataFrame for plotting
    plot_df = pd.DataFrame({
        'batch': obs[batch_key].astype(str),
    }, index=obs.index)

    if 'doublet_score' in obs.columns:
        plot_df['doublet_score'] = pd.to_numeric(obs['doublet_score'], errors='coerce')
    if 'predicted_doublet' in obs.columns:
        plot_df['predicted_doublet'] = obs['predicted_doublet'].astype(bool)

    # Drop rows with no data at all
    if 'doublet_score' in plot_df.columns:
        plot_df = plot_df.dropna(subset=['doublet_score'], how='all')

    if len(plot_df) == 0:
        print('No doublet-score rows available after merging/cleanup. Skipping doublet plot.')
        return

    sns.set(style='whitegrid')
    # Create two panels: left = score distribution per batch; right = percent predicted doublet per batch
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax0 = axes[0]
    if 'doublet_score' in plot_df.columns:
        sns.violinplot(x='batch', y='doublet_score', data=plot_df, inner='quartile', ax=ax0)
        ax0.set_ylabel('doublet_score')
    else:
        ax0.text(0.5, 0.5, 'No doublet_score available', ha='center', va='center')
        ax0.set_ylabel('doublet_score')
    ax0.set_xlabel(batch_key)
    ax0.set_xticklabels(ax0.get_xticklabels(), rotation=45, ha='right')

    ax1 = axes[1]
    if 'predicted_doublet' in plot_df.columns:
        pct = plot_df.groupby('batch')['predicted_doublet'].mean().reset_index()
        pct['percent'] = pct['predicted_doublet'] * 100.0
        sns.barplot(x='batch', y='percent', data=pct, ax=ax1)
        ax1.set_ylabel('% predicted doublets')
        for i, row in pct.iterrows():
            ax1.text(i, row['percent'] + 1.0, f"{row['percent']:.1f}%", ha='center')
    else:
        ax1.text(0.5, 0.5, 'No predicted_doublet available', ha='center', va='center')
        ax1.set_ylabel('% predicted doublets')
    ax1.set_xlabel(batch_key)
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')

    plt.tight_layout()
    if out_path is not None:
        outp = Path(out_path)
        outp.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outp, dpi=150)
        print(f"Saved doublet plot to {outp}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
