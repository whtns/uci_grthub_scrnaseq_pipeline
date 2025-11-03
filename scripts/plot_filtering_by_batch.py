#!/usr/bin/env python3
"""
Plot number of genes detected per cell grouped by batch at key filtering steps.

Stages shown:
 - raw: original AnnData
 - after_filter_cells: after sc.pp.filter_cells(min_genes=...)
 - after_filter_cells_and_genes: after sc.pp.filter_genes(min_cells=...) applied after filter_cells

Usage:
  python scripts/plot_filtering_by_batch.py --input path/to/adata.h5ad --batch-key batch --min-genes 200 --min-cells 3 --out plot.png

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
        batch = adata.obs.loc[s.index, batch_key].astype(str)
    else:
        batch = pd.Series(["unknown"] * len(s), index=s.index)
    df = pd.DataFrame({"batch": batch.values, "n_genes": s.values})
    df["stage"] = stage_label
    return df


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Plot genes detected (n_genes) per cell by batch at filtering steps"
    )
    p.add_argument("--input", "-i", required=True, help="Path to input AnnData (.h5ad)")
    p.add_argument("--batch-key", default="batch", help="Column in adata.obs to use as batch (default: 'batch')")
    p.add_argument("--min-genes", type=int, default=200, help="min_genes for sc.pp.filter_cells (default: 200)")
    p.add_argument("--min-cells", type=int, default=3, help="min_cells for sc.pp.filter_genes (default: 3)")
    p.add_argument("--out", "-o", default="filtering_by_batch.png", help="Output figure path (png/svg/pdf)")
    # use violin plots by default; allow opting out with --no-violin
    p.add_argument("--no-violin", action="store_true", help="Use boxplot instead of violin (default: violin)")
    args = p.parse_args(argv)

    input_path = Path(args.input)
    if not input_path.exists():
        p.error(f"Input file not found: {input_path}")

    # Read AnnData (scanpy.read handles .h5ad)
    try:
        adata = sc.read(input_path)
    except Exception as e:
        # fallback to anndata read
        adata = ad.read_h5ad(input_path)

    # Ensure obs index is set
    if adata.obs_names is None or len(adata.obs_names) == 0:
        raise ValueError("AnnData contains no observations (cells)")

    # Stage 1: raw
    df_raw = build_stage_df(adata, args.batch_key, "raw")

    # Stage 2: after filter_cells
    ad_fc = adata.copy()
    sc.pp.filter_cells(ad_fc, min_genes=args.min_genes)
    df_fc = build_stage_df(ad_fc, args.batch_key, "after_filter_cells")

    # Stage 3: after filter_cells + filter_genes
    ad_fcg = ad_fc.copy()
    sc.pp.filter_genes(ad_fcg, min_cells=args.min_cells)
    df_fcg = build_stage_df(ad_fcg, args.batch_key, "after_filter_cells_and_genes")

    df_all = pd.concat([df_raw, df_fc, df_fcg], axis=0)

    # Order stages for plotting
    stage_order = ["raw", "after_filter_cells", "after_filter_cells_and_genes"]
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
    plt.tight_layout()
    out_path = Path(args.out)
    plt.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    main()
