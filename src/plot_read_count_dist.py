#!/usr/bin/env python3
"""
Compute per-cell total counts from CellRanger filtered_feature_bc_matrix.h5
and plot distribution across samples using Scanpy (preferred).

Saves:
 - output/plots/read_count_distribution.png
 - output/plots/read_count_distribution_summary.csv

This script prefers Scanpy's `read_10x_h5`. If Scanpy isn't available it will
try to read the HDF5 directly using h5py/scipy as a fallback.
"""

from pathlib import Path
from typing import Tuple, List, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

try:
    import scanpy as sc  # type: ignore
except Exception:
    sc = None


def read_counts_from_h5(h5_path: Path) -> Tuple[np.ndarray, np.ndarray, Optional[List[str]]]:
    """Read per-cell counts and genes-detected from a 10x `filtered_feature_bc_matrix.h5` file.

    Tries Scanpy's reader first (recommended). If Scanpy isn't available or
    fails, falls back to reading the HDF5 structure directly using h5py.
    Returns (counts_array, genes_detected_array, barcodes_or_none).
    """
    # Try scanpy when available
    if sc is not None:
        try:
            adata = sc.read_10x_h5(str(h5_path))
            # adata.X is cells x genes; sum over genes for per-cell counts
            counts = np.array(adata.X.sum(axis=1)).ravel()
            # genes detected per cell (non-zero genes)
            try:
                genes = np.array((adata.X > 0).sum(axis=1)).ravel()
            except Exception:
                genes = np.array((adata.X.astype(bool)).sum(axis=1)).ravel()
            barcodes = list(adata.obs_names) if adata.obs_names is not None else None
            return counts, genes, barcodes
        except Exception as e:
            print(f"scanpy failed to read {h5_path}: {e}; falling back to h5py reader")

    # Fallback: use h5py + scipy to reconstruct matrix
    try:
        import h5py
        import scipy.sparse as sp
    except Exception as e:
        raise ImportError(
            "Neither scanpy nor h5py/scipy are available to read the 10x HDF5 file."
        ) from e

    with h5py.File(h5_path, 'r') as f:
        if 'matrix' in f:
            grp = f['matrix']
            data = grp['data'][:]
            indices = grp['indices'][:]
            indptr = grp['indptr'][:]
            shape = tuple(grp['shape'][:])
            mat = sp.csr_matrix((data, indices, indptr), shape=shape)
            # 10x HDF5 stores matrix as features x barcodes (genes x cells)
            counts = np.array(mat.sum(axis=0)).ravel()
            genes = np.array((mat > 0).sum(axis=0)).ravel()
            barcodes = None
            if 'barcodes' in grp:
                barcodes = [b.decode('utf-8') if isinstance(b, bytes) else str(b) for b in grp['barcodes'][:]]
            return counts, genes, barcodes
        else:
            if all(k in f for k in ('data', 'indices', 'indptr', 'shape')):
                data = f['data'][:]
                indices = f['indices'][:]
                indptr = f['indptr'][:]
                shape = tuple(f['shape'][:])
                mat = sp.csr_matrix((data, indices, indptr), shape=shape)
                counts = np.array(mat.sum(axis=0)).ravel()
                genes = np.array((mat > 0).sum(axis=0)).ravel()
                return counts, genes, None
            raise ValueError(f"Unrecognized HDF5 structure in {h5_path}")


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    base = repo_root / 'output' / 'cellranger'
    if not base.exists():
        raise SystemExit(f"Expected directory not found: {base}")

    results = []
    for sample_dir in sorted(base.iterdir()):
        outs = sample_dir / 'outs'
        h5 = outs / 'filtered_feature_bc_matrix.h5'
        if h5.exists():
            print(f"Reading counts for sample: {sample_dir.name}")
            try:
                counts, genes_detected, barcodes = read_counts_from_h5(h5)
            except Exception as e:
                print(f"Failed to read {h5}: {e}")
                continue
            df = pd.DataFrame({'sample': sample_dir.name, 'counts': counts, 'genes_detected': genes_detected})
            results.append(df)
        else:
            print(f"No filtered_feature_bc_matrix.h5 for {sample_dir.name}, skipping")

    if not results:
        raise SystemExit("No sample matrices found to process.")

    all_df = pd.concat(results, ignore_index=True)

    out_dir = repo_root / 'output' / 'plots'
    out_dir.mkdir(parents=True, exist_ok=True)

    # Summary statistics for both counts and genes detected
    summary = all_df.groupby('sample').agg(
        total_cells=('counts', 'count'),
        counts_mean=('counts', 'mean'),
        counts_median=('counts', 'median'),
        counts_std=('counts', 'std'),
        counts_min=('counts', 'min'),
        counts_max=('counts', 'max'),
        genes_mean=('genes_detected', 'mean'),
        genes_median=('genes_detected', 'median'),
        genes_std=('genes_detected', 'std'),
        genes_min=('genes_detected', 'min'),
        genes_max=('genes_detected', 'max'),
    )
    summary_path = out_dir / 'read_and_genes_detected_distribution_summary.csv'
    summary.to_csv(summary_path)
    print(f"Wrote summary to {summary_path}")

    # Create a two-panel figure: left = total counts (log), right = genes detected (linear)
    sns.set(style='whitegrid')
    order = list(all_df['sample'].unique())
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax0 = axes[0]
    sns.violinplot(x='sample', y='counts', data=all_df, order=order, inner=None, color='0.9', ax=ax0)
    sns.boxplot(x='sample', y='counts', data=all_df, order=order, width=0.2, ax=ax0)
    ax0.set_yscale('log')
    ax0.set_ylabel('Total counts per cell (log scale)')
    ax0.set_xlabel('Sample')
    ax0.set_title('Per-cell total counts distribution (UMI counts)')
    # annotate number of cells per sample on the violin
    ylim0 = ax0.get_ylim()
    y_text0 = ylim0[1] * 0.95
    for i, sample in enumerate(order):
        n = int((all_df['sample'] == sample).sum())
        ax0.text(i, y_text0, f"n={n}", ha='center', va='top', fontsize=9, color='black')

    ax1 = axes[1]
    sns.violinplot(x='sample', y='genes_detected', data=all_df, order=order, inner=None, color='0.9', ax=ax1)
    sns.boxplot(x='sample', y='genes_detected', data=all_df, order=order, width=0.2, ax=ax1)
    ax1.set_ylabel('Number of genes detected per cell')
    ax1.set_xlabel('Sample')
    ax1.set_title('Per-cell genes-detected distribution across samples')
    # annotate number of cells per sample on the genes plot as well
    ylim1 = ax1.get_ylim()
    y_text1 = ylim1[1] * 0.95
    for i, sample in enumerate(order):
        n = int((all_df['sample'] == sample).sum())
        ax1.text(i, y_text1, f"n={n}", ha='center', va='top', fontsize=9, color='black')

    plt.tight_layout()
    out_png = out_dir / 'read_and_genes_detected_distribution.png'
    plt.savefig(out_png, dpi=150)
    print(f"Saved combined plot to {out_png}")

    # Filtered view: keep cells with >= 10000 reads per cell
    threshold = 10000
    filtered_df = all_df[all_df['counts'] >= threshold].copy()
    if filtered_df.empty:
        print(f"No cells remain after filtering with counts >= {threshold}")
    else:
        filtered_summary = filtered_df.groupby('sample').agg(
            total_cells=('counts', 'count'),
            counts_mean=('counts', 'mean'),
            counts_median=('counts', 'median'),
            counts_std=('counts', 'std'),
            counts_min=('counts', 'min'),
            counts_max=('counts', 'max'),
            genes_mean=('genes_detected', 'mean'),
            genes_median=('genes_detected', 'median'),
            genes_std=('genes_detected', 'std'),
            genes_min=('genes_detected', 'min'),
            genes_max=('genes_detected', 'max'),
        )
        filtered_summary_path = out_dir / f'read_and_genes_detected_distribution_summary_filtered_{int(threshold)}.csv'
        filtered_summary.to_csv(filtered_summary_path)
        print(f"Wrote filtered summary to {filtered_summary_path}")

        # Plot filtered
        sns.set(style='whitegrid')
        order_f = list(filtered_df['sample'].unique())
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        ax0 = axes[0]
        sns.violinplot(x='sample', y='counts', data=filtered_df, order=order_f, inner=None, color='0.9', ax=ax0)
        sns.boxplot(x='sample', y='counts', data=filtered_df, order=order_f, width=0.2, ax=ax0)
        ax0.set_yscale('log')
        ax0.set_ylabel('Total counts per cell (log scale)')
        ax0.set_xlabel('Sample')
        ax0.set_title(f'Per-cell total counts (>= {threshold} reads)')
        ylim0 = ax0.get_ylim()
        y_text0 = ylim0[1] * 0.95
        for i, sample in enumerate(order_f):
            n = int((filtered_df['sample'] == sample).sum())
            ax0.text(i, y_text0, f"n={n}", ha='center', va='top', fontsize=9, color='black')

        ax1 = axes[1]
        sns.violinplot(x='sample', y='genes_detected', data=filtered_df, order=order_f, inner=None, color='0.9', ax=ax1)
        sns.boxplot(x='sample', y='genes_detected', data=filtered_df, order=order_f, width=0.2, ax=ax1)
        ax1.set_ylabel('Number of genes detected per cell')
        ax1.set_xlabel('Sample')
        ax1.set_title(f'Per-cell genes-detected (cells >= {threshold} reads)')
        ylim1 = ax1.get_ylim()
        y_text1 = ylim1[1] * 0.95
        for i, sample in enumerate(order_f):
            n = int((filtered_df['sample'] == sample).sum())
            ax1.text(i, y_text1, f"n={n}", ha='center', va='top', fontsize=9, color='black')

        plt.tight_layout()
        out_png_f = out_dir / f'read_and_genes_detected_distribution_filtered_{int(threshold)}.png'
        plt.savefig(out_png_f, dpi=150)
        print(f"Saved filtered combined plot to {out_png_f}")


if __name__ == '__main__':
    main()
