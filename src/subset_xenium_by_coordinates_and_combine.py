#!/usr/bin/env python
"""Subset Xenium dataset by coordinates and combine filtered tables.

This script was converted from a notebook. Running it will:
- read or convert a Xenium dataset into a .zarr store
- load polygon selections from coordinates.csv
- query the SpatialData object for each polygon and write per-sample filtered tables
- concatenate filtered tables and write a combined h5ad

Usage:
    python src/subset_xenium_by_coordinates_and_combine.py.py /path/to/xenium_dir [--overwrite]
"""

from __future__ import annotations

import os
import logging
from pathlib import Path
import shutil
import argparse
import numpy as np
from tifffile import tiffcomment
import spatialdata_io
import spatialdata as sd
import spatialdata_plot
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity
from shapely import MultiPolygon, Polygon
from geopandas import GeoDataFrame
import pandas as pd
import anndata as ad
import grthub_tools as gt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Subset a Xenium dataset by coordinates and combine filtered tables. "
            "Provide the path to the Xenium dataset directory (the folder that contains coordinates.csv)."
        )
    )
    parser.add_argument(
        "-x",
        "--xenium_dir",
        default="../data/output-XETG00221__0069979__Adult_cohort__20251003__181017",
        required=True,
        help=("Path to the Xenium dataset directory. If omitted, a repository-relative default will be used."),
    )
    parser.add_argument(
        "-m",
        "--metadata",
        help=("A csv with sample metadata."),
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="When writing zarr or outputs, overwrite existing files if present.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    xenium_dir = Path(args.xenium_dir)
    zarr_path = Path(str(xenium_dir) + ".zarr")

        # Read polygon IDs and clean names
    polygon_ids = pd.read_csv(os.path.join(xenium_dir, "coordinates.csv"), skiprows=2).iloc[:, 0].unique()

    clean_name = lambda s: (
        s.str.lower()
        .str.replace("[^a-z0-9_ ]", "", regex=True)
        .str.strip()
        .str.replace(" ", "_", regex=False)
    )

    polygon_ids = clean_name(pd.Series(polygon_ids)).tolist()

    # Read polygons selection
    polygons = spatialdata_io.xenium_explorer_selection(os.path.join(xenium_dir, "coordinates.csv"), pixel_size=0.2125, return_list=False)

    if args.metadata is None:
        args.metadata = xenium_dir / "metadata.csv"

    metadata = args.metadata if args.metadata else None

    if os.path.exists(metadata):
        sample_metadata = pd.read_csv(metadata)
        sample_metadata['condition'] = sample_metadata['sample_id'].str.replace(r'_[0-9]+', '', regex=True)

    # Build GeoDataFrame of selections
    gdf = GeoDataFrame({"geometry": polygons})
    gdf = pd.concat([gdf, sample_metadata], axis=1)
    gdf = ShapesModel.parse(gdf, transformations={"aligned": Identity()})
    # gdf["polygon_id"] = polygon_ids

    # Convert to zarr if missing
    if not zarr_path.exists():
        sdata = spatialdata_io.xenium(xenium_dir)
        sdata.write(zarr_path, overwrite=args.overwrite)

    sdata = sd.read_zarr(zarr_path)
    sdata.shapes["polygons"] = gdf

    # Optional plotting (requires display backend)
    try:
        sdata.shapes["polygons"].plot()
    except Exception:
        pass

    from spatialdata import polygon_query
    from collections import defaultdict

    filtered_tables_unmerged: dict[str, list[ad.AnnData]] = defaultdict(list)
    filtered_tables: dict[str, ad.AnnData] = {}

    # Query per sample and write filtered tables
    out_dir = xenium_dir.parents[1] / Path("output/scanpy") / xenium_dir.name
    out_dir.mkdir(parents=True, exist_ok=True)

    for sample_id in sample_metadata.sample_id.unique().tolist():
        for i, polygon in sdata["polygons"][sdata["polygons"].sample_id == sample_id].geometry.items():
            table = polygon_query(sdata, polygon=polygon, target_coordinate_system="global")["table"]

            table.obs = table.obs.merge(sample_metadata.iloc[[i]], how="cross")

            filtered_tables_unmerged[sample_id].append(table)

            polygon_id = sdata["polygons"].loc[i,'polygon_id']

        if filtered_tables_unmerged[sample_id]:
            filtered_tables[sample_id] = ad.concat(filtered_tables_unmerged[sample_id])
            filtered_tables[sample_id].write(out_dir / f"filtered_table_{sample_id}.h5ad", compression="gzip")

    # Discover filtered tables (in output dir)
    filtered_tables_paths: dict[str, Path] = {}
    for file_path in out_dir.rglob("filtered_table_*.h5ad"):
        key = file_path.stem.replace("filtered_table_", "")
        filtered_tables_paths[key] = file_path

    # Load and report
    loaded_tables: dict[str, ad.AnnData] = {}
    for sample, path in filtered_tables_paths.items():
        loaded_tables[sample] = ad.read_h5ad(path)
        print(f"Sample: {sample}, Number of cells in filtered table: {loaded_tables[sample].n_obs}")

    if loaded_tables:
        combined_adata = ad.concat([adata for adata in loaded_tables.values()])
        combined_adata.obs["batch"] = combined_adata.obs["polygon_id"].astype("category")
        combined_adata.write(out_dir / "combined_filtered_tables.h5ad", compression="gzip")
    else:
        logging.getLogger(__name__).warning("No filtered tables found to combine.")
    # combined_adata has already been written in the branch above when tables exist


if __name__ == "__main__":
    main()
