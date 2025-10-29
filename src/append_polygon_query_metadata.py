#!/usr/bin/env python
import pandas as pd
import scanpy as sc
from pathlib import Path
import re
import anndata as ad

xenium_dir = Path("data/output-XETG00221__0069979__Adult_cohort__20251003__181017")

meta = pd.read_csv(xenium_dir / "metadata.csv")

meta['condition'] = meta['sample_id'].str.replace(r'_[0-9]+', '', regex=True)

scanpy_dir = Path("output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017")

adata_paths = list(scanpy_dir.glob("**/filtered_table_*.h5ad"))

def prep_adata(adata_path: Path, meta: pd.DataFrame) -> sc.AnnData:
    adata = sc.read(adata_path)
    adata.obs = adata.obs.drop(columns=["polygon_id"])
    sample_id = adata_path.stem.split("filtered_table_")[-1]
    adata.obs["sample_id"] = sample_id
    return adata

adatas = [prep_adata(adata_path, meta) for adata_path in adata_paths]

combined_adata = ad.concat(adatas, axis=0)

pd.merge(combined_adata.obs, meta, on = "sample_id", how= "inner")
