"""AnnData helper utilities.

Provides a small helper to map gene symbols (in `adata.var` index) to
Ensembl gene IDs and append them as a new column `gene_ids`.

This implementation uses the Ensembl REST API and is intentionally
lightweight (per-symbol lookup). For large gene lists prefer a batched
approach or use `mygene` / local annotation tables.
"""
from typing import Optional
import time
import re

import pandas as pd


def add_ensembl_ids(adata, species: str = "human", var_name: str = "gene_ids", overwrite: bool = False, sleep: float = 0.05, timeout: float = 5.0) -> None:
    """Map `adata.var` index (gene symbols) to Ensembl IDs and add column.

    Parameters
    ----------
    adata
        An ``anndata.AnnData`` object with gene symbols as ``adata.var`` index.
    species
        Species string understood by Ensembl REST (e.g. ``"human"`` or ``"mouse"``).
    var_name
        Name of the new column to write into ``adata.var`` (default: ``gene_ids``).
    overwrite
        If False and ``var_name`` already exists, do nothing. If True overwrite.
    sleep
        Seconds to sleep between REST calls to be polite to Ensembl.
    timeout
        Per-request timeout (seconds).

    Notes
    -----
    - Uses the Ensembl REST endpoint ``/xrefs/symbol/{species}/{symbol}``.
    - For symbols with multiple mappings the first Ensembl gene ID is used.
    - If the index already appears to contain Ensembl IDs (many values
      starting with ``ENS``) the function will copy the index into
      ``adata.var[var_name]`` instead of querying the REST API.
    """
    try:
        import requests
    except Exception as e:
        raise ImportError("requests is required for add_ensembl_ids: install requests") from e

    if var_name in adata.var.columns and not overwrite:
        return

    index = adata.var.index.astype(str)

    # Quick heuristic: if most index values already look like Ensembl IDs, copy them
    ens_like = sum(1 for s in index if isinstance(s, str) and re.match(r"^ENS[A-Z]*G", s))
    if ens_like > len(index) * 0.5:
        adata.var[var_name] = index.values
        return

    headers = {"Accept": "application/json"}
    mapped = []

    for sym in index:
        if sym is None or sym == "" or (isinstance(sym, float) and pd.isna(sym)):
            mapped.append(None)
            continue

        url = f"https://rest.ensembl.org/xrefs/symbol/{species}/{sym}"
        try:
            r = requests.get(url, headers=headers, timeout=timeout)
            if r.status_code == 200:
                data = r.json()
                # Prefer entries that look like gene IDs (ENS*G)
                gene_ids = [d.get("id") for d in data if isinstance(d.get("id"), str) and re.match(r"^ENS[A-Z]*G", d.get("id"))]
                mapped.append(gene_ids[0] if gene_ids else None)
            else:
                mapped.append(None)
        except Exception:
            mapped.append(None)

        time.sleep(sleep)

    adata.var[var_name] = pd.Series(mapped, index=adata.var.index)


def add_ensembl_ids_if_symbols(adata, species: str = "human", var_name: str = "gene_ids", **kwargs) -> None:
    """Convenience wrapper: only add if `adata.var` index looks like gene symbols.

    This checks that few index values start with 'ENS' before attempting mapping.
    """
    index = adata.var.index.astype(str)
    ens_like = sum(1 for s in index if isinstance(s, str) and s.startswith("ENS"))
    if ens_like > len(index) * 0.5:
        # already Ensembl-like, copy index if column missing
        if var_name not in adata.var.columns:
            adata.var[var_name] = index.values
        return

    add_ensembl_ids(adata, species=species, var_name=var_name, **kwargs)
