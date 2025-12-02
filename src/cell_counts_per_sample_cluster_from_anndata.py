def cell_counts_per_sample_cluster_from_anndata(adata, sample_key, cluster_key):
    """
    Calculate the number of cells per sample and cluster from an AnnData object.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    sample_key (str): The key in adata.obs that contains sample identifiers.
    cluster_key (str): The key in adata.obs that contains cluster identifiers.

    Returns:
    pd.DataFrame: A DataFrame with samples as rows, clusters as columns,
                  and cell counts as values.
    """
    import pandas as pd

    # Create a DataFrame from the relevant obs columns
    df = adata.obs[[sample_key, cluster_key]].copy()

    # Create a contingency table (cross-tabulation)
    cell_counts = pd.crosstab(df[sample_key], df[cluster_key])

    return cell_counts
    