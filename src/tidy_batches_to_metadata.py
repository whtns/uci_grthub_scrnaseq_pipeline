#!/usr/bin/env python 
import pandas as pd
from typing import List
import scanpy as sc

def create_sequential_timepoint_table(batch_list: List[str]) -> pd.DataFrame:
    """
    Converts a list of batchs formatted as {patient}_{date} into a DataFrame
    with sequential timepoints based on the date order for each patient.

    Args:
        batch_list: A list of strings containing patient batchs and dates.

    Returns:
        A pandas DataFrame with 'batch', 'patient', and 'timepoint' columns.
    """
    # 1. Create the initial DataFrame from the list of batchs
    df = pd.DataFrame(batch_list, columns=['batch'])

    # 2. Split the batch string into patient and Date columns
    # We use 'expand=True' to create new columns directly
    df[['patient', 'Date_Str']] = df['batch'].str.split('_', expand=True)
    df['patient'] = df['patient'].astype("category")

    # 3. Convert the date string (MMDDYY) into a proper datetime object
    # This is crucial for correct chronological sorting.
    df['timepoint_Date'] = pd.to_datetime(df['Date_Str'], format='%m%d%y')

    # 4. Sort the DataFrame by patient batch and then by the Date
    df_sorted = df.sort_values(by=['patient', 'timepoint_Date']).reset_index(drop=True)

    # 5. Group by patient and assign a sequential timepoint (rank)
    # The 'rank' function assigns sequential numbers (1, 2, 3...) based on the order
    # within each group (method='first' ensures consistent tie-breaking, though dates should be unique here).
    df_sorted['timepoint'] = df_sorted.groupby('patient')['timepoint_Date'].rank(
        method='first', 
        ascending=True
    ).astype("int64").astype('category')

    # 6. Select and return the final desired columns
    return df_sorted[['batch', 'patient', 'timepoint']]

adata = sc.read("output/scanpy/combined_harmony_integrated.h5ad")

meta = create_sequential_timepoint_table(adata.obs['batch'].unique().tolist())

meta.to_csv("metadata/metadata.csv", index=True)

adata.obs = pd.merge(adata.obs, meta, on="batch")

adata.write("output/scanpy/copy_combined_harmony_integrated.h5ad")


