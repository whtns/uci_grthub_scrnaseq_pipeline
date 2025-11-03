"""Small helpers for writing pandas DataFrames / dict-of-DataFrames to Excel files.

Provides `save_dfs_to_excel` which writes each DataFrame in a dict to its own
worksheet and optionally auto-sizes columns. Sheet names are sanitized to
meet Excel limits.

Usage:
    from scripts.pandas_helpers import save_dfs_to_excel
    save_dfs_to_excel(res, "out.xlsx")
"""
from typing import Dict
import re

import pandas as pd


def sanitize_sheet_name(name: str) -> str:
    """Sanitize a sheet name for Excel: remove forbidden characters and trim to 31 chars."""
    # Forbidden characters: : \ / ? * [ ] and leading/trailing single quotes are problematic
    s = re.sub(r"[:\\\\\/?*\[\]]", "_", str(name))
    s = s.strip()
    if len(s) == 0:
        s = "sheet"
    return s[:31]


def save_dfs_to_excel(
    dfs: Dict[str, pd.DataFrame],
    out_path: str,
    engine: str = "openpyxl",
    auto_width: bool = True,
    sanitize_names: bool = True,
) -> None:
    """Write a dict of DataFrames to an Excel file with one sheet per key.

    Parameters
    ----------
    dfs
        Mapping of sheet_name -> DataFrame.
    out_path
        Output .xlsx path.
    engine
        pandas Excel writer engine. 'openpyxl' or 'xlsxwriter' recommended.
    auto_width
        If True, set basic column widths based on the longest cell in each column.
    sanitize_names
        If True, sanitize sheet names to be valid Excel names.

    Notes
    -----
    - For auto_width with xlsxwriter, column widths use set_column.
    - For openpyxl, column widths are adjusted on the worksheet.column_dimensions.
    """
    if not isinstance(dfs, dict):
        raise TypeError("dfs must be a dict of sheet_name -> pandas.DataFrame")

    # Use pandas ExcelWriter which handles engine-specific details
    with pd.ExcelWriter(out_path, engine=engine) as writer:
        for sheet_name, df in dfs.items():
            if sanitize_names:
                sheet = sanitize_sheet_name(sheet_name)
            else:
                sheet = str(sheet_name)[:31]

            # Ensure DataFrame is written
            df.to_excel(writer, sheet_name=sheet, index=False)

            if not auto_width:
                continue

            # Adjust column widths
            try:
                if engine == "xlsxwriter":
                    workbook = writer.book
                    worksheet = writer.sheets[sheet]
                    for idx, col in enumerate(df.columns):
                        col_data = df[col].astype(str).fillna("")
                        max_len = max(col_data.map(len).max() if len(col_data) else 0, len(str(col)))
                        # set_column takes zero-indexed columns
                        worksheet.set_column(idx, idx, max_len + 2)
                else:
                    # openpyxl path
                    try:
                        from openpyxl.utils import get_column_letter

                        workbook = writer.book
                        worksheet = writer.sheets[sheet]
                        for idx, col in enumerate(df.columns, start=1):
                            col_letter = get_column_letter(idx)
                            col_data = df[col].astype(str).fillna("")
                            max_len = max(col_data.map(len).max() if len(col_data) else 0, len(str(col)))
                            worksheet.column_dimensions[col_letter].width = max_len + 2
                    except Exception:
                        # If something goes wrong adjusting openpyxl widths, ignore but keep file
                        pass
            except Exception:
                # Be tolerant: writing should not fail due to width adjustments
                pass
