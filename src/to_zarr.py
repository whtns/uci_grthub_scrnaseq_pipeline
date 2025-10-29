#!/usr/bin/env python


# import ome_types
# from tifffile import tiffcomment
import spatialdata_io
import spatialdata as sd
import matplotlib.pyplot as plt
import spatialdata_plot
import logging
from pathlib import Path
import shutil
import argparse
import sys



# ome_xml_string = tiffcomment(filename)
# ome_object = ome_types.from_xml(ome_xml_string)

# %%
spatialdata_io.__version__

# %% [markdown]
# ## Loading the data

# %% [markdown]
# A reader for Xenium data is available in `spatialdata-io`. We used it to parse and convert to Zarr a [Xenium dataset of Human Lung Cancer](https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard).
# 
# You can download the data from the link above, or from this [convenience python script](https://github.com/giovp/spatialdata-sandbox/blob/main/xenium_2.0.0_io/download.py) and convert it to spatialdata format with [this script](https://github.com/giovp/spatialdata-sandbox/blob/main/xenium_2.0.0_io/to_zarr.py), rename the `.zarr` store to `xenium.zarr` and place it in the current folder (in alternatively you can use symlinks to make the data visible).
# 
# The dataset used in this notebook is the latest Xenium Multimodal Cell Segmentation extension of the Xenium technology, available for data processed using Xenium Analyzer version 2.0.0 and when Cell Segmentation Kit was used. Nevertheless, the `xenium()` reader supports all the Xenium versions.

# %%
def write_zarr_if_missing(crop0, out_path, overwrite: bool = False, logger: logging.Logger | None = None) -> bool:
    """Write a SpatialData-like object to a zarr store only if it does not already exist.

    Returns True on success or if the output already exists and overwrite is False.
    Returns False on failure.

    Parameters
    - crop0: object with a .write(path) method (e.g. SpatialData)
    - out_path: path-like to write (will be converted to Path)
    - overwrite: if True, remove existing out_path and write anew
    - logger: logger instance to use for messages
    """
    out_path = Path(out_path)
    if logger is None:
        logger = logging.getLogger(__name__)

    if out_path.exists():
        if overwrite:
            logger.info("Overwriting existing output: %s", out_path)
            try:
                shutil.rmtree(out_path)
            except Exception as exc:  # pragma: no cover - filesystem dependent
                logger.exception("Failed to remove existing output %s: %s", out_path, exc)
                return False
        else:
            logger.warning("Output already exists, skipping (use --overwrite to replace): %s", out_path)
            return True

    logger.info("Writing the data to %s", out_path)
    try:
        crop0.write(out_path)
        logger.info("Write complete")
        return True
    except Exception as exc:  # pragma: no cover - runtime dependent
        logger.exception("Failed to write zarr to %s: %s", out_path, exc)
        return False


# %%
def main(argv: list[str] | None = None) -> int:
    """Command-line entry point.

    Accepts a xenium directory path and an optional --overwrite flag.

    Returns exit code 0 on success, non-zero on failure.
    """
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(description="Convert Xenium dataset to Zarr using spatialdata-io")
    parser.add_argument("xenium_dir", help="Path to the Xenium dataset directory (input)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing zarr output if present")
    args = parser.parse_args(argv)

    xenium_dir = args.xenium_dir

    # read xenium data directly in using spatialdata_io
    sdata = spatialdata_io.xenium(xenium_dir)

    zarr_path = str(xenium_dir) + ".zarr"
    success = write_zarr_if_missing(sdata, zarr_path, overwrite=args.overwrite)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
