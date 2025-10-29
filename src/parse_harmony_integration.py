#!/usr/bin/env python


import scanpy as sc
import scanpy.external as sce
from pathlib import Path
import anndata as ad
import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import grthub_tools as gt
import scvi

gt.integrate_w_harmony(combined_adata, "all_combined_harmony")
