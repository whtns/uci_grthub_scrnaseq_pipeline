from pathlib import Path
from pyhere import here
import json
import os
import re
import scanpy as sc

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes, Url

#spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot
img_channels = "rgb"
default_channels = {"blue": "DAPI"}
default_gene = "SLC17A6"

#   Names of continuous features expected to be columns in the observation 
#   data (colData) of the AnnData
# spe_cont_features = ['PpTau', 'PAbeta']

#   Diagnosis by brain number (not included in the sample_info sheet)
# sample_dx = {
#     'Br3854': 'AD', 'Br3873': 'AD', 'Br3880': 'AD', 'Br3874': 'control'
# }

sample_info_path = here(
    'xenium_HYP','raw-data', 'sample-info', 'demos_by_sampleid-xen.xlsx'
)

spe_path = here('xenium_HYP',"processed-data", "12-samui", "spe.h5ad")
# notes_path = str(Path(here('code', '16_samui', 'feature_notes.md')).resolve())
img_path = here('xenium_HYP','processed-data', '11-XeniumSPE_imageunderlays', '{}')
#json_path = here(
#    'processed-data', '09_spaceranger_reorg', 'spg', '{}', 'outs', 'spatial',
#    'scalefactors_json.json'
#)

out_dir = here('xenium_HYP','processed-data', '12-samui', '{}')

################################################################################
#   Read in sample info and clean
################################################################################

#   Read in sample info, subset to relevant columns, and clean
# sample_info = (pd.read_excel(sample_info_path)
#     .query('`Sequenced? ` == "Yes"')
#     .filter(["Br####", "Slide SN #", "Array #", "Sample #"])
#     #   Clean up column names
#     .rename(
#         columns = {
#             "Br####": "br_num",
#             "Slide  #": "sample_id",
#             "Array #": "array_num",
#             "Sample #": "sample_num"
#         }
#     )
# )
sample_info = pd.read_excel(sample_info_path)
sample_info = sample_info.dropna(subset=['sample_id', 'Tissue'])

#   Prepend "Br" tp brain number and make a string
sample_info['br_num'] = sample_info['brnum'].astype(str)

#   Add diagnosis using brain number
# sample_info['diagnosis'] = sample_info['br_num'].replace(sample_dx)

#   Fix the experiment number column (use strings of integers), create experiment number for each slide
sample_info['sample_num'] = sample_info.index + 1
sample_info['experiment_num'] = ((sample_info['sample_num'] - 1) // 4  + 1).astype(str)

# create column called sample_id from 'Sample #'
#pattern = '_shk'
sample_info['sample_id'] = sample_info['sample_id']

#   Different forms of sample IDs appear to be used for spaceranger outputs
#   and raw images
# sample_info = (sample_info
#     .assign(
#         spaceranger_id = sample_info['sample_id'].transform(lambda x: x.replace('-', '')) +
#             '_' + sample_info['array_num'] + '_' + sample_info['br_num'],
#         image_id = 'VIFAD' + sample_info['experiment_num'] + '_' + sample_info['sample_id'] + '_' + sample_info['array_num']
#     )
# )
sample_info['spaceranger_id'] = sample_info['obj_id']
sample_info['image_id'] = sample_info['image'] 
# sample_id_spaceranger = sample_info['Slide #'] + '_' + sample_info['Array #']
# sample_id_image = sample_info['Slide #'] + '_' + sample_info['Array #']

#   Subset all types of IDs to this sample only
sample_id_spaceranger = sample_info['spaceranger_id'].iloc[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]
sample_id_image = sample_info['image_id'].iloc[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]
sample_id_samui = sample_info['spaceranger_id'].iloc[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]


#sample_id_spaceranger = sample_info['spaceranger_id'].iloc[int()]
# #   Update paths for this sample ID
out_dir = Path(str(out_dir).format(sample_id_samui))
#json_path = Path(str(json_path).format(sample_id_spaceranger))
img_path = Path(str(img_path).format(sample_id_image))

# #find size of dataframe and assign slide + array info to sample_id_spaceranger
# samplenums = sample_info['sample_num']
# samplerows = list(range(len(sample_info)))
# sample_id_spaceranger = sample_info['Slide #'] + '_' + sample_info['Array #']

# for i in samplerows:
#     sample_id_img = sample_id_spaceranger[i]
    
#     out_dir = Path(str(out_dir).format(sample_id_img))
#     json_path = Path(str(json_path).format(sample_id_img))
#     img_path = Path(str(img_path).format(sample_id_img))

out_dir.mkdir(exist_ok = True)

#   All paths should exist
assert all([x.exists() for x in [out_dir, img_path]])

################################################################################
#   Read in scale-factors info
################################################################################

#   Read in the spaceranger JSON to calculate meters per pixel for
#   the full-resolution image
#with open(json_path, 'r') as f:
#    spaceranger_json = json.load(f)

# m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']
# pixel / size(m)

################################################################################
#   Gather gene-expression data into a DataFrame to later as a feature
################################################################################

#   Read in AnnData and subset to this sample
spe = sc.read(spe_path)
#path_groups = spe.obs['path_groups'].cat.categories
spe = spe[spe.obs['sample_id'] == sample_id_spaceranger, :]
spe.obs.index.name = "barcode"

#   Convert the sparse gene-expression matrix to pandas DataFrame, with the
#   gene symbols as column names
gene_df = pd.DataFrame(
    spe.X.toarray(),
    index = spe.obs.index,
    columns = spe.var['Symbol']
)

#   Some gene symbols are actually duplicated. Just take the first column in
#   any duplicated cases
gene_df = gene_df.loc[: , ~gene_df.columns.duplicated()].copy()

#   Samui seems to break when using > ~ 5,000 genes. Take just the genes where
#   at least 10% of spots have nonzero counts
# gene_df = gene_df.loc[:, np.sum(gene_df > 0, axis = 0) > (gene_df.shape[0] * 0.1)].copy()

assert default_gene in gene_df.columns, "Default gene not in AnnData"

print('Using {} genes as features.'.format(gene_df.shape[1]))

################################################################################
#   Split 'path_groups' column into binary columns for each of its values
################################################################################

#   Circumvent a Samui bug (https://github.com/chaichontat/samui/issues/84);
#   turn the categorical column 'path_groups' into several numeric columns with
#   just values of 0 and 1
# path_df = pd.DataFrame()
# for path_group in path_groups:
#     path_df[path_group] = (spe.obs['path_groups'] == path_group).astype(int)

################################################################################
#   Use the Samui API to create the importable directory for this sample
################################################################################

this_sample = Sample(name = sample_id_samui, path = out_dir)

this_sample.add_coords(
    spe.obsm["spatial"].rename(
        columns = {"x_centroid": "x", "y_centroid": "y"}
    ),
    name = "coords",  size = 2e-5 , mPerPx = 1e-6
)

#   Add the IF image for this sample
#this_sample.add_image(
#    tiff = img_path, channels = img_channels, scale = 1e-6,
#    defaultChannels = default_channels
#)

#   Add gene expression results (multiple columns) as a feature
this_sample.add_chunked_feature(
    gene_df, name = "Genes", coordName = "coords", dataType = "quantitative"
)

#   Add additional requested observational columns (colData columns)
this_sample.add_csv_feature(
    spe.obs["Banksy_cluster"], name = "Banksy-Cluster", coordName = "coords",
    dataType = "categorical"
)

this_sample.add_csv_feature(
    spe.obs["banksy_cluster_group"], name = "Banksy Cluster Grouped", coordName = "coords",
    dataType = "categorical"
)

this_sample.add_csv_feature(
    spe.obs["domain_of_cell"], name = "Domain of Cell", coordName = "coords",
    dataType = "categorical"
)

this_sample.add_csv_feature(
    spe.obs["nCounts"], name = "nCounts", coordName = "coords",
    dataType = "quantitative"
)

this_sample.add_csv_feature(
    spe.obs["nGenes"], name = "nGenes", coordName = "coords",
    dataType = "quantitative"
)

this_sample.add_csv_feature(
    spe.obs["cell_area"], name = "Cell Area", coordName = "coords",
    dataType = "quantitative"
)

this_sample.add_csv_feature(
    spe.obs["nucleus_area"], name = "Nucleus Area", coordName = "coords",
    dataType = "quantitative"
)

#   Add pathology groups
# this_sample.add_csv_feature(
#     path_df, name = "Pathology Group", coordName = "coords",
#     dataType = "quantitative"
# )

this_sample.set_default_feature(group = "Genes", feature = default_gene)

this_sample.write()


#  /dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/xenium_HYP/code/13-samui/../../../raw-data/