#!/usr/bin/env python
# coding: utf-8

# ## Generate Morpheus Input Data
# 
# **Gregory Way, 2019**
# 
# Use this script to concatenate all of the cell painting data into one `.gct` file for input into morpheus.

# In[1]:


import os
import pandas as pd

from pycytominer import write_gct
from pycytominer.cyto_utils import infer_cp_features


# In[2]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[3]:


batch_id = "CRISPR_PILOT_B1"
backend_dir = os.path.join("..", "1.generate-profiles", "data", "profiles")

plate_dirs = [os.path.join(backend_dir, x) for x in os.listdir(backend_dir) if x != ".DS_Store"]
plate_dirs


# In[4]:


# Build full cell painting dataset
df_list = []
all_plate_files = {}
for plate_dir in plate_dirs:
    plate_files = os.listdir(plate_dir)
    for plate_file in plate_files:
        if "normalized_feature_select.csv" in plate_file:
            full_plate_file = os.path.join(plate_dir, plate_file)
            
            plate = plate_dir.split("/")[-1]
            all_plate_files[plate] = full_plate_file
            
            df = pd.read_csv(full_plate_file)
            print("reading {} with profile count: {}".format(plate_file, df.shape[0]))
            df_list.append(df)


# In[5]:


# Combine into a single file
cp_df = pd.concat(df_list, sort=True)#.reset_index(drop=True)
cp_features = infer_cp_features(cp_df)
meta_features = cp_df.drop(cp_features, axis="columns").columns.tolist()

cp_df = cp_df.loc[:, meta_features + cp_features]

print(cp_df.shape)
cp_df.head()


# # Extract Data for the Most Replicable Genes
# 
# The genes include: *ITGAV*, *KIF11*, *MYC*, *POLR2D*, and *PSMA1*. (Plus one control LacZ)

# In[6]:


# We are interested here in the most replicable CRISPR'd genes
genes = ['ITGAV', 'KIF11', 'MYC', 'POLR2D', 'PSMA1', 'LacZ']
cp_genes_df = cp_df.query("Metadata_gene_name in @genes").reset_index(drop=True)

cp_genes_df = (
    cp_genes_df
    .groupby(
        ['Metadata_cell_line', 'Metadata_gene_name', 'Metadata_pert_name']
    )
    .mean()
    .reset_index()
)

print(cp_genes_df.shape)
cp_genes_df.head(2)


# ## Create `.gct` files for Morpheus heatmap inputs

# In[7]:


# Build and output gct file for all genes
output_file = os.path.join(
    "results", "morpheus", "full_genes_morpheus.gct"
)

write_gct(
    profiles=cp_df,
    output_file=output_file,
    features=cp_features
)


# In[8]:


# Build and output gct file for select genes
output_file = os.path.join(
    "results", "morpheus", "reproducible_genes.gct"
)

write_gct(
    profiles=cp_genes_df,
    output_file=output_file,
    features=cp_features
)


# In[9]:


# Write a gct file for all plates
for plate in all_plate_files:

    df = pd.read_csv(all_plate_files[plate])
    
    output_file = os.path.join(
        "results", "morpheus", "{}_morpheus.gct".format(plate)
    )
    
    write_gct(
        profiles=df,
        output_file=output_file,
        features=cp_features
    )

