#!/usr/bin/env python
# coding: utf-8

# # Normalize Cell Health Labels
# 
# **Gregory Way, 2019**

# In[1]:


import os
import numpy as np
from scipy.stats import median_absolute_deviation
import pandas as pd

from pycytominer import write_gct


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


# Function to scale cell health target variables
def scale(x):
    x_median = np.nanmedian(x)
    x_std = np.nanstd(x)
    x_scale = (x - x_median) / x_std
    return x_scale


# In[4]:


file = os.path.join("data", "labels", "cell_health_labels.tsv")
label_df = pd.read_csv(file, sep='\t')

print(label_df.shape)
label_df.head(2)


# In[5]:


# Some infinite values are present, replace them with NA
label_df = label_df.replace([np.inf, -np.inf], np.nan)

# Apply normalization by plate
normalized_label_df = (
    label_df
    .drop(["guide", "well_col", "well_row"], axis="columns")
    .groupby(["cell_id", "plate_name"])
    .transform(scale)
)

normalized_label_df = pd.concat(
    [
        label_df.loc[:, ["cell_id", "guide", "plate_name", "well_col", "well_row"]],
        normalized_label_df
    ],
    axis="columns"
)

print(normalized_label_df.shape)
normalized_label_df.head(2)


# In[6]:


# Write to file
file = os.path.join("data", "labels", "normalized_cell_health_labels.tsv")
normalized_label_df.to_csv(file, index=False, sep='\t')


# ## Build Cell Health Target Variable GCT Files
# 
# For viewing heatmaps in Morpheus.

# In[7]:


# Recode metadata variables
normalized_label_df = (
    normalized_label_df
    .rename(
        {
            "cell_id": "Metadata_cell_id",
            "guide": "Metadata_guide",
            "plate_name": "Metadata_plate_name",
            "well_col": "Metadata_well_col",
            "well_row": "Metadata_well_row"
        },
        axis="columns"
    )
)

print(normalized_label_df.shape)
normalized_label_df.head(2)


# In[8]:


# Load feature map
file = os.path.join("data", "labels", "feature_mapping_annotated.csv")
feature_map = (
    pd.read_csv(file, index_col=0)
    .rename(
        {
            "id": "variable_name"
        },
        axis="columns"
    )
    .transpose()
    .reset_index()
    .transpose()
    .rename(
        {
            "index": "id"
        },
        axis="rows"
    )
)

feature_map.head()


# In[9]:


# Build and output gct file
cell_health_features = [x for x in normalized_label_df.columns if not x.startswith("Metadata_")]
output_file = os.path.join("data", "labels", "normalized_cell_health_labels.gct")

write_gct(
    profiles=normalized_label_df,
    output_file=output_file,
    features=cell_health_features,
    meta_features="infer",
    feature_metadata=feature_map
)

