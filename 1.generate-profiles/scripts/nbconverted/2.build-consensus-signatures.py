#!/usr/bin/env python
# coding: utf-8

# # Generate Consensus Signatures
# 
# **Gregory Way, 2019**
# 
# We do not have well-level information for the cell health data.
# Therefore, we cannot map to cell painting replicates.
# 
# Instead, we generate consensus signatures for each treatment.
# We generate consensus signatures in two ways.
# 
# 1. Median consensus
# 2. MODZ (moderated z-score) transform used in the L1000 analysis paper ([Subramanian et al. 2017](https://doi.org/10.1016/j.cell.2017.10.049)).
# 
# We apply these transformations to both:
# 
# * Cell Painting Data
# * Cell Health Assay Readout Data

# In[1]:


import os
import numpy as np
import pandas as pd

from pycytominer.consensus import modz
from pycytominer import get_na_columns, aggregate


# ## Load Cell Painting Data
# 
# This will be our x matrix in machine learning appications.

# In[2]:


profile_dir = os.path.join("data", "profiles")

all_profile_files = []
for plate in os.listdir(profile_dir):
    plate_dir = os.path.join(profile_dir, plate)
    
    if plate == '.DS_Store':
        continue

    for profile_file in os.listdir(plate_dir):
        if "feature_select" in profile_file:
            all_profile_files.append(os.path.join(plate_dir, profile_file))


# In[3]:


all_profile_files


# In[4]:


# Concatentate all cell painting datasets
x_df = (
    pd.concat(
        [pd.read_csv(x) for x in all_profile_files],
        sort=True
    )
    .rename(
        {
            "Image_Metadata_Plate": "Metadata_Plate",
            "Image_Metadata_Well": "Metadata_Well"
        },
        axis="columns")
    .drop(["Metadata_broad_sample"], axis="columns")
)

# Realign metadata column names
x_metadata_cols = x_df.columns[x_df.columns.str.startswith("Metadata")]
x_metadata_df = x_df.loc[:, x_metadata_cols]

x_df = x_df.drop(x_metadata_cols, axis="columns")
x_df = pd.concat([x_metadata_df, x_df], axis="columns")

# Drop columns with na values
na_cols_to_drop = get_na_columns(x_df, cutoff=0)
print("Dropping {} columns because of missing data".format(len(na_cols_to_drop)))
x_df = x_df.drop(na_cols_to_drop, axis="columns")

# Also drop Costes features
costes_cols_to_drop = [x for x in x_df.columns if "costes" in x.lower()]
print("Dropping {} costes features".format(len(costes_cols_to_drop)))
x_df = x_df.drop(costes_cols_to_drop, axis="columns")

print(x_df.shape)
x_df.head(2)


# ## Load Cell Health Assay Data
# 
# This will be the y matrix in machine learning applications.

# In[5]:


file = os.path.join("data", "labels", "normalized_cell_health_labels.tsv")
y_df = pd.read_csv(file, sep='\t').drop(["plate_name", "well_col", "well_row"], axis="columns")

print(y_df.shape)
y_df.head(2)


# ## Determine how many Cell Painting profiles have Cell Health status labels

# In[6]:


x_groupby_cols = ["Metadata_gene_name", "Metadata_pert_name", "Metadata_cell_line"]

x_metacount_df = (
    x_df
    .loc[:, x_groupby_cols]
    .assign(n_measurements=1)
    .groupby(x_groupby_cols)
    .count()
    .reset_index()
    .assign(data_type="cell_painting")
    .merge(x_df.loc[:, x_groupby_cols + ["Metadata_Well", "Metadata_Plate"]],
           how="left",
           on=x_groupby_cols)
)

print(x_metacount_df.shape)
x_metacount_df.head(2)


# In[7]:


y_groupby_cols = ["guide", "cell_id"]

y_metacount_df = (
    y_df
    .loc[:, y_groupby_cols]
    .assign(n_measurements=1)
    .groupby(y_groupby_cols)
    .count()
    .reset_index()
    .assign(data_type="cell_health")
)

print(y_metacount_df.shape)
y_metacount_df.head(2)


# In[8]:


all_measurements_df = (
    x_metacount_df
    .merge(
        y_metacount_df,
        left_on=["Metadata_pert_name", "Metadata_cell_line"],
        right_on=["guide", "cell_id"],
        suffixes=["_paint", "_health"],
        how="inner")
    .sort_values(by=["Metadata_cell_line", "Metadata_pert_name"])
    .reset_index(drop=True)
    .drop(["Metadata_Well", "guide", "cell_id"], axis="columns")
)

file = os.path.join("results", "all_profile_metadata.tsv")
all_measurements_df.to_csv(file, sep='\t', index=False)

print(all_measurements_df.shape)
all_measurements_df.head()


# ## A. Apply Median Consensus Aggregation
# 
# ### 1) To the Cell Painting Data

# In[9]:


x_median_df = aggregate(
    x_df,
    strata=["Metadata_cell_line", "Metadata_pert_name"],
    features="infer",
    operation="median"
)


x_median_df = (
    x_median_df
    .query("Metadata_pert_name in @all_measurements_df.Metadata_pert_name.unique()")
    .query("Metadata_cell_line in @all_measurements_df.Metadata_cell_line.unique()")
    .reset_index(drop=True)
    .reset_index()
    .rename({"index": "Metadata_profile_id"}, axis='columns')
)
x_median_df.Metadata_profile_id = ["profile_{}".format(x) for x in x_median_df.Metadata_profile_id]

print(x_median_df.shape)
x_median_df.head()


# In[10]:


# Output Profile Mapping for Downstream Analysis
profile_id_mapping_df = x_median_df.loc[:, x_median_df.columns.str.startswith("Metadata")]
file = os.path.join("data", "profile_id_metadata_mapping.tsv")
profile_id_mapping_df.to_csv(file, sep='\t', index=False)

profile_id_mapping_df.head()


# ### 2) To the Cell Health Assay Data

# In[11]:


cell_health_meta_features = ["cell_id", "guide"]
cell_health_features = y_df.drop(cell_health_meta_features, axis="columns").columns.tolist()
y_meta_merge_cols = ["Metadata_profile_id", "Metadata_pert_name", "Metadata_cell_line"]


# In[12]:


y_median_df = aggregate(
    y_df,
    strata=cell_health_meta_features,
    features=cell_health_features,
    operation="median"
)

print(y_median_df.shape)
y_median_df.head()


# In[13]:


y_median_df = (
    y_median_df
    .reset_index(drop=True)
    .merge(
        x_median_df.loc[:, y_meta_merge_cols],
        left_on=["guide", "cell_id"],
        right_on=["Metadata_pert_name", "Metadata_cell_line"],
        how="right"
    )
)

# Get columns in correct order
y_columns = (
    y_meta_merge_cols +
    y_median_df
    .loc[:, ~y_median_df.columns.str.startswith("Metadata_")]
    .columns
    .tolist()
)

y_median_df = (
    y_median_df
    .loc[:, y_columns]
    .drop(["guide", "cell_id"], axis="columns")
)

print(y_median_df.shape)
y_median_df.head(5)


# In[14]:


# Confirm that matrices are aligned
pd.testing.assert_series_equal(x_median_df.Metadata_profile_id,
                               y_median_df.Metadata_profile_id, check_names=True)

# Are the guides aligned?
pd.testing.assert_series_equal(x_median_df.Metadata_pert_name,
                               y_median_df.Metadata_pert_name, check_names=True)

# Are the cells aligned?
pd.testing.assert_series_equal(x_median_df.Metadata_cell_line,
                               y_median_df.Metadata_cell_line, check_names=True)


# ## B. Apply the MODZ Consensus Aggregation
# 
# ### 1) To the Cell Painting Data

# In[15]:


x_consensus_df = modz(
    x_df,
    replicate_columns=["Metadata_cell_line", "Metadata_pert_name"],
    precision=5
)

x_consensus_df.head()


# In[16]:


x_consensus_df = (
    x_consensus_df
    .reset_index()
    .query("Metadata_pert_name in @all_measurements_df.Metadata_pert_name.unique()")
    .query("Metadata_cell_line in @all_measurements_df.Metadata_cell_line.unique()")
    .reset_index(drop=True)
    .reset_index()
    .rename({"index": "Metadata_profile_id"}, axis='columns')
)
x_consensus_df.Metadata_profile_id = ["profile_{}".format(x) for x in x_consensus_df.Metadata_profile_id]

print(x_consensus_df.shape)
x_consensus_df.head(5)


# ### 2) To the Cell Health Assay Data

# In[17]:


y_consensus_df = modz(
    y_df,
    features=cell_health_features,
    replicate_columns=cell_health_meta_features,
    precision=5)

print(y_consensus_df.shape)
y_consensus_df.head()


# In[18]:


y_consensus_df = (
    y_consensus_df
    .reset_index()
    .reset_index(drop=True)
    .merge(
        x_consensus_df.loc[:, y_meta_merge_cols],
        left_on=["guide", "cell_id"],
        right_on=["Metadata_pert_name", "Metadata_cell_line"],
        how="right"
    )
    .loc[:, y_columns]
    .drop(["guide", "cell_id"], axis="columns")
)

print(y_consensus_df.shape)
y_consensus_df.head(5)


# In[19]:


# Confirm that matrices are aligned
pd.testing.assert_series_equal(x_consensus_df.Metadata_profile_id,
                               y_consensus_df.Metadata_profile_id, check_names=True)

# Are the guides aligned?
pd.testing.assert_series_equal(x_consensus_df.Metadata_pert_name,
                               y_consensus_df.Metadata_pert_name, check_names=True)

# Are the cells aligned?
pd.testing.assert_series_equal(x_consensus_df.Metadata_cell_line,
                               y_consensus_df.Metadata_cell_line, check_names=True)


# ## Output Median and MODZ Consensus Signatures

# In[20]:


consensus_dir = os.path.join("data", "consensus")

file = os.path.join(consensus_dir, "cell_painting_median.tsv.gz")
x_median_df.to_csv(file, sep="\t", index=False)

file = os.path.join(consensus_dir, "cell_health_median.tsv.gz")
y_median_df.to_csv(file, sep="\t", index=False)

file = os.path.join(consensus_dir, "cell_painting_modz.tsv.gz")
x_consensus_df.to_csv(file, sep="\t", index=False)

file = os.path.join(consensus_dir, "cell_health_modz.tsv.gz")
y_consensus_df.to_csv(file, sep="\t", index=False)

