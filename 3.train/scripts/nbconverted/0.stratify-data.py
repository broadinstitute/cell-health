#!/usr/bin/env python
# coding: utf-8

# # Stratify Data into Training and Testing Sets
# 
# **Gregory Way, 2019**
# 
# Split the input data into training and testing sets balanced by guide infection.

# In[1]:


import os
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from pycytominer.get_na_columns import get_na_columns


# In[2]:


np.random.seed(123)


# ## Load X Matrices

# In[3]:


batch = "CRISPR_PILOT_B1"
data_dir = os.path.join("..", "0.generate-profiles", "data")
profile_dir = os.path.join(data_dir, "profiles", batch)

all_profile_files = []
for plates in os.listdir(profile_dir):
    plate_dir = os.path.join(profile_dir, plates)
    for profile_file in os.listdir(plate_dir):
        if "feature_select" in profile_file:
            all_profile_files.append(os.path.join(plate_dir, profile_file))


# In[4]:


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
)

# Drop all features that have missing values
additional_exclude_features = get_na_columns(x_df, features="infer", cutoff=0)
x_df = x_df.drop(additional_exclude_features, axis="columns")

print(x_df.shape)
x_df.head(2)


# ## Load Y Matrix

# In[5]:


file = os.path.join(data_dir, "labels", "normalized_cell_health_labels.tsv")
y_df = pd.read_csv(file, sep='\t').drop(["plate_name", "well_col", "well_row"], axis="columns")

print(y_df.shape)
y_df.head(2)


# ## Determine how many profiles have status labels

# In[6]:


x_groupby_cols = ["Metadata_gene_name", "Metadata_pert_name", "Metadata_cell_line"]


# In[7]:


x_meta_df = (
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

print(x_meta_df.shape)
x_meta_df.head(8)


# In[8]:


y_groupby_cols = ["guide", "cell_id"]


# In[9]:


y_meta_df = (
    y_df
    .loc[:, y_groupby_cols]
    .assign(n_measurements=1)
    .groupby(y_groupby_cols)
    .count()
    .reset_index()
    .assign(data_type="cell_health")
)

print(y_meta_df.shape)
y_meta_df.head(8)


# In[10]:


all_measurements_df = (
    x_meta_df
    .merge(
        y_meta_df,
        left_on=["Metadata_pert_name", "Metadata_cell_line"],
        right_on=["guide", "cell_id"],
        suffixes=["_paint", "_health"],
        how="inner")
    .sort_values(by=["Metadata_cell_line", "Metadata_pert_name"])
    .reset_index(drop=True)
)

file = os.path.join("results", "all_profile_metadata.tsv")
all_measurements_df.to_csv(file, sep='\t', index=False)

print(all_measurements_df.shape)
all_measurements_df.head()


# ## Aggregate Profiles and Outcomes Further
# 
# Because the plates do not match (no way to map wells across experiments), we must aggregate the ~6 cell painting replicates per guide and ~4 cell health replicates per guide together to form a single profile and single outcome.

# In[11]:


x_columns = x_groupby_cols + x_df.loc[:, ~x_df.columns.str.startswith("Metadata_")].columns.tolist()


# In[12]:


x_agg_df = (
    x_df
    .loc[:, x_columns]
    .groupby(x_groupby_cols)
    .median()
    .reset_index()
    .query("Metadata_gene_name in @all_measurements_df.Metadata_gene_name.unique()")
    .query("Metadata_pert_name in @all_measurements_df.Metadata_pert_name.unique()")
    .query("Metadata_cell_line in @all_measurements_df.Metadata_cell_line.unique()")
    .sort_values(by=["Metadata_cell_line", "Metadata_pert_name"])
    .reset_index(drop=True)
    .reset_index()
    .rename({"index": "Metadata_profile_id"}, axis='columns')
)

x_agg_df.Metadata_profile_id = ["profile_{}".format(x) for x in x_agg_df.Metadata_profile_id]


print(x_agg_df.shape)
x_agg_df.head(5)


# In[13]:


y_meta_cols = ["Metadata_profile_id", "Metadata_gene_name", "Metadata_pert_name", "Metadata_cell_line"]

y_agg_df = (
    y_df
    .groupby(y_groupby_cols)
    .median()
    .reset_index()
    .query("guide in @all_measurements_df.Metadata_pert_name.unique()")
    .query("cell_id in @all_measurements_df.Metadata_cell_line.unique()")
    .sort_values(by=["cell_id", "guide"])
    .reset_index(drop=True)
    .merge(
        x_agg_df.loc[:, y_meta_cols],
        left_on=["guide", "cell_id"],
        right_on=["Metadata_pert_name", "Metadata_cell_line"]
    )
)

y_columns = y_meta_cols + y_agg_df.loc[:, ~y_agg_df.columns.str.startswith("Metadata_")].columns.tolist()
y_agg_df = y_agg_df.loc[:, y_columns].drop(["guide", "cell_id"], axis="columns")

print(y_agg_df.shape)
y_agg_df.head(2)


# In[14]:


# Confirm that matrices are aligned
pd.testing.assert_series_equal(x_agg_df.Metadata_profile_id, y_agg_df.Metadata_profile_id, check_names=False)

# Are the guides aligned?
pd.testing.assert_series_equal(x_agg_df.Metadata_pert_name, y_agg_df.Metadata_pert_name, check_names=False)

# Are the cells aligned?
pd.testing.assert_series_equal(x_agg_df.Metadata_cell_line, y_agg_df.Metadata_cell_line, check_names=False)


# ## Split into Training and Testing

# In[15]:


test_proportion = 0.15


# In[16]:


x_train_df, x_test_df, y_train_df, y_test_df = train_test_split(
    x_agg_df,
    y_agg_df,
    test_size=test_proportion,
    random_state=42)


# In[17]:


print(x_train_df.shape)
print(x_test_df.shape)


# In[18]:


file = os.path.join("data", "x_train.tsv.gz")
x_train_df.to_csv(file, sep="\t", index=False)

file = os.path.join("data", "y_train.tsv.gz")
y_train_df.to_csv(file, sep="\t", index=False)

file = os.path.join("data", "x_test.tsv.gz")
x_test_df.to_csv(file, sep="\t", index=False)

file = os.path.join("data", "y_test.tsv.gz")
y_test_df.to_csv(file, sep="\t", index=False)

