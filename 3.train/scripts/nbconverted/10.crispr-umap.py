#!/usr/bin/env python
# coding: utf-8

# ## Apply UMAP to CRISPR perturbations
# 
# These are the Cell Painting profiles of the CRISPR perturbations used to train and test each Cell Health model

# In[1]:


import umap
import pathlib
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


np.random.seed(123)


# In[3]:


# Set constants and file names
consensus = "modz"

data_dir = pathlib.Path("..", "1.generate-profiles", "data", "consensus")
cell_process_dir = pathlib.Path("..", "1.generate-profiles", "tables")
results_dir = pathlib.Path("results")
shiny_app_dir = pathlib.Path("..", "4.apply", "repurposing_cellhealth_shiny", "data")

profile_file = pathlib.Path(data_dir, f"cell_painting_{consensus}.tsv.gz")
cell_health_file = pathlib.Path(data_dir, f"cell_health_{consensus}.tsv.gz")

cell_process_file = pathlib.Path(cell_process_dir, "supplementary_table_1_perturbation_details.tsv")
output_file = pathlib.Path(shiny_app_dir, f"profile_umap_with_cell_health_{consensus}.tsv")


# In[4]:


# Load profile data
df = (
    pd.read_csv(profile_file, sep="\t")
    .sort_values(by="Metadata_profile_id")
    .reset_index(drop=True)
)

cp_features = infer_cp_features(df)

print(df.shape)
df.head()


# In[5]:


# Load cell health data
cell_health_df = (
    pd.read_csv(cell_health_file, sep="\t")
    .sort_values(by="Metadata_profile_id")
    .reset_index(drop=True)
)

print(cell_health_df.shape)
cell_health_df.head()


# In[6]:


# Load cell process annotation file
cell_process_df = pd.read_csv(cell_process_file, sep="\t")
cell_process_df.columns = [f"Metadata_{x}" for x in cell_process_df.columns]

print(cell_process_df.shape)
cell_process_df.head()


# In[7]:


# Ensure data and predictions are aligned
assert df.Metadata_profile_id.tolist() == cell_health_df.Metadata_profile_id.tolist()


# In[8]:


# Apply UMAP
reducer = umap.UMAP(random_state=1234, n_components=2)

predict_embedding_df = pd.DataFrame(
    reducer.fit_transform(df.loc[:, cp_features]),
    columns=["umap_x", "umap_y"]
)

print(predict_embedding_df.shape)
predict_embedding_df.head()


# In[9]:


# Combine data to form a single output file
output_df = (
    predict_embedding_df
    .merge(
        cell_health_df,
        left_index=True,
        right_index=True
    )
    .merge(
        cell_process_df,
        left_on="Metadata_pert_name",
        right_on="Metadata_pert_name",
        how="left"
    )
    # Drops 3 redundant "Empty" pert IDs
    .drop_duplicates(subset="Metadata_profile_id")
)

print(output_df.shape)
output_df.head()


# In[10]:


# Output to file
output_df.to_csv(output_file, sep="\t", index=False)

