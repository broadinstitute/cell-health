#!/usr/bin/env python
# coding: utf-8

# # Apply Cell Health Models to Repurposing Set
# 
# **Gregory Way, 2019**
# 
# The models are trained to predict cell health phenotypes.
# Here, I apply the models to Cell Painting data from the repurposing set.
# 
# I will use these predictions to identify compound perturbation signatures of cell health impact.

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
from joblib import load
import umap

from pycytominer.consensus import modz

sys.path.append("../3.train")
from scripts.ml_utils import load_train_test, load_models


# In[2]:


np.random.seed(123)


# ## 1) Load Models and Training Data

# In[3]:


consensus = "modz"
output_dir = "data"


# In[4]:


model_dir = os.path.join("..", "3.train", "models")

model_dict, model_coef = load_models(
    model_dir=model_dir,
    consensus=consensus
)


# In[5]:


data_dir = os.path.join("..", "3.train", "data")

x_train_df, x_test_df, y_train_df, y_test_df = load_train_test(
    data_dir=data_dir,
    consensus=consensus,
    drop_metadata=True
)


# ## 2) Extract Repurposing Data Files
# 
# **NOTE** - these files are not yet public!

# In[6]:


# List drug repurposing data
repurposing_project_id = "2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad"

repurposing_profile_dir = os.path.join(
    "/home",
    "ubuntu",
    "efs",
    repurposing_project_id,
    "workspace",
    "software",
    repurposing_project_id,
    "subsampling",
    "full_profile_data"
)

all_plates = list(set([x.split("_")[0] for x in os.listdir(repurposing_profile_dir)]))


# In[7]:


# Build a single data frame that holds all profiles
plate_info = {}
all_dfs = []
all_metadata_dfs = []
for plate in all_plates:
    norm_file = os.path.join(repurposing_profile_dir, 
                             "{}_subsample_all_normalized.csv".format(plate))

    plate_info[plate] = norm_file
    
    if os.path.exists(norm_file):
        df = pd.read_csv(norm_file)

        feature_df = df.reindex(x_test_df.columns, axis="columns").fillna(0)
        metadata_df = df.loc[:, df.columns.str.contains("Metadata_")]
        
        all_dfs.append(feature_df)
        all_metadata_dfs.append(metadata_df)


# In[8]:


# Merge feature data and metadata
all_df = pd.concat(all_dfs, sort=True)
all_metadata_df = pd.concat(all_metadata_dfs, sort=True)

complete_df = pd.concat([all_metadata_df, all_df], axis="columns").reset_index(drop=True)

# Fill in NaN in Metadata_broad_sample as DMSO
complete_df.Metadata_broad_sample = complete_df.Metadata_broad_sample.fillna("DMSO")

print(complete_df.shape)
complete_df.head()


# In[9]:


# Confirm that all plates are loaded
assert (
    sorted(list(complete_df.Image_Metadata_Plate.unique())) == sorted(all_plates)
)


# ## Recode Dose Information

# In[10]:


def recode_dose(x, doses, return_level=False):
    closest_index = np.argmin([np.abs(dose - x) for dose in doses])
    if np.isnan(x):
        return 0
    if return_level:
        return closest_index + 1
    else:
        return doses[closest_index]


# In[11]:


primary_dose_mapping = [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]


# In[12]:


complete_df = complete_df.assign(
    Metadata_dose_recode=(
        complete_df
        .Metadata_mmoles_per_liter
        .apply(
            lambda x: recode_dose(x, primary_dose_mapping, return_level=True)
        )
    )
)

print(complete_df.shape)
complete_df.head()


# In[13]:


complete_df.Metadata_dose_recode.value_counts()


# ## Create Consensus Profiles
# 
# ### a) Generate different consensus profiles for DMSO
# 
# Include Well Level Information

# In[14]:


replicate_cols = ["Metadata_broad_sample", "Metadata_dose_recode", "Image_Metadata_Well"]

dmso_consensus_df = modz(
    complete_df.query("Metadata_broad_sample == 'DMSO'"),
    features="infer",
    replicate_columns=replicate_cols,
    precision=5
)

dmso_consensus_df = dmso_consensus_df.reset_index()

print(dmso_consensus_df.shape)
dmso_consensus_df.head(2)


# ### b) Generate consensus profiles for all treatments

# In[15]:


replicate_cols = ["Metadata_broad_sample", "Metadata_dose_recode"]

complete_consensus_df = modz(
    complete_df.query("Metadata_broad_sample != 'DMSO'"),
    features="infer",
    replicate_columns=replicate_cols,
    precision=5
)

complete_consensus_df = complete_consensus_df.reset_index()
complete_consensus_df = complete_consensus_df.assign(Image_Metadata_Well="collapsed")

print(complete_consensus_df.shape)
complete_consensus_df.head(2)


# ### c) Merge Together

# In[16]:


repurp_cp_cols = (
    complete_consensus_df
    .columns
    [~complete_consensus_df.columns.str.contains("Metadata")]
    .tolist()
)

meta_cols = (
    complete_consensus_df
    .drop(repurp_cp_cols, axis="columns")
    .columns
    .tolist()
)


# In[17]:


complete_consensus_df = (
    pd.concat(
        [
            complete_consensus_df,
            dmso_consensus_df
        ],
        sort=True
    )
    .reset_index(drop=True)
)

complete_consensus_df = complete_consensus_df.loc[:, meta_cols + repurp_cp_cols]

print(complete_consensus_df.shape)
complete_consensus_df.head()


# ### d) Output Profiles

# In[18]:


# Output consensus profiles
output_file = os.path.join(output_dir, "repurposing_{}_consensus.tsv.gz".format(consensus))
complete_consensus_df.to_csv(output_file, sep='\t', compression="gzip", index=False)


# In[19]:


# Extract cell profiler and metadata features
cp_features = x_test_df.columns[~x_test_df.columns.str.startswith("Metadata")].tolist()


# ## 3) Apply all Regression Models to all Repurposing Plates

# In[20]:


feature_df = complete_consensus_df.reindex(x_test_df.columns, axis="columns")
metadata_df = complete_consensus_df.loc[:, meta_cols]

all_scores = {}
for cell_health_feature in model_dict.keys():
    # Apply Real Model Classifiers
    model_clf = model_dict[cell_health_feature]
    pred_df = model_clf.predict(feature_df)
    all_scores[cell_health_feature] = pred_df


# ## 4) Output Results

# In[21]:


# Output scores
all_score_df = pd.DataFrame.from_dict(all_scores)
full_df = (
    metadata_df
    .merge(all_score_df,
           left_index=True,
           right_index=True)
)

output_real_file = os.path.join(output_dir,
                                "repurposing_transformed_real_models_{}.tsv.gz".format(consensus))
full_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(full_df.shape)
full_df.head()


# ## 5) Apply UMAP
# 
# ### Part 1: Apply UMAP to Cell Health Transformed Repurposing Hub Features

# In[22]:


cell_health_features = list(model_dict.keys())


# In[23]:


reducer = umap.UMAP(random_state=1234, n_components=2)

metadata_df = full_df.drop(cell_health_features, axis="columns")

real_embedding_df = pd.DataFrame(
    reducer.fit_transform(full_df.loc[:, cell_health_features]),
    columns=["umap_x", "umap_y"]
)

real_embedding_df = (
    metadata_df
    .merge(real_embedding_df,
           left_index=True,
           right_index=True)
)

output_real_file = os.path.join(output_dir,
                                "repurposing_umap_transformed_real_models_{}.tsv.gz".format(consensus))
real_embedding_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")


# ### Part 2: Apply UMAP to All Repurposing Hub Cell Painting Profiles

# In[24]:


reducer = umap.UMAP(random_state=1234, n_components=2)

complete_metadata_df = complete_consensus_df.drop(cp_features, axis="columns")

complete_embedding_df = pd.DataFrame(
    reducer.fit_transform(complete_consensus_df.loc[:, cp_features]),
    columns=["umap_x", "umap_y"]
)

complete_embedding_df = (
    complete_metadata_df
    .merge(complete_embedding_df,
           left_index=True,
           right_index=True)
)

output_real_file = os.path.join(output_dir,
                                "repurposing_umap_transformed_cell_painting_{}.tsv.gz".format(consensus))
complete_embedding_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")


# ## Merge Data Together for Shiny App Exploration

# In[25]:


shiny_merge_cols = ["Metadata_broad_sample", "Metadata_dose_recode", "Image_Metadata_Well"]

shiny_df = real_embedding_df.merge(
    full_df,
    left_on=shiny_merge_cols,
    right_on=shiny_merge_cols,
    how="inner"
)

print(shiny_df.shape)
shiny_df.head()


# In[26]:


shiny_file = os.path.join("repurposing_cellhealth_shiny",
                          "data",
                          "moa_cell_health_{}.tsv.gz".format(consensus))

shiny_df.to_csv(shiny_file, sep='\t', index=False, compression="gzip")

