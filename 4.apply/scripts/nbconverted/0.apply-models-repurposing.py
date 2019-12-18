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


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


np.random.seed(123)


# ## 1) Load Models and Training Data

# In[4]:


model_dir = os.path.join("..", "3.train", "models")

model_dict, model_coef = load_models(model_dir=model_dir)
shuffle_model_dict, shuffle_model_coef = load_models(model_dir=model_dir, shuffle=True)


# In[5]:


data_dir = os.path.join("..", "3.train", "data")
x_train_df, x_test_df, y_train_df, y_test_df = load_train_test(data_dir=data_dir, drop_metadata=True)


# ## 2) Extract Repurposing Data Files
# 
# **NOTE** - these files are not yet public!

# In[6]:


# List drug repurposing data
repurposing_project_id = "2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad"

repurposing_profile_dir = os.path.join(
    "/Users",
    "gway",
    "work",
    "projects",
    repurposing_project_id,
    "workspace",
    "software",
    repurposing_project_id,
    "subsampling",
    "data",
    "profiles"
)


# In[7]:


# Build a single data frame that holds all profiles
plate_info = {}
all_dfs = []
all_metadata_dfs = []
all_plates = os.listdir(repurposing_profile_dir)
for plate in all_plates:
    plate_dir = os.path.join(repurposing_profile_dir, plate, "n_all")
    norm_file = os.path.join(plate_dir, "{}_subsample_all_normalized.csv".format(plate))
    plate_info[plate] = norm_file
    
    if os.path.exists(norm_file):
        df = pd.read_csv(norm_file)

        feature_df = df.reindex(x_test_df.columns, axis="columns").fillna(0)
        metadata_df = df.loc[:, df.columns.str.contains("Metadata_")]
        
        all_dfs.append(feature_df)
        all_metadata_dfs.append(metadata_df)


# In[8]:


metadata_df.head()


# In[9]:


# Merge feature data and metadata
all_df = pd.concat(all_dfs)
all_metadata_df = pd.concat(all_metadata_dfs)

complete_df = pd.concat([all_metadata_df, all_df], axis="columns").reset_index(drop=True)

# Fill in NaN in Metadata_broad_sample as DMSO
complete_df.Metadata_broad_sample = complete_df.Metadata_broad_sample.fillna("DMSO")

print(complete_df.shape)
complete_df.head()


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


repurp_cp_cols = complete_consensus_df.columns[~complete_consensus_df.columns.str.contains("Metadata")].tolist()
meta_cols = complete_consensus_df.drop(repurp_cp_cols, axis="columns").columns.tolist()


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
output_file = os.path.join("data", "repurposing_modz_consensus.tsv.gz")
complete_consensus_df.to_csv(output_file, sep='\t', compression="gzip", index=False)


# In[19]:


# Extract cell profiler and metadata features
cp_features = x_test_df.columns[~x_test_df.columns.str.startswith("Metadata")].tolist()


# ## 3) Apply all real and shuffled Models to all Repurposing Plates

# In[20]:


feature_df = complete_consensus_df.reindex(x_test_df.columns, axis="columns")
metadata_df = complete_consensus_df.loc[:, meta_cols]

all_scores = {}
all_shuffle_scores = {}
for cell_health_feature in model_dict.keys():
    # Apply Real Model Classifiers
    model_clf = model_dict[cell_health_feature]
    pred_df = model_clf.predict(feature_df)
    all_scores[cell_health_feature] = pred_df

    # Apply Shuffled Model Classifiers
    shuffle_model_clf = shuffle_model_dict[cell_health_feature]
    shuffle_pred_df = shuffle_model_clf.predict(feature_df)
    all_shuffle_scores[cell_health_feature] = shuffle_pred_df


# ## 4) Output Results

# In[21]:


output_dir = os.path.join("data", "repurposing_transformed")


# In[22]:


# Output scores
all_score_df = pd.DataFrame.from_dict(all_scores)
full_df = (
    metadata_df
    .merge(all_score_df,
           left_index=True,
           right_index=True)
)

output_real_file = os.path.join(output_dir, "repurposing_transformed_real_models.tsv.gz")
full_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(full_df.shape)
full_df.head()


# In[23]:


shuff_score_df = pd.DataFrame.from_dict(all_shuffle_scores)
full_shuff_df = (
    metadata_df
    .merge(shuff_score_df,
           left_index=True,
           right_index=True)
)

output_real_file = os.path.join(output_dir, "repurposing_transformed_shuffled_models.tsv.gz")
full_shuff_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(full_shuff_df.shape)
full_shuff_df.head()


# ## 5) Apply UMAP
# 
# ### Part 1: Apply UMAP to Cell Health Transformed Repurposing Hub Features

# In[24]:


cell_health_features = list(model_dict.keys())


# In[25]:


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

output_real_file = os.path.join(output_dir, "repurposing_umap_transformed_real_models.tsv.gz")
real_embedding_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")


# ### Part 2: Apply UMAP to All Repurposing Hub Cell Painting Profiles

# In[26]:


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

output_real_file = os.path.join(output_dir, "repurposing_umap_transformed_cell_painting.tsv.gz")
complete_embedding_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

