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
import scipy.stats
from joblib import load
import umap

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features

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


# ## 2) Load Cell Painting Repurposing Data Files
# 
# These files are available from https://github.com/broadinstitute/lincs-cell-painting

# In[6]:


batch = "2016_04_01_a549_48hr_batch1"
commit_hash = "27a2d7dd74067b5754c2c045e9b1a9cfb0581ae4"

# We have noticed particular technical issues with this platemap
# remove it from downstream consideration
# https://github.com/broadinstitute/lincs-cell-painting/issues/43
filter_platemap = "C-7161-01-LM6-011"


# In[7]:


# Load data
base_url = "https://media.githubusercontent.com/media/broadinstitute/lincs-cell-painting/"
repurp_url = f"{base_url}/{commit_hash}/consensus/{batch}/{batch}_consensus_{consensus}.csv.gz"

complete_consensus_df = pd.read_csv(repurp_url)

print(complete_consensus_df.shape)
complete_consensus_df.head()


# In[8]:


# Apply feature selection to the consensus profiles
feature_ops = [
    "variance_threshold",
    "drop_na_columns",
    "blacklist",
    "drop_outliers"
]

consensus_feature_select_df = feature_select(
    complete_consensus_df,
    operation=feature_ops,
    na_cutoff=0
)

print(consensus_feature_select_df.shape)
consensus_feature_select_df.head()


# In[9]:


# Split metadata and CP Features
cp_features = infer_cp_features(x_test_df)
meta_features = infer_cp_features(complete_consensus_df, metadata=True)

# Realign LINCS data to the same feature ordering as the test dataset
feature_df = complete_consensus_df.reindex(cp_features, axis="columns")
metadata_df = complete_consensus_df.loc[:, meta_features]

print(feature_df.shape)
feature_df.head()


# ## 3) Apply all Regression Models to all Repurposing Plates

# In[10]:


cell_health_features = list(model_dict.keys())

all_scores = {}
for cell_health_feature in cell_health_features:
    model_clf = model_dict[cell_health_feature]
    pred_df = model_clf.predict(feature_df)
    all_scores[cell_health_feature] = pred_df


# ## 4) Output Results

# In[11]:


# Output scores
all_score_df = pd.DataFrame.from_dict(all_scores)
repurp_predict_df = (
    metadata_df
    .merge(
        all_score_df,
        left_index=True,
        right_index=True
    )
    .query("Metadata_Plate_Map_Name != @filter_platemap")
)

output_real_file = os.path.join(
    output_dir,
    "repurposing_transformed_real_models_{}.tsv.gz".format(consensus)
)
repurp_predict_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(repurp_predict_df.shape)
repurp_predict_df.head()


# ## 5) Apply UMAP
# 
# ### Part 1: Apply UMAP to Cell Health Transformed Repurposing Hub Features

# In[12]:


reducer = umap.UMAP(random_state=1234, n_components=2)

predict_embedding_df = pd.DataFrame(
    reducer.fit_transform(repurp_predict_df.loc[:, cell_health_features]),
    columns=["umap_x", "umap_y"]
)

predict_embedding_df = (
    metadata_df
    .merge(
        predict_embedding_df,
        left_index=True,
        right_index=True
    )
    .query("Metadata_Plate_Map_Name != @filter_platemap")
)

output_real_file = os.path.join(
    output_dir,
    "repurposing_umap_transformed_real_models_{}.tsv.gz".format(consensus)
)

predict_embedding_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(predict_embedding_df.shape)
predict_embedding_df.head()


# ### Part 2: Apply UMAP to All Repurposing Hub Cell Painting Profiles

# In[13]:


reducer = umap.UMAP(random_state=1234, n_components=2)

repurp_embedding_df = pd.DataFrame(
    reducer.fit_transform(
        consensus_feature_select_df.loc[:, infer_cp_features(consensus_feature_select_df)]
    ),
    columns=["umap_x", "umap_y"]
)

repurp_embedding_df = (
    metadata_df
    .merge(
        repurp_embedding_df,
        left_index=True,
        right_index=True
    )
    .query("Metadata_Plate_Map_Name != @filter_platemap")
)

output_real_file = os.path.join(
    output_dir,
    "repurposing_umap_transformed_cell_painting_{}.tsv.gz".format(consensus)
)
repurp_embedding_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(repurp_embedding_df.shape)
repurp_embedding_df.head()


# ## Merge Data Together for Shiny App Exploration

# In[14]:


# Load MOA file
moa_url = "https://raw.githubusercontent.com/broadinstitute/lincs-cell-painting/"
moa_url = f"{moa_url}/{commit_hash}/metadata/moa/repurposing_info_external_moa_map_resolved.tsv"

moa_df = pd.read_csv(moa_url, sep="\t")

print(moa_df.shape)
moa_df.head(3)


# In[15]:


core_id = [
    "{}-{}".format(
        x.split("-")[0],
        x.split("-")[1]
    ) if x != "DMSO"
    else x
    for x in repurp_embedding_df.Metadata_broad_sample
]

repurp_embedding_with_pert_df = (
    repurp_embedding_df
    .assign(Metadata_broad_core_id=core_id)
    .sort_index(axis="columns")
    .merge(
        moa_df,
        left_on="Metadata_broad_core_id",
        right_on="broad_id",
        how="left"
    )
)

print(repurp_embedding_with_pert_df.shape)
repurp_embedding_with_pert_df.head()


# In[16]:


shiny_merge_cols = [
    "Metadata_Plate_Map_Name",
    "Metadata_broad_sample",
    "Metadata_dose_recode",
    "Metadata_mmoles_per_liter",
    "Metadata_pert_well"
]

shiny_df = (
    repurp_embedding_with_pert_df.merge(
        repurp_predict_df,
        left_on=shiny_merge_cols,
        right_on=shiny_merge_cols,
        how="inner"
    )
    .drop(["broad_sample"], axis="columns")
    .query("Metadata_Plate_Map_Name != @filter_platemap")
)

print(shiny_df.shape)
shiny_df.head()


# In[17]:


shiny_file = os.path.join(
    "repurposing_cellhealth_shiny",
    "data",
    "moa_cell_health_{}.tsv.gz".format(consensus)
)

shiny_df.to_csv(shiny_file, sep='\t', index=False, compression="gzip")


# In[18]:


shiny_combined_df = shiny_df.merge(
    complete_consensus_df,
    on=infer_cp_features(complete_consensus_df, metadata=True),
    how="inner"
)


# ## Output Correlation Matrix

# In[19]:


shiny_features = infer_cp_features(consensus_feature_select_df)
cell_health_features = [x for x in shiny_df if x.startswith("cell_health")]


# In[20]:


all_results = []
for cell_health_feature in cell_health_features:
    cell_health = shiny_combined_df.loc[:, cell_health_feature]
    for cp_feature in shiny_features:
        feature = shiny_combined_df.loc[:, cp_feature]
        cor_result, pval = scipy.stats.pearsonr(cell_health, feature)
        all_results.append([cell_health_feature, cp_feature, cor_result, pval])


# In[21]:


# Output correlation matrix for cell health predictions and CellProfiler features
cor_results_df = (
    pd.DataFrame(
        np.array(all_results), columns=["cell_health", "cp_feature", "pearson_cor", "pval"]
    )
    .sort_values(by="pearson_cor", ascending=False)
    .reset_index(drop=True)
)

cor_results_df.pearson_cor = cor_results_df.pearson_cor.astype(float)

cor_results_df = (
    cor_results_df
    .pivot_table(columns=["cell_health"], index=["cp_feature"], values="pearson_cor")
)

print(cor_results_df.shape)
cor_results_df.head(3)

