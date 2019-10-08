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

sys.path.append("../3.train")
from scripts.ml_utils import load_train_test, load_models


# In[2]:


np.random.seed(123)


# ## 1) Load Models and Training Data

# In[3]:


model_dir = os.path.join("..", "3.train", "models")

model_dict, model_coef = load_models(model_dir=model_dir)
shuffle_model_dict, shuffle_model_coef = load_models(model_dir=model_dir, shuffle=True)


# In[4]:


data_dir = os.path.join("..", "3.train", "data")
x_train_df, x_test_df, y_train_df, y_test_df = load_train_test(data_dir=data_dir, drop_metadata=True)


# ## 2) Extract Repurposing Data Files

# In[5]:


# List drug repurposing data
repurposing_profile_dir = os.path.join(
    "/Users",
    "gway",
    "work",
    "projects",
    "2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad",
    "workspace",
    "software",
    "2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad",
    "subsampling",
    "data",
    "profiles"
)


# In[6]:


plate_info = {}
all_plates = os.listdir(repurposing_profile_dir)
for plate in all_plates:
    plate_dir = os.path.join(repurposing_profile_dir, plate, "n_all")
    norm_file = os.path.join(plate_dir, "{}_subsample_all_normalized.csv".format(plate))
    plate_info[plate] = norm_file


# ## 3) Apply all real and shuffled Models to all Repurposing Plates

# In[7]:


output_dir = os.path.join("data", "repurposing_transformed")


# In[8]:


real_models = []
shuffled_models = []
all_dfs = []
all_metadata_dfs = []
for plate in all_plates:
    norm_file = plate_info[plate]
    if os.path.exists(norm_file):
        df = pd.read_csv(norm_file)

        feature_df = df.reindex(x_test_df.columns, axis="columns").fillna(0)
        metadata_df = df.loc[:, df.columns.str.startswith("Metadata_")]
        
        all_dfs.append(feature_df)
        all_metadata_dfs.append(metadata_df)
        
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
    
        # Output scores
        all_score_df = pd.DataFrame.from_dict(all_scores)
        full_df = (
            metadata_df
            .merge(all_score_df,
                   left_index=True,
                   right_index=True)
            .assign(Metadata_plate=plate)
        )
        real_models.append(full_df)
            
        shuff_score_df = pd.DataFrame.from_dict(all_shuffle_scores)
        full_shuff_df = (
            metadata_df
            .merge(shuff_score_df,
                   left_index=True,
                   right_index=True)
            .assign(Metadata_plate=plate)
        )
        shuffled_models.append(full_shuff_df)
        
    else:
        print(plate)


# ## 4) Output Results

# In[9]:


all_df = pd.concat(all_dfs)
all_metadata_df = pd.concat(all_metadata_dfs)

complete_df = pd.concat([all_metadata_df, all_df], axis="columns").reset_index(drop=True)

print(complete_df.shape)
complete_df.head()


# In[10]:


full_real_df = pd.concat(real_models)

# Determine proper alignment of columns
output_cols = full_real_df.columns.tolist()
output_cols.insert(0, output_cols.pop(output_cols.index("Metadata_plate")))
full_real_df = full_real_df.loc[:, output_cols].reset_index(drop=True)

output_real_file = os.path.join(output_dir, "repurposing_transformed_real_models.tsv.gz")
full_real_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(full_real_df.shape)
full_real_df.head()


# In[11]:


full_shuffled_df = pd.concat(shuffled_models)
full_shuffled_df = full_shuffled_df.loc[:, output_cols].reset_index(drop=True)

output_real_file = os.path.join(output_dir, "repurposing_transformed_shuffled_models.tsv.gz")
full_shuffled_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(full_shuffled_df.shape)
full_shuffled_df.head()


# ## 5) Apply UMAP
# 
# ### Part 1: Apply UMAP to Cell Health Transformed Repurposing Hub Features

# In[12]:


cell_health_features = list(model_dict.keys())


# In[13]:


reducer = umap.UMAP(random_state=1234, n_components=2)

metadata_df = full_real_df.drop(cell_health_features, axis="columns")

real_embedding_df = pd.DataFrame(
    reducer.fit_transform(full_real_df.loc[:, cell_health_features]),
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

# In[14]:


cell_painting_features = [x for x in complete_df.columns if not x.startswith("Metadata_") ]


# In[15]:


reducer = umap.UMAP(random_state=1234, n_components=2)

complete_metadata_df = complete_df.drop(cell_painting_features, axis="columns")

complete_embedding_df = pd.DataFrame(
    reducer.fit_transform(complete_df.loc[:, cell_painting_features]),
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

