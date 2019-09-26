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
import pandas as pd
from joblib import load

sys.path.append("../3.train")
from scripts.ml_utils import load_train_test, load_models


# ## 1) Load Models and Training Data

# In[2]:


model_dir = os.path.join("..", "3.train", "models")

model_dict, model_coef = load_models(model_dir=model_dir)
shuffle_model_dict, shuffle_model_coef = load_models(model_dir=model_dir, shuffle=True)


# In[3]:


data_dir = os.path.join("..", "3.train", "data")
x_train_df, x_test_df, y_train_df, y_test_df = load_train_test(data_dir=data_dir, drop_metadata=True)


# ## 2) Extract Repurposing Data Files

# In[4]:


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


# In[5]:


plate_info = {}
all_plates = os.listdir(repurposing_profile_dir)
for plate in all_plates:
    plate_dir = os.path.join(repurposing_profile_dir, plate, "n_all")
    norm_file = os.path.join(plate_dir, "{}_subsample_all_normalized.csv".format(plate))
    plate_info[plate] = norm_file


# ## 3) Apply all real and shuffled Models to all Repurposing Plates

# In[6]:


output_dir = os.path.join("data", "repurposing_transformed")


# In[7]:


real_models = []
shuffled_models = []
for plate in all_plates:
    norm_file = plate_info[plate]
    if os.path.exists(norm_file):
        df = pd.read_csv(norm_file)
        feature_df = df.reindex(x_test_df.columns, axis="columns").fillna(0)
        metadata_df = df.loc[:, df.columns.str.startswith("Metadata_")]
        
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

# In[8]:


full_real_df = pd.concat(real_models)

output_real_file = os.path.join(output_dir, "repurposing_transformed_real_models.tsv.gz")
full_real_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(full_real_df.shape)
full_real_df.head()


# In[9]:


shuffled_real_df = pd.concat(shuffled_models)

output_real_file = os.path.join(output_dir, "repurposing_transformed_shuffled_models.tsv.gz")
shuffled_real_df.to_csv(output_real_file, sep="\t", index=False, compression="gzip")

print(shuffled_real_df.shape)
shuffled_real_df.head()

