#!/usr/bin/env python
# coding: utf-8

# ## Apply all Cell-Health Models to Training and Testing Sets
# 
# **Gregory Way, 2019**

# In[1]:


import os
import pandas as pd
from joblib import load

from scripts.ml_utils import load_train_test


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


# Load all models
model_dir = "models"
model_string = "shuffle_False_transform_raw"

model_dict = {}
for model_file in os.listdir(model_dir):

    if model_string not in model_file:
        continue

    model_file_full = os.path.join(model_dir, model_file)
    cell_health_var = (
        model_file_full
        .split("/")[1]
        .replace("cell_health_target_", "")
        .replace("{}.joblib".format(model_string), "")
    )

    model_dict[cell_health_var] = load(model_file_full)


# In[4]:


# Load training and testing matrices
x_train_df, x_test_df, y_train_df, y_test_df = load_train_test(drop_metadata=True)


# In[5]:


# Apply the models
all_scores = []
for cell_health_feature in model_dict.keys():
    
    model_clf = model_dict[cell_health_feature]
    pred_train_df = (
        pd.DataFrame(model_clf.predict(x_train_df), columns=["score"])
        .assign(profiles=x_train_df.index,
            data_type="train",
                model=cell_health_feature)
    )
    pred_test_df = (
        pd.DataFrame(model_clf.predict(x_test_df), columns=["score"])
        .assign(profiles=x_test_df.index,
                data_type="test",
                model=cell_health_feature)
    )

    pred_df = pd.concat([pred_train_df, pred_test_df]).reset_index(drop=True)
    all_scores.append(pred_df)


# In[6]:


# Concatenate scores
all_scores = pd.concat(all_scores).reset_index(drop=True)

print(all_scores.shape)
all_scores.head(2)


# In[7]:


# Spread data into wide format
all_scores = (
    all_scores
    .pivot_table(index=["profiles", "data_type"],
                 columns="model",
                 values="score")
    .reset_index()
)

all_scores.head(2)


# In[8]:


# Load Metadata Mapping File
file = os.path.join("data", "profile_id_metadata_mapping.tsv")
metadata_df = pd.read_csv(file, sep='\t')

metadata_df.head()


# In[9]:


# Merge together
all_scores = (
    metadata_df.merge(all_scores,
                      left_on="Metadata_profile_id",
                      right_on="profiles")
    .drop("profiles", axis="columns")
)

all_scores.head()


# In[10]:


# Output file
file = os.path.join("results", "all_model_predictions.tsv")
all_scores.to_csv(file, sep='\t', index=False)

